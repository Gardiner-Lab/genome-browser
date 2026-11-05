"""
The conversion toolchain and per-job outcome logic for the Genome Browser.

An uploaded BAM file or Wig_File is never stored as-is: a BAM holds every
individual sequencing read and can be tens to hundreds of gigabytes, so the
module converts it to a far smaller per-base coverage BigWig and keeps only
that (Requirements 13.1, 13.2, 13.9). This module is the one place that
conversion happens. The web app accepts the upload, writes the original into
``scratch/<job_id>/`` and commits a ``conversion_jobs`` row (task 6/8); the
single Conversion_Worker loop (task 9.2) claims a job and calls
:func:`convert_job` here. This module does the toolchain run, the before/after
database writes, and the scratch cleanup; the worker owns only the queue
polling and the claim transaction.

The design splits cleanly into two layers so the outcome logic is testable
without the real tools installed (the Testing Strategy stubs the toolchain in
property tests 20-26):

* **The toolchain layer** — :func:`_run_samtools_sort`, :func:`_run_samtools_index`,
  :func:`_run_bamcoverage`, :func:`_run_bedtools_genomecov`,
  :func:`_run_bedgraph_to_bigwig`, and :func:`_run_wig_to_bigwig`. Each is a
  thin wrapper around one external-tool command line, isolated so a test
  monkeypatches it (or passes a :class:`Toolchain` of stubs) and never shells
  out. The exact command lines each runs are documented on the function.

* **The outcome layer** — :func:`convert_job`. It validates the input header
  against the genome's ``chrom_sizes`` *before* any heavy work (Requirements
  13.7, 13.8), runs the toolchain for the kind, and then makes the terminal
  decision: on success it moves the BigWig into ``tracks/<id>/``, sets the
  track row ``track_status='ready'`` with its ``file_path``, and deletes the
  original (Requirements 13.4, 13.9); on failure it sets ``track_status='failed'``
  with a descriptive, reason-naming error and persists no BigWig (Requirement
  13.5). Either way it removes the whole job scratch directory in a ``finally``,
  leaving no orphans (Requirement 13.10).

The free-disk floor gate (:func:`has_free_disk_floor`) is a separate checkable
function the worker calls *before* claiming/starting a job: if free disk in the
scratch area is below ``CONVERT_MIN_FREE_GB`` the worker leaves the job queued
and waits, so the gate never starts heavy work it cannot finish (Requirement
13.13). Keeping it a pure predicate (with an injectable free-bytes reader) lets
a property test drive it with any reported free-disk value.

All database writes go through ``db.connect(data_root / 'genome.db')`` so the
app and the worker coordinate through the one WAL database, as everywhere else
in the module.
"""

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Optional

import db
import readers
from align import check_alignment

# One gigabyte in bytes, for the free-disk floor comparison.
_BYTES_PER_GB = 1024 ** 3


class ConversionError(Exception):
    """
    A conversion step failed with a member-facing reason.

    Raised by the toolchain wrappers and the validation step so
    :func:`convert_job` has a single, descriptive exception type to catch and
    record as the track's failure reason (Requirements 13.5, 13.8). The message
    is written verbatim into ``tracks.error``.
    """


# ---------------------------------------------------------------------------
# Free-disk floor gate (Requirement 13.13, Property 24)
# ---------------------------------------------------------------------------


def free_bytes(path) -> int:
    """
    Bytes of free space on the filesystem holding ``path``.

    Factored out (rather than inlined into the gate) so a test injects a fake
    free-bytes reader and drives the floor check with any reported value
    without touching a real disk. Uses ``shutil.disk_usage`` on the nearest
    existing ancestor of ``path`` so the check works before the per-job scratch
    directory itself has been created.
    """
    probe = Path(path)
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    return shutil.disk_usage(probe).free


def has_free_disk_floor(
    scratch_root,
    min_free_gb: int,
    free_bytes_reader: Callable[[object], int] = free_bytes,
) -> bool:
    """
    Whether free disk in the scratch area is at or above the safety floor.

    The worker calls this *before* it starts a queued job: ``True`` means there
    is at least ``min_free_gb`` free in ``scratch_root`` and the job may start;
    ``False`` means the worker must leave the job queued and wait, starting no
    heavy work (Requirement 13.13, Property 24). A floor of zero (or less) is
    always satisfied — the gate is effectively disabled — which keeps tests and
    tiny dev volumes from being blocked.

    ``free_bytes_reader`` defaults to :func:`free_bytes` and is injectable so a
    property test can report any free-disk value deterministically.
    """
    if min_free_gb <= 0:
        return True
    return free_bytes_reader(scratch_root) >= min_free_gb * _BYTES_PER_GB


# ---------------------------------------------------------------------------
# The toolchain layer — one wrapper per external command (Requirement 13.11)
# ---------------------------------------------------------------------------
#
# Each wrapper runs exactly one bundled tool and raises ConversionError on a
# non-zero exit, capturing the tool's stderr so the failure reason names the
# step that failed. They are gathered on a Toolchain dataclass so a caller (the
# worker, or a test) can swap the whole set for stubs in one place, while the
# default instance shells out to the real binaries. The command lines below are
# the exact invocations the worker performs.


def _run(cmd, step: str) -> None:
    """
    Run one toolchain command, raising :class:`ConversionError` on failure.

    Shared by the real wrappers: it runs ``cmd`` (a list, never a shell string,
    so member-supplied paths cannot inject), and on a non-zero exit or a missing
    binary raises a descriptive error naming the ``step`` and the tool's stderr,
    which becomes the track's failure reason (Requirement 13.8).
    """
    try:
        completed = subprocess.run(cmd, capture_output=True, text=True)
    except (OSError, ValueError) as exc:
        raise ConversionError(f'{step} could not be run: {exc}') from exc
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or '').strip()
        raise ConversionError(
            f'{step} failed (exit {completed.returncode}): {detail}'.strip())


def _run_samtools_sort(input_bam, sorted_bam, threads: int, mem_mb: int) -> None:
    """
    ``samtools sort -@ <threads> -m <mem_mb>M -o <sorted_bam> <input_bam>``

    Sort the uploaded BAM by coordinate, which ``samtools index`` and
    ``bamCoverage`` require. ``-@`` caps worker threads at ``SAMTOOLS_THREADS``
    (default 2) and ``-m`` bounds per-thread memory so a conversion cannot
    starve the other suite tools on the shared server (Requirement 13.15).
    """
    _run(['samtools', 'sort', '-@', str(threads), '-m', f'{mem_mb}M',
          '-o', str(sorted_bam), str(input_bam)], 'samtools sort')


def _run_samtools_index(sorted_bam) -> None:
    """
    ``samtools index <sorted_bam>``

    Build the ``.bai`` index that ``bamCoverage`` needs to read the sorted BAM.
    """
    _run(['samtools', 'index', str(sorted_bam)], 'samtools index')


def _run_bamcoverage(sorted_bam, out_bigwig, threads: int) -> None:
    """
    ``bamCoverage -b <sorted_bam> -o <out_bigwig> -of bigwig -p <threads>``

    deepTools' ``bamCoverage`` summarizes the sorted BAM into a per-base
    coverage BigWig directly — the preferred BAM path (Requirement 13.1). ``-p``
    caps its processes at the same thread budget as the sort step.
    """
    _run(['bamCoverage', '-b', str(sorted_bam), '-o', str(out_bigwig),
          '-of', 'bigwig', '-p', str(threads)], 'bamCoverage')


def _run_bedtools_genomecov(sorted_bam, chrom_sizes_path, out_bedgraph) -> None:
    """
    ``bedtools genomecov -bga -ibam <sorted_bam> -g <chrom_sizes>`` > bedGraph

    The documented fallback when deepTools' ``bamCoverage`` is unavailable:
    produce a per-base coverage bedGraph (``-bga`` reports zero-coverage
    regions too), which :func:`_run_bedgraph_to_bigwig` then packs into a
    BigWig. ``genomecov`` writes the bedGraph to stdout, so this wrapper
    redirects it to ``out_bedgraph``.
    """
    try:
        with open(out_bedgraph, 'w') as handle:
            completed = subprocess.run(
                ['bedtools', 'genomecov', '-bga', '-ibam', str(sorted_bam),
                 '-g', str(chrom_sizes_path)],
                stdout=handle, stderr=subprocess.PIPE, text=True)
    except (OSError, ValueError) as exc:
        raise ConversionError(
            f'bedtools genomecov could not be run: {exc}') from exc
    if completed.returncode != 0:
        raise ConversionError(
            f'bedtools genomecov failed (exit {completed.returncode}): '
            f'{(completed.stderr or "").strip()}'.strip())


def _run_bedgraph_to_bigwig(bedgraph, chrom_sizes_path, out_bigwig) -> None:
    """
    ``bedGraphToBigWig <bedgraph> <chrom_sizes> <out_bigwig>``

    The UCSC utility that packs the fallback bedGraph into a BigWig using the
    genome's chrom.sizes, completing the ``bedtools genomecov`` fallback path.
    """
    _run(['bedGraphToBigWig', str(bedgraph), str(chrom_sizes_path),
          str(out_bigwig)], 'bedGraphToBigWig')


def _run_wig_to_bigwig(input_wig, chrom_sizes_path, out_bigwig) -> None:
    """
    ``wigToBigWig <input_wig> <chrom_sizes> <out_bigwig>``

    The UCSC utility that converts an uploaded wig/wiggle (fixedStep,
    variableStep, or bedGraph) into a BigWig against the genome's chrom.sizes —
    the whole of the wig conversion path (Requirement 13.2).
    """
    _run(['wigToBigWig', str(input_wig), str(chrom_sizes_path),
          str(out_bigwig)], 'wigToBigWig')


@dataclass
class Toolchain:
    """
    The set of toolchain wrappers :func:`convert_job` calls.

    Bundling them here means a caller swaps the entire toolchain for a stub set
    in one place: the worker uses the default (real binaries), while a property
    test passes a :class:`Toolchain` whose members write a tiny valid BigWig on
    "success" or raise :class:`ConversionError` on "failure", so the outcome
    logic (status transitions, persist-nothing-on-failure, scratch cleanup)
    runs 100+ cheap iterations without the real tools installed.

    ``prefer_bamcoverage`` chooses the BAM coverage step: the deepTools
    ``bamCoverage`` path by default, or the ``bedtools genomecov`` →
    ``bedGraphToBigWig`` fallback when it is set false (or when ``bamCoverage``
    is not on ``PATH``, decided in :func:`convert_job`).
    """

    samtools_sort: Callable[..., None] = _run_samtools_sort
    samtools_index: Callable[..., None] = _run_samtools_index
    bamcoverage: Callable[..., None] = _run_bamcoverage
    bedtools_genomecov: Callable[..., None] = _run_bedtools_genomecov
    bedgraph_to_bigwig: Callable[..., None] = _run_bedgraph_to_bigwig
    wig_to_bigwig: Callable[..., None] = _run_wig_to_bigwig
    prefer_bamcoverage: bool = True


# ---------------------------------------------------------------------------
# The outcome layer — convert one job end to end (Requirements 13.4-13.10)
# ---------------------------------------------------------------------------


@dataclass
class ConversionOutcome:
    """
    The terminal result of one conversion job, for the worker to record on the
    ``conversion_jobs`` row.

    ``status`` is ``'ready'`` on success or ``'failed'`` on failure — the same
    value written to the track row. ``file_path`` is the stored BigWig's path on
    success (``None`` on failure, since no BigWig is persisted). ``error`` is the
    descriptive reason on failure (empty on success). The track row itself is
    already updated by :func:`convert_job`; this is the summary the worker uses
    to close out the job and surface status to the member (Requirement 13.19).
    """

    status: str
    file_path: Optional[str] = None
    error: str = ''

    @property
    def succeeded(self) -> bool:
        return self.status == 'ready'


def convert_job(
    *,
    track_id: str,
    genome_id: str,
    kind: str,
    input_path,
    data_root,
    samtools_threads: int,
    sort_mem_mb: int = 768,
    toolchain: Optional[Toolchain] = None,
) -> ConversionOutcome:
    """
    Convert one claimed job end to end, writing the track's terminal status.

    This is what the worker calls once it has claimed a ``queued`` job and the
    free-disk floor (:func:`has_free_disk_floor`) has passed. It does, in order:

    1. **Validate before heavy work.** Read the input header (BAM ``@SQ`` lines
       or wig ``chrom=`` declarations) and check it against the genome's
       ``chrom_sizes`` with :func:`align.check_alignment`. A chromosome name
       that does not match case-sensitively, or a coordinate past a matched
       chromosome's length, fails the job *before* any sort/coverage work, with
       an error naming the offenders and no BigWig persisted (Requirements 13.7,
       13.8).
    2. **Run the toolchain for the kind.** BAM: ``samtools sort`` → ``samtools
       index`` → ``bamCoverage`` (or the ``bedtools genomecov`` →
       ``bedGraphToBigWig`` fallback). wig: ``wigToBigWig`` against a chrom.sizes
       file written from the ``chrom_sizes`` table (Requirements 13.1, 13.2,
       13.15). All intermediates live under the job scratch directory.
    3. **Record the terminal outcome.** On success, move the BigWig into
       ``tracks/<track_id>/``, set the track ``track_status='ready'`` with its
       ``file_path``, and delete the original input (Requirements 13.4, 13.9).
       On any :class:`ConversionError`, set ``track_status='failed'`` with the
       reason and persist no BigWig (Requirement 13.5).
    4. **Always clean scratch.** The whole ``scratch/<job_id>/`` directory is
       removed in a ``finally``, on both success and failure, leaving no
       orphaned conversion files (Requirement 13.10).

    ``input_path`` is the uploaded original inside ``scratch/<job_id>/``; its
    parent directory *is* that job scratch directory, which is what gets cleaned.
    ``sort_mem_mb`` bounds ``samtools sort`` per-thread memory (the ``-m``
    budget). All DB writes go through ``db.connect(data_root / 'genome.db')``.
    Returns a :class:`ConversionOutcome` for the worker to close out the job.
    """
    tools = toolchain or Toolchain()
    data_root = Path(data_root)
    input_path = Path(input_path)
    scratch_dir = input_path.parent
    db_path = data_root / 'genome.db'

    try:
        try:
            # 1. Validate the input header against the reference first.
            chrom_sizes = _load_chrom_sizes(db_path, genome_id)
            _validate_input(kind, input_path, chrom_sizes)

            # 2. Run the toolchain, producing a BigWig in scratch.
            scratch_bigwig = scratch_dir / f'{track_id}.bw'
            if kind == 'bam':
                _convert_bam(tools, input_path, scratch_bigwig, scratch_dir,
                             chrom_sizes, samtools_threads, sort_mem_mb)
            elif kind == 'wig':
                _convert_wig(tools, input_path, scratch_bigwig, scratch_dir,
                             chrom_sizes)
            else:
                raise ConversionError(
                    f'Unknown conversion kind {kind!r}; expected "bam" or "wig".')

            if not scratch_bigwig.exists():
                raise ConversionError(
                    'Conversion produced no BigWig output file.')

            # 3a. Success: move the BigWig into the track's store, mark ready,
            #     delete the original input.
            final_path = _store_bigwig(data_root, track_id, scratch_bigwig)
            _mark_ready(db_path, track_id, str(final_path))
            _delete_original(input_path)
            return ConversionOutcome(status='ready', file_path=str(final_path))

        except ConversionError as exc:
            # 3b. Failure: mark failed with the reason, persist no BigWig.
            _mark_failed(db_path, track_id, str(exc))
            return ConversionOutcome(status='failed', error=str(exc))
        except Exception as exc:  # pragma: no cover - defensive catch-all
            # Any unexpected error is still a conversion failure: record a
            # descriptive reason and persist nothing (Requirement 13.5).
            message = f'Conversion failed: {exc}'
            _mark_failed(db_path, track_id, message)
            return ConversionOutcome(status='failed', error=message)
    finally:
        # 4. Always remove the job scratch directory, success or failure,
        #    leaving no orphaned conversion files (Requirement 13.10).
        _remove_scratch(scratch_dir)


# ---------------------------------------------------------------------------
# Step helpers
# ---------------------------------------------------------------------------


def _validate_input(kind: str, input_path, chrom_sizes: Dict[str, int]) -> None:
    """
    Check the input header against the reference before heavy work.

    Reads the chromosome names/extents the input references (BAM ``@SQ`` lines
    or wig declarations) and runs the shared alignment rule
    (:func:`align.check_alignment`) against the genome's ``chrom_sizes``. On a
    rejection it raises :class:`ConversionError` with the alignment result's
    full, offender-naming message, so the job fails before any sort/coverage
    work and persists nothing (Requirements 13.7, 13.8).
    """
    fmt = 'bam' if kind == 'bam' else 'wig' if kind == 'wig' else kind
    header = readers.read_track_header(fmt, str(input_path))
    result = check_alignment(header, chrom_sizes)
    if not result.accepted:
        raise ConversionError(result.error_message
                              or 'Conversion input does not align to the '
                                 'reference genome.')


def _convert_bam(tools: Toolchain, input_bam, out_bigwig, scratch_dir,
                 chrom_sizes: Dict[str, int], threads: int,
                 sort_mem_mb: int) -> None:
    """
    BAM path: sort → index → coverage BigWig (Requirements 13.1, 13.15).

    ``samtools sort`` (threads capped, ``-m`` memory budget) then ``samtools
    index`` produce the sorted/indexed BAM the coverage step reads. deepTools'
    ``bamCoverage`` is preferred; when it is unavailable (not on ``PATH``, or the
    toolchain is configured to skip it) the documented fallback runs ``bedtools
    genomecov`` → ``bedGraphToBigWig`` against a chrom.sizes written from the
    genome. All intermediates live under ``scratch_dir`` and are cleaned with it.
    """
    sorted_bam = scratch_dir / 'sorted.bam'
    tools.samtools_sort(input_bam, sorted_bam, threads, sort_mem_mb)
    tools.samtools_index(sorted_bam)

    if tools.prefer_bamcoverage and shutil.which('bamCoverage') is not None:
        tools.bamcoverage(sorted_bam, out_bigwig, threads)
        return

    # Fallback: bedtools genomecov -> bedGraphToBigWig, which needs chrom.sizes.
    chrom_sizes_path = scratch_dir / 'chrom.sizes'
    _write_chrom_sizes_file(chrom_sizes, chrom_sizes_path)
    bedgraph = scratch_dir / 'coverage.bedgraph'
    tools.bedtools_genomecov(sorted_bam, chrom_sizes_path, bedgraph)
    tools.bedgraph_to_bigwig(bedgraph, chrom_sizes_path, out_bigwig)


def _convert_wig(tools: Toolchain, input_wig, out_bigwig, scratch_dir,
                 chrom_sizes: Dict[str, int]) -> None:
    """
    wig path: ``wigToBigWig`` against the genome's chrom.sizes (Req 13.2).

    Writes a chrom.sizes file from the ``chrom_sizes`` table (the UCSC utility
    needs the per-chromosome lengths to bound the BigWig) and runs
    ``wigToBigWig`` to produce the BigWig.
    """
    chrom_sizes_path = scratch_dir / 'chrom.sizes'
    _write_chrom_sizes_file(chrom_sizes, chrom_sizes_path)
    tools.wig_to_bigwig(input_wig, chrom_sizes_path, out_bigwig)


def _store_bigwig(data_root, track_id: str, scratch_bigwig) -> Path:
    """
    Move the produced BigWig into ``tracks/<track_id>/`` and return its path.

    The stored BigWig is the only artifact the module keeps from a conversion
    (Requirement 13.4). Moving (not copying) off the scratch volume means the
    scratch cleanup afterwards leaves only the stored file behind.
    """
    track_dir = Path(data_root) / 'tracks' / track_id
    track_dir.mkdir(parents=True, exist_ok=True)
    final_path = track_dir / f'{track_id}.bw'
    shutil.move(str(scratch_bigwig), str(final_path))
    return final_path


def _delete_original(input_path) -> None:
    """
    Delete the uploaded original BAM/wig on success (Requirement 13.9).

    The module never retains the original or makes it re-downloadable. The whole
    scratch directory is removed in the ``finally`` regardless; this makes the
    intent explicit and keeps success semantics clear.
    """
    try:
        Path(input_path).unlink()
    except FileNotFoundError:
        pass


def _remove_scratch(scratch_dir) -> None:
    """
    Remove the entire job scratch directory, leaving no orphans (Req 13.10).

    Called in :func:`convert_job`'s ``finally`` on both success and failure.
    ``ignore_errors`` keeps a half-written or already-removed directory from
    masking the job's real outcome.
    """
    shutil.rmtree(Path(scratch_dir), ignore_errors=True)


def _write_chrom_sizes_file(chrom_sizes: Dict[str, int], path) -> None:
    """
    Write a UCSC chrom.sizes file (``<chrom>\\t<length>`` per line).

    ``wigToBigWig`` and ``bedGraphToBigWig`` take the genome's chromosome
    lengths as a file; this materializes the ``chrom_sizes`` table rows into one
    inside the job scratch directory.
    """
    with open(path, 'w', encoding='utf-8') as handle:
        for chrom, length in chrom_sizes.items():
            handle.write(f'{chrom}\t{length}\n')


# ---------------------------------------------------------------------------
# Database access — all through db.connect(data_root / 'genome.db')
# ---------------------------------------------------------------------------


def _load_chrom_sizes(db_path, genome_id: str) -> Dict[str, int]:
    """
    The genome's ``chrom_sizes`` as a name -> length map for validation.

    Raises :class:`ConversionError` when the genome has no chromosome rows (a
    genome that cannot be validated against is a conversion failure, not a
    silent pass).
    """
    conn = db.connect(db_path)
    try:
        rows = conn.execute(
            'SELECT chrom, length_bp FROM chrom_sizes WHERE genome_id = ?',
            (genome_id,)).fetchall()
    finally:
        conn.close()
    if not rows:
        raise ConversionError(
            f'Reference genome {genome_id} has no chromosome sizes to '
            f'validate against.')
    return {row['chrom']: int(row['length_bp']) for row in rows}


def _mark_ready(db_path, track_id: str, file_path: str) -> None:
    """
    Set the track ``track_status='ready'`` with its stored BigWig path (13.4).

    Clears any prior error. The track becomes viewable in the Track_Library only
    once this runs (Requirement 13.6, enforced by the library query).
    """
    conn = db.connect(db_path)
    try:
        conn.execute(
            "UPDATE tracks SET track_status = 'ready', file_path = ?, "
            "error = NULL WHERE id = ?",
            (file_path, track_id))
        conn.commit()
    finally:
        conn.close()


def _mark_failed(db_path, track_id: str, error: str) -> None:
    """
    Set the track ``track_status='failed'`` with a descriptive error (13.5).

    Persists no BigWig: ``file_path`` is cleared so a failed conversion can
    never leave a stale path pointing at a partial or absent file. The returning
    member sees ``failed`` and this message (Requirement 13.19).
    """
    conn = db.connect(db_path)
    try:
        conn.execute(
            "UPDATE tracks SET track_status = 'failed', error = ?, "
            "file_path = NULL WHERE id = ?",
            (error, track_id))
        conn.commit()
    finally:
        conn.close()
