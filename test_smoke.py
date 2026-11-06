"""
Build/smoke check for the Genome Browser image and its conversion toolchain
(task 18.1).

This test has two groups with deliberately different runtime requirements:

* **Group A — static assertions about the image.** These read the
  ``Dockerfile`` and inspect ``conversion.py``'s command construction. They need
  neither Docker nor the conversion toolchain, so they ALWAYS run and must pass.
  They pin the settings a future maintainer could silently break: the bundled
  toolchain and gunicorn entrypoint (Req 13.11), and the samtools-sort thread cap
  and ``-m`` budget (Req 13.15).

* **Group B — toolchain/convert smoke.** These exercise the real conversion
  binaries (``samtools``/``bamCoverage``/``wigToBigWig``/``bedGraphToBigWig``).
  The binaries live in the built image, not on the dev/CI host, so each Group B
  test SKIPS with a clear reason when its tools (or Docker) are unavailable, and
  runs the real conversion only where they are present. A full in-image build is
  additionally blocked here by the Dockerfile's intentionally-unpinned UCSC
  checksum placeholder, so the Docker-build test skips until that is pinned.

How this tool is wired into the suite — the worker's CPU and memory caps, its
concurrency of one, the published image both containers share, and the nginx
block that streams large uploads — is asserted in the suite's own repository
(``tests/test_genome_wiring.py`` in Trellis Tools), where ``docker-compose.yml``
and ``nginx/default.conf`` live.

Run it as ``GENOME_DATA_ROOT="$(mktemp -d)" python -m pytest test_smoke.py``. In
an environment without the toolchain or a running Docker daemon, Group A passes
and Group B skips.
"""

import re
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest

import conversion

GENOME_DIR = Path(__file__).resolve().parent
DOCKERFILE_PATH = GENOME_DIR / 'Dockerfile'

# The conversion tools that must be on PATH inside the built image.
CONVERSION_TOOLS = ('samtools', 'bamCoverage', 'wigToBigWig', 'bedGraphToBigWig')


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def dockerfile():
    """The raw text of ``genome/Dockerfile``."""
    assert DOCKERFILE_PATH.is_file(), f'missing {DOCKERFILE_PATH}'
    return DOCKERFILE_PATH.read_text(encoding='utf-8')


# ===========================================================================
# GROUP A - static assertions about the image (ALWAYS run; no Docker/toolchain)
# ===========================================================================


def test_dockerfile_installs_toolchain_and_gunicorn(dockerfile):
    """
    The image installs the toolchain, fetches UCSC utils, exposes 5008, and the
    default CMD is gunicorn (Req 13.11).
    """
    # apt-installed toolchain packages.
    for pkg in ('samtools', 'bedtools', 'deeptools', 'python3-pybigwig'):
        assert pkg in dockerfile, f'Dockerfile must install {pkg}'

    # UCSC static binaries fetched at build time.
    assert 'wigToBigWig' in dockerfile, 'Dockerfile must fetch wigToBigWig'
    assert 'bedGraphToBigWig' in dockerfile, \
        'Dockerfile must fetch bedGraphToBigWig'

    # Port and entrypoint.
    assert re.search(r'EXPOSE\s+5008', dockerfile), \
        'Dockerfile must EXPOSE 5008'
    cmd_match = re.search(r'^CMD\s+(.+)$', dockerfile, flags=re.MULTILINE)
    assert cmd_match, 'Dockerfile must define a default CMD'
    assert 'gunicorn' in cmd_match.group(1), \
        f'default CMD should be gunicorn, got {cmd_match.group(1)!r}'


def test_samtools_sort_caps_threads_and_memory():
    """
    samtools sort uses the configured thread cap (-@ <threads>) and -m budget
    (Req 13.15). Inspects conversion.py's command construction directly, so this
    needs no tools installed.
    """
    captured = {}

    def fake_run(cmd, step):
        captured['cmd'] = cmd
        captured['step'] = step

    original = conversion._run
    conversion._run = fake_run
    try:
        conversion._run_samtools_sort(
            'in.bam', 'sorted.bam', threads=2, mem_mb=768)
    finally:
        conversion._run = original

    cmd = captured['cmd']
    assert cmd[0:2] == ['samtools', 'sort'], \
        f'expected a samtools sort invocation, got {cmd!r}'

    # Thread cap: -@ immediately followed by the configured thread count.
    assert '-@' in cmd, f'samtools sort must pass -@, got {cmd!r}'
    assert cmd[cmd.index('-@') + 1] == '2', \
        f'samtools sort -@ must carry the thread cap, got {cmd!r}'

    # Memory budget: -m immediately followed by the per-thread budget in M.
    assert '-m' in cmd, f'samtools sort must pass -m, got {cmd!r}'
    assert cmd[cmd.index('-m') + 1] == '768M', \
        f'samtools sort -m must carry the memory budget, got {cmd!r}'


# ===========================================================================
# GROUP B — toolchain/convert smoke (SKIPPED when the tools/Docker are absent)
# ===========================================================================


def _missing_tools():
    """Conversion tools not found on PATH."""
    return [tool for tool in CONVERSION_TOOLS if shutil.which(tool) is None]


def _assert_valid_bigwig(path):
    """
    Assert ``path`` is a valid BigWig.

    Prefers pyBigWig (also bundled in the image) when importable; otherwise
    falls back to checking the BigWig magic number in the file header, so the
    check still means something where pyBigWig is not installed.
    """
    try:
        import pyBigWig  # noqa: WPS433 (optional dependency, bundled in image)
    except ImportError:
        header = Path(path).read_bytes()[:4]
        # BigWig magic 0x888FFC26, little- or big-endian on disk.
        magic = {struct.pack('<I', 0x888FFC26), struct.pack('>I', 0x888FFC26)}
        assert header in magic, \
            f'{path} is not a valid BigWig (bad magic {header!r})'
        return

    handle = pyBigWig.open(str(path))
    try:
        assert handle.isBigWig(), f'{path} is not a valid BigWig'
    finally:
        handle.close()


def test_conversion_tools_on_path():
    """
    samtools/bamCoverage/wigToBigWig/bedGraphToBigWig are on PATH.

    In the built image all four are present; on a bare host they are not, so
    this skips with a reason rather than failing. This is the assertion that
    the image ships the toolchain (Req 13.11).
    """
    missing = _missing_tools()
    if missing:
        pytest.skip(
            'conversion toolchain not installed on host '
            f'(missing {", ".join(missing)}); runs in the built image')
    # All present: nothing more to assert here beyond their presence.
    assert not missing


def test_wig_converts_to_valid_bigwig(tmp_path):
    """
    A tiny real wig converts to a valid BigWig via wigToBigWig.

    Skipped unless wigToBigWig is on PATH (the built image). Builds a minimal
    chrom.sizes and a one-track fixedStep wig, runs the real wigToBigWig, and
    checks the output is a valid BigWig.
    """
    if shutil.which('wigToBigWig') is None:
        pytest.skip(
            'wigToBigWig not installed on host; runs in the built image')

    chrom_sizes = tmp_path / 'chrom.sizes'
    chrom_sizes.write_text('chr1\t1000\n', encoding='utf-8')

    wig = tmp_path / 'tiny.wig'
    wig.write_text(
        'fixedStep chrom=chr1 start=1 step=1 span=1\n'
        + '\n'.join(str(v) for v in (1, 2, 3, 4, 5)) + '\n',
        encoding='utf-8')

    out = tmp_path / 'tiny.bw'
    result = subprocess.run(
        ['wigToBigWig', str(wig), str(chrom_sizes), str(out)],
        capture_output=True, text=True)
    assert result.returncode == 0, \
        f'wigToBigWig failed: {result.stderr or result.stdout}'
    assert out.is_file(), 'wigToBigWig produced no output file'
    _assert_valid_bigwig(out)


def test_bam_converts_to_valid_bigwig(tmp_path):
    """
    A tiny real BAM converts to a valid BigWig through conversion.convert_job.

    Skipped unless samtools AND bamCoverage are on PATH (the built image).
    Builds a tiny SAM, encodes it to BAM with samtools, seeds a temp genome.db
    with the matching chrom size, and runs the REAL toolchain end to end via
    conversion.convert_job, asserting a ready outcome and a valid BigWig.
    """
    if shutil.which('samtools') is None or shutil.which('bamCoverage') is None:
        pytest.skip(
            'samtools/bamCoverage not installed on host; '
            'runs in the built image')

    import db  # local import: only needed on the toolchain path.

    # 1. A tiny coordinate-ordered SAM with one chromosome and two reads.
    sam = tmp_path / 'tiny.sam'
    sam.write_text(
        '@HD\tVN:1.6\tSO:coordinate\n'
        '@SQ\tSN:chr1\tLN:1000\n'
        'r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n'
        'r2\t0\tchr1\t20\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n',
        encoding='utf-8')
    bam = tmp_path / 'input.bam'
    enc = subprocess.run(
        ['samtools', 'view', '-b', '-o', str(bam), str(sam)],
        capture_output=True, text=True)
    assert enc.returncode == 0, \
        f'samtools view failed to build the test BAM: {enc.stderr}'

    # 2. A genome.db with the matching chromosome so validation passes.
    data_root = tmp_path / 'data'
    data_root.mkdir()
    db_path = data_root / 'genome.db'
    db.init_db(db_path)
    conn = db.connect(db_path)
    conn.execute(
        'INSERT INTO reference_genomes '
        '(id, name, assembly_id, owner_lab, genome_visibility, created_by) '
        "VALUES ('g1', 'G', 'asm', 'labA', 'owner_lab_only', 'alice')")
    conn.execute(
        'INSERT INTO chrom_sizes (genome_id, chrom, length_bp) '
        "VALUES ('g1', 'chr1', 1000)")
    conn.execute(
        'INSERT INTO tracks (id, genome_id, name, type, owner_lab, '
        'track_status, created_by) '
        "VALUES ('t1', 'g1', 'Track', 'bigwig', 'labA', 'converting', 'alice')")
    conn.commit()
    conn.close()

    # 3. The upload lives inside its own job scratch dir, as in production.
    scratch = data_root / 'scratch' / 'job1'
    scratch.mkdir(parents=True)
    staged = scratch / 'input.bam'
    shutil.move(str(bam), str(staged))

    # 4. Run the REAL toolchain end to end.
    outcome = conversion.convert_job(
        track_id='t1', genome_id='g1', kind='bam',
        input_path=staged, data_root=data_root, samtools_threads=2)

    assert outcome.succeeded, f'conversion failed: {outcome.error}'
    assert outcome.file_path and Path(outcome.file_path).is_file(), \
        'conversion reported ready but stored no BigWig'
    _assert_valid_bigwig(outcome.file_path)


def test_docker_build_smoke():
    """
    A full `docker build` of the image.

    Skipped when Docker is unavailable, the daemon is down, OR the Dockerfile's
    UCSC checksum placeholder is still unpinned (which makes a real build fail
    verification by design). A maintainer who pins the checksums and runs this
    where a daemon is up gets a real end-to-end build check.
    """
    if shutil.which('docker') is None:
        pytest.skip('docker CLI not installed; build smoke runs in CI/release')

    ping = subprocess.run(
        ['docker', 'info'], capture_output=True, text=True)
    if ping.returncode != 0:
        pytest.skip('docker daemon not running; build smoke runs in CI/release')

    dockerfile_text = DOCKERFILE_PATH.read_text(encoding='utf-8')
    if re.search(r'ARG\s+UCSC_WIGTOBIGWIG_SHA256=\s*$',
                 dockerfile_text, flags=re.MULTILINE):
        pytest.skip(
            'Dockerfile UCSC checksum placeholder is unpinned; a real build '
            'fails verification by design until a maintainer pins the SHA-256s')

    build = subprocess.run(
        ['docker', 'build', '-t', 'trellis-genome:smoke', '.'],
        cwd=str(GENOME_DIR), capture_output=True, text=True)
    assert build.returncode == 0, \
        f'docker build failed: {build.stderr[-2000:]}'

    # The four conversion tools must be on PATH inside the built image.
    for tool in CONVERSION_TOOLS:
        which = subprocess.run(
            ['docker', 'run', '--rm', 'trellis-genome:smoke',
             'which', tool],
            capture_output=True, text=True)
        assert which.returncode == 0, \
            f'{tool} not on PATH in the built image'
