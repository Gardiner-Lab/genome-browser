# Changelog

Releases of the Genome Browser. The [Trellis Tools](https://github.com/Gardiner-Lab/Trellis-Tools)
suite pins a version of this tool as `GENOME_VERSION`, so a release here changes no installation
until the suite is updated to use it.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [semantic versioning](https://semver.org/).

## [Unreleased]

## [0.1.0] - 2026-11-06

The Genome Browser's first release of its own. It shipped inside the suite repository as `genome/`
in Trellis Tools 2.33.0 to 2.34.4, pinned to the suite's own version; it now lives here, with its
history, versioned and released like the suite's other tools. Nothing about how it runs changes:
the same two containers, the same `/genome/` route, the same `genome-data` volume. Only the image
name changes, from `ghcr.io/gardiner-lab/lab-suite-genome` to
`ghcr.io/gardiner-lab/genome-browser`.

### Added
- Reference genomes uploaded as FASTA (up to 4 GB), and a library of aligned tracks over them:
  BigWig, BED and GFF, with BAM and wig files converted to BigWig on the server so only the
  coverage summary is kept. A track is accepted only when every chromosome name and coordinate
  matches its reference genome exactly.
- An IGV-style browser, navigated by locus or region, with the current view exportable as an SVG
  figure or sent straight to a Lab Notebook page.
- Tags linking a track or genome to a plasmid (MoClo or the Sequence Editor) or a seed stock; a
  line or plant tag feeds that line's family tree in the Seed Bank (`CONTRACT_seed_stock_tags.md`).
- Sharing on the suite's model: a reference genome defaults to its owner lab's chair group, a track
  to its owner lab, and either can be shared further by a member of the owner lab.
- Conversion in a single worker with its own CPU and memory limits, coordinating through a durable
  on-disk job queue, so an accepted conversion survives a restart.
- The conversion toolchain baked into the image per architecture (samtools, bedtools, deepTools and
  the UCSC utilities), with the UCSC binaries verified against pinned checksums at build time.

### Fixed
Carried over from the suite releases that fixed this tool:
- The containers would not start: the image provides neither `gunicorn` nor `python` by those
  names (suite 2.34.1).
- Every request answered 404 through the proxy, which passed the `/genome` prefix on to an app
  serving at its root (suite 2.34.2).
- The worker could crash on a fresh install, querying the job table before the web container had
  created the schema; it now creates its directories and schema itself before recovering. Its
  health probe, inherited from a web image it does not serve web requests from, is disabled
  (suite 2.34.2).
- The health endpoint reported version `0.1.0` while the image was published as 2.34.2
  (suite 2.34.3).
- `genome.db` was in no backup at all: its volume reached neither the suite's nightly backup nor
  its restore service (suite 2.34.4). The genome and track files are still deliberately left out,
  being far too large for a nightly archive; see the suite's `docs/BACKUP.md`.
