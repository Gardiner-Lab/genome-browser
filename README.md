# Genome Browser

Trellis Tools' reference genomes and the aligned data tracks over them, viewed in an IGV-style
browser.

- **Reference genomes** are uploaded as FASTA, up to 4 GB. Each belongs to a lab and is shared with
  that lab's chair group by default.
- **Tracks** are BigWig, BED or GFF. BAM and wig files are converted to BigWig on the server, so
  only the coverage summary is kept rather than the alignment.
  - A track is accepted only when **every chromosome name and coordinate matches its reference
    genome exactly** — a mismatch is refused with the names that disagree, rather than drawn in the
    wrong place.
  - Conversion runs in a single worker with its own CPU and memory limits, pulling from a durable
    on-disk queue, so an accepted job survives a restart and never ties up the web service.
- **Navigation** is by locus or region, and a view can be exported as an **SVG figure** or sent
  straight to a **Lab Notebook** page.
- **Tags** link a track or genome to a plasmid (from MoClo or the Sequence Editor) or to a seed
  stock. A line or plant tag feeds that line's family tree in the Seed Bank; the shape of that data
  is set out in [`CONTRACT_seed_stock_tags.md`](CONTRACT_seed_stock_tags.md).
- **Who sees what** follows the suite: a genome defaults to its owner lab's chair group, a track to
  its owner lab, and either can be shared further by a member of the owner lab. Administrators see
  all.

## Running it

It runs inside [Trellis Tools](https://github.com/Gardiner-Lab/Trellis-Tools), behind its proxy at
`/genome/`, signed in through the suite. It has no accounts of its own: the proxy says who is signed
in in an `X-Suite-User` header. **Don't expose it outside the proxy** — reachable directly, anyone
could set that header themselves.

The suite runs **two containers from this one image**: the web service (the image's own command) and
the conversion worker (`python3 worker.py`). The worker is the only process that runs the heavy
toolchain, and there is exactly one of it.

| Setting | Default | |
| --- | --- | --- |
| `GENOME_DATA_ROOT` | `/genome` | its own volume: the database, the genome and track files, the conversion scratch |
| `SUITE_DATABASE_PATH` | `/suite/suite.db` | labs and chair groups, read-only |
| `MOCLO_DB` | `/moclo/moclo.db` | the suite's account directory, read-only |
| `URL_PREFIX` | `/genome` | where the proxy mounts it; the prefix is stripped before requests arrive |
| `SUITE_AUTH` | — | `proxy` to trust the suite's header; the only supported mode inside the suite |
| `NOTEBOOK_URL`, `MOCLO_URL`, `EDITOR_URL`, `SEEDS_URL` | `/notebook`, `/`, `/editor`, `/seeds` | the other tools, for tags, links and notebook export |
| `MAX_BAM_GB` | `50` | the largest BAM an upload may be, checked before any job starts |
| `CONVERT_MIN_FREE_GB` | `50` | a conversion only starts with at least this much free disk |
| `SAMTOOLS_THREADS` | `2` | how many threads `samtools sort` may use |

Inside the container the health route is `/health`, **not** `/genome/health`: the proxy strips the
prefix. Through the suite it is at `/genome/health`.

## The viewer, and citing it

The interactive browser is [igv.js](https://github.com/igvteam/igv.js) 3.8.10, vendored unchanged
into `static/igv.min.js` rather than loaded from a CDN, so the suite has no outbound dependency at
run time. It is MIT licensed; the notice travels with it in `static/igv.js-LICENSE.txt`, and
`NOTICE` records it along with the libraries its bundle carries.

**If a view from here goes into a paper, cite igv.js:** James T Robinson, Helga Thorvaldsdóttir,
Douglass Turner, Jill P Mesirov, "igv.js: an embeddable JavaScript implementation of the Integrative
Genomics Viewer (IGV)", *Bioinformatics*, [doi:10.1093/bioinformatics/btac830][igvjs].

[igvjs]: https://doi.org/10.1093/bioinformatics/btac830

## The conversion toolchain

The image is Debian-based rather than Alpine because the toolchain is not reliably packaged for
Alpine: `samtools`, `bedtools` and deepTools (`bamCoverage`) come from Debian packages, and the two
UCSC utilities (`wigToBigWig`, `bedGraphToBigWig`) are fetched per architecture at build time and
**verified against pinned SHA-256 checksums**, so the build fails rather than silently accepting a
substituted binary. UCSC rebuilds those binaries periodically; if a build fails the checksum, see
the note in the `Dockerfile` before changing anything.

## Tests

```sh
GENOME_DATA_ROOT="$(mktemp -d)" pytest -q
```

`GENOME_DATA_ROOT` matters: the application creates its data root when the module is imported, and
the default is `/genome`, which an ordinary user cannot make.

The conversion tests that need the real binaries skip where the toolchain is absent, so a plain
checkout runs everything else. To run them too, run the suite inside this image:

```sh
docker run --rm --user 0 -v "$PWD":/src:ro ghcr.io/gardiner-lab/genome-browser:0.1.0 sh -c \
  'mkdir -p /work && cp -r /src/. /work/ && apt-get update -qq \
   && apt-get install -y -qq python3-pytest python3-yaml \
   && cd /work && GENOME_DATA_ROOT=/tmp/g python3 -m pytest -q'
```

## Releases

The image is `ghcr.io/gardiner-lab/genome-browser`, built for Intel and ARM servers. Pushing a
`v0.1.0` tag publishes it; the tag must match `APP_VERSION` in `app.py` or the release stops. The
suite pins an exact version as `GENOME_VERSION`, so publishing here changes no installation until
the suite is updated to use it.

This tool shipped inside the suite repository as `genome/` for its first releases (Trellis Tools
2.33.0 to 2.34.4) and moved here, with its history, at 0.1.0. See `CHANGELOG.md`.
