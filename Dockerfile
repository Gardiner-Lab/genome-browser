# Genome Browser module image.
#
# Debian base (not Alpine): the conversion toolchain — samtools, deepTools
# (bamCoverage), bedtools, and the UCSC wigToBigWig/bedGraphToBigWig utilities —
# is not reliably available via apk. We install the toolchain from Debian
# packages plus the two static UCSC binaries so conversions run entirely inside
# this image rather than against host tools. (Requirement 13.11.)
#
# This same image runs BOTH containers defined in docker-compose:
#   * genome-browser (web)  -> default CMD below (gunicorn)
#   * genome-worker         -> compose overrides the command with
#                              `python3 worker.py` (single Conversion_Worker,
#                              concurrency 1). No separate Dockerfile is needed.
# CPU/memory limits for the worker are applied by compose (Requirement 13.16).
FROM debian:bookworm-slim

# ---------------------------------------------------------------------------
# The genome the browser ships with (bundled.py): TAIR10 with the Araport11
# annotation, from TAIR (CC BY 4.0; see NOTICE). Kept as assets of this repo's
# data-araport11 release rather than in git, and checked against the SHA-256s
# pinned here, so a download cannot be silently substituted. First, so the
# 50 MB is cached apart from the packages below and not fetched again when they
# change. The worker adds them as the suite's genome on a volume's first start.
ADD --checksum=sha256:d43e19a85e6ae5dd18cf6adea2d21ed15ad8a58ae374f58dbfdbfcff1b9b13eb --chmod=644     https://github.com/Gardiner-Lab/genome-browser/releases/download/data-araport11/TAIR10_chr_all.fasta.gz /opt/bundled/TAIR10_chr_all.fasta.gz
ADD --checksum=sha256:cedd155953964ad06b492741569d4e5a61d8db7a1e52a0def164f9f32169d03d --chmod=644     https://github.com/Gardiner-Lab/genome-browser/releases/download/data-araport11/Araport11_GFF3_genes_transposons.20250813.gff.gz /opt/bundled/Araport11_GFF3_genes_transposons.20250813.gff.gz
# ADD --chmod gives the directory it creates the same 644, which the unprivileged
# worker cannot enter; the files stay 644.
RUN chmod 755 /opt/bundled

# ---------------------------------------------------------------------------
# System packages + bundled conversion toolchain.
#
# The suite prefers Debian packages over pip, so Flask/gunicorn/requests and the
# bio toolchain all come from apt. python3-requests is required because the
# module imports `requests` (e.g. notebook.py export to the Lab Notebook).
# tabix brings bgzip and tabix, which index BED/GFF tracks for the viewer
# (serving.py); samtools faidx indexes the reference sequences. ncbi-blast+
# checks primers against a genome (primer_check.py): blastn, and makeblastdb for
# the database the worker makes of a genome the first time it is checked.
#
# Supply chain: the UCSC static binaries are verified against pinned SHA-256s.
# They are fetched at build time from the canonical UCSC download host, per
# architecture (linux.x86_64 for amd64, linux.aarch64 for arm64), and each is
# checked against the SHA-256 pinned below (captured 2026-11-05). The
# ``sha256sum -c`` check fails the build on any mismatch, so a download cannot
# be silently substituted. UCSC periodically rebuilds these binaries; if a
# future build fails the check, re-fetch each for its architecture, confirm its
# provenance, run ``sha256sum <file>``, and update the matching ARG below.
#
# TARGETARCH is provided by BuildKit (amd64 / arm64); the per-arch URL and
# checksums are selected from it at build time.
ARG TARGETARCH
ARG UCSC_WIGTOBIGWIG_SHA256_amd64=5b9eccc4929b339be26bb69bce27de63070195cb0a431b410627c0bed8ab4f15
ARG UCSC_BEDGRAPHTOBIGWIG_SHA256_amd64=fe9ed8ed858731d22d9a3cefa99d192f51bff1496693599b2f8054fa350adb48
ARG UCSC_WIGTOBIGWIG_SHA256_arm64=fb347a6b90a7d4cbae4ac035a3c35c1700589e1a69f9553468d29754efbdacbf
ARG UCSC_BEDGRAPHTOBIGWIG_SHA256_arm64=3c61e496f446d9965f4b945ca3c5ebfa1a669e920635fbb89a2b705c8a50f5d6

RUN apt-get update && apt-get install -y --no-install-recommends \
      python3 python3-flask python3-gunicorn python3-requests \
      samtools tabix bedtools python3-pybigwig python3-deeptools ca-certificates wget \
      circos liblist-moreutils-perl ncbi-blast+ \
 && case "$TARGETARCH" in \
      amd64) ucsc_dir=linux.x86_64;  wig_sha="$UCSC_WIGTOBIGWIG_SHA256_amd64";  bgr_sha="$UCSC_BEDGRAPHTOBIGWIG_SHA256_amd64" ;; \
      arm64) ucsc_dir=linux.aarch64; wig_sha="$UCSC_WIGTOBIGWIG_SHA256_arm64";  bgr_sha="$UCSC_BEDGRAPHTOBIGWIG_SHA256_arm64" ;; \
      *) echo "unsupported TARGETARCH: $TARGETARCH" >&2; exit 1 ;; \
    esac \
 && base="https://hgdownload.soe.ucsc.edu/admin/exe/$ucsc_dir" \
 && wget -O /usr/local/bin/wigToBigWig      "$base/wigToBigWig" \
 && wget -O /usr/local/bin/bedGraphToBigWig "$base/bedGraphToBigWig" \
 && echo "$wig_sha  /usr/local/bin/wigToBigWig"      | sha256sum -c - \
 && echo "$bgr_sha  /usr/local/bin/bedGraphToBigWig" | sha256sum -c - \
 && chmod +x /usr/local/bin/wigToBigWig /usr/local/bin/bedGraphToBigWig \
 && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Unprivileged runtime user (uid 1000) owning the data mount.
# The named volume `genome-data` is mounted at /genome and holds the SQLite DB,
# stored BigWig/feature files, and the conversion scratch area.
RUN useradd -u 1000 -m genome \
 && mkdir -p /genome \
 && chown genome:genome /genome

# ---------------------------------------------------------------------------
# Module source. `COPY . /srv/genome` with the accompanying .dockerignore keeps
# this maintainable as modules are added/removed; caches and tests are excluded.
WORKDIR /srv/genome
COPY . /srv/genome
RUN chown -R genome:genome /srv/genome

ENV URL_PREFIX=/genome

USER genome

EXPOSE 5008

# Health probe hits the module's unauthenticated health endpoint (the only
# route reachable without proxy auth). wget is already present from apt.
#
# NOTE: the path is /health, NOT /genome/health. The proxy mounts the module at
# /genome and strips that prefix before the request reaches gunicorn, so inside
# the container Flask serves the route at /health (see app.py register_health).
# A /genome/health probe would 404 here. This differs from the design sketch,
# which showed the externally-visible proxied path.
HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
  CMD wget -qO- http://127.0.0.1:5008/health > /dev/null || exit 1

# Web container entrypoint. The worker container overrides this in compose with
# `python3 worker.py`.
# Circos plots (circos_plot.py) can take a minute on a big genome: longer than gunicorn's 30 s default.
CMD ["python3","-m","gunicorn","--bind","0.0.0.0:5008","--workers","2","--threads","4","--timeout","180","app:app"]
