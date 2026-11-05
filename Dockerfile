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
# System packages + bundled conversion toolchain.
#
# The suite prefers Debian packages over pip, so Flask/gunicorn/requests and the
# bio toolchain all come from apt. python3-requests is required because the
# module imports `requests` (e.g. notebook.py export to the Lab Notebook).
#
# Supply chain: the UCSC static binaries are verified against pinned SHA-256s.
# The two UCSC static binaries are fetched at build time from the canonical
# UCSC download host and checked against the SHA-256 values pinned in the ARGs
# below (captured 2026-11-05); the sha256sum -c check fails the build on any
# mismatch, so the download cannot be silently substituted. UCSC periodically
# rebuilds these binaries; if a future build fails the check, re-fetch each,
# confirm its provenance, run `sha256sum <file>`, and update the matching ARG.
ARG UCSC_BASE_URL=https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64
ARG UCSC_WIGTOBIGWIG_SHA256=5b9eccc4929b339be26bb69bce27de63070195cb0a431b410627c0bed8ab4f15
ARG UCSC_BEDGRAPHTOBIGWIG_SHA256=fe9ed8ed858731d22d9a3cefa99d192f51bff1496693599b2f8054fa350adb48

RUN apt-get update && apt-get install -y --no-install-recommends \
      python3 python3-flask python3-gunicorn python3-requests \
      samtools bedtools python3-pybigwig python3-deeptools ca-certificates wget \
 && wget -O /usr/local/bin/wigToBigWig "${UCSC_BASE_URL}/wigToBigWig" \
 && wget -O /usr/local/bin/bedGraphToBigWig "${UCSC_BASE_URL}/bedGraphToBigWig" \
 && echo "${UCSC_WIGTOBIGWIG_SHA256}  /usr/local/bin/wigToBigWig"         | sha256sum -c - \
 && echo "${UCSC_BEDGRAPHTOBIGWIG_SHA256}  /usr/local/bin/bedGraphToBigWig" | sha256sum -c - \
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
CMD ["gunicorn","--bind","0.0.0.0:5008","--workers","2","--threads","4","app:app"]
