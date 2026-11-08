# Changelog

Releases of the Genome Browser. The [Trellis Tools](https://github.com/Gardiner-Lab/Trellis-Tools)
suite pins a version of this tool as `GENOME_VERSION`, so a release here changes no installation
until the suite is updated to use it.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project uses [semantic versioning](https://semver.org/).

## [Unreleased]

## [0.8.0] - 2026-11-08

### Added
- **Analyse.** A panel beside Circos plot for what the tracks say across the whole genome, each result opening
  like a figure (Lab Notebook, SVG, and its data):
  - **profiles around features**, with a **heatmap** of every feature (sorted or clustered): the coverage
    tracks chosen across genes, transposons or any feature type of a gene annotation (or a BED track's regions),
    scaled or lined up on starts or ends, optionally only a list of genes;
  - **difference tracks**: A minus B, log2(A / B), A / B or their mean, kept as a new coverage track of the lab;
  - **replicate checks**: a correlation heatmap and a PCA of the tracks along the genome;
  - **signal per feature**: each track's mean over every gene (or other feature), as a CSV table, with the
    features that differ most between two tracks shown;
  - **bar charts** of the tracks over the selection, and **synteny dot plots** of a links track.
  The first four are drawn by deepTools on the server, as jobs the genome worker runs between conversions (up
  to three waiting per person; results kept seven days); the page shows how they are getting on.

- **The code to make a figure again, in the notebook.** Saving any figure to the Lab Notebook adds a page under it,
  How to make it again: the inputs (links that download each track under the name the code uses), the software
  and versions, and the code with its files attached. A browser view comes with an igv.js page and its session
  (and the commands to index the reference and annotation); a Circos plot with Circos's configuration and data;
  a profile, heatmap or replicate check with the deepTools commands as run (`run.sh`) and the regions used; a bar
  chart or dot plot with a Python script. The tests run a table's recipe and get the same table back.
- The genome's track list gives each track's download name (`export_name`).

### Fixed
- A page opened straight into the browser (a link, a reload) now reads the genome's chromosomes when a feature
  needs them.

## [0.7.0] - 2026-11-08

### Added
- **Circos plots.** **Circos plot**, beside Save to notebook, draws the genome round a circle with every
  track shown, outermost first, in the tracks' own colours: coverage as its mean signal along the genome,
  a gene annotation as genes per step, other feature tracks as features per step, links as arcs across
  the circle. Choose:
  - each ring's style (histogram, line or heatmap) and height (thin, normal, tall);
  - the whole genome (small chromosomes left out or not), chosen chromosomes, or a region;
  - coarse, medium or fine detail, or an exact step size;
  - a shared scale for the coverage rings, so their heights compare, and a log scale;
  - the selection highlighted through the rings, genes labelled by name, and a title.
  The plot opens like a figure, with a legend and notes (a gene not found, links left out): save it to the
  Lab Notebook, download the SVG, or download Circos's own files to change it and run Circos again. The
  choices are kept for the next plot of the genome. Drawn by Circos (circos.ca, GPL), now installed in
  the image; credited and cited on the About page and in NOTICE; explained on the Help page.
- **Links tracks (BEDPE).** A track of links between two regions (synteny blocks, duplications,
  structural variants): `chrom1 start1 end1 chrom2 start2 end2 [name score]`, up to 100 MB. Both ends
  are checked against the genome and their chromosome names matched, as for other tracks. Drawn as
  arcs in the browser and across a Circos plot (the 25,000 highest-scoring links, at most).

### Changed
- The bundled genome is fetched first in the image build, checked by `ADD --checksum`, so a change to
  the packages no longer downloads it again.

## [0.6.0] - 2026-11-08

### Changed
- **Sharing can be changed at any time, not only set once.** A genome's own chair group, which sees it
  from the moment it is uploaded, now has **Remove** in the Share dialog like any other share: the
  genome is then seen by its owner lab alone (plus anything shared on purpose), and the chair group
  can be added back from the list. Shares with other labs and chair groups were already added and
  removed one at a time; the help says so.

## [0.5.0] - 2026-11-08

### Demo
- **The teaching demo sees only the demo's.** Someone in the suite's demo sees the suite's reference genomes and tracks and the demo labs' own; not a real lab's genome or track, even one shared with them, and is offered only the demo's labs to share with.

### Added
- **Comes with Arabidopsis thaliana.** On its first start the worker adds TAIR10 with the Araport11 genes and transposons (2025-08-13) as the suite's reference genome, for every lab (bundled.py). An installation that has a suite TAIR10 already keeps its own, a deleted one is not added back, and `GENOME_BUNDLED_ARAPORT11=0` leaves it out. The data is TAIR's, under CC BY 4.0, credited on the About page, in the README and in NOTICE.

## [0.4.0] - 2026-11-07

### Added
- **Track search.** A search box above the track list filters it as you type. Every word must
  appear in a track's name or kind (*rnaseq rep2*, *bigwig*, *peaks*), and the matches are
  highlighted, with a count (*3 of 19 tracks*). **Tick shown** ticks everything found. Matches on
  the member's other genomes are listed below, each with **Go to that genome**, which switches there
  with that track ticked. Track settings has its own search over the open tracks, and its *Select
  all* ticks only what that search shows.

- **A track's colour as a hex code** (`#1f5fa8`, `1f5fa8` or the short `#15a`) in Bulk track
  settings, beside the swatches and the colour picker, which follow it.
- **A click on the tracks outside the selection clears it**, as clicking away does anywhere else.
  A click inside the selection, or one that ends a drag along the genome, leaves it.

- **Sharing between labs, from the page.**
  - **Share genome** beside the genome list, and **Share** on every track, open a dialog listing
    who sees it and why: owner lab, its chair group by default, everyone, shared, or *through a
    shared track*.
  - A member of the owning lab can share with another lab or a whole chair group, and remove a
    share. Anyone else sees the list but cannot change it.
  - `GET /api/labs`, `GET /api/genomes/<id>/shares` and `GET /api/tracks/<id>/shares` feed the
    dialog. A lab's registration code is never read.
- **A genome can be shared with a single lab** (`labs` in `POST|DELETE /api/genomes/<id>/share`,
  new table `genome_lab_shares`), not only with chair groups.
- **A shared track takes its genome with it.** Sharing a track with a lab or chair group that
  cannot see its genome shares the genome with them too (`auto = 1`, *through a shared track*), and
  the track share's answer says so (`genome_shared_with`).
  - Once no shared track on the genome needs it, that share goes again.
  - If the genome's lab also shares it on purpose, it stays, and un-sharing the genome on purpose
    leaves the automatic share for as long as a track needs it (`kept`).
  - Suite genomes need none of this.

### Changed
- A long track list scrolls inside itself, so *Start at* and *Open in browser* stay within reach
  below it.
- *Track settings* is now **Bulk track settings**, and *Figure* is now **Save to notebook**. Both
  moved from the far right of the toolbar to the left of the selection bar, before *Select region*,
  where they are seen. The settings panel opens under that bar.

## [0.3.0] - 2026-11-07

A browser that works end to end: choose what to show, select sequence and send it on, add
genomes and tracks (many at once), change several tracks together, and save figures to the Lab
Notebook, with every upload checked against its genome.

### Added
- **A landing page, and a browser that shows something.** `/genome/` used to start igv.js with
  no genome (`genome: null`) and offer no way to choose one, so it stayed blank. Now it opens on
  *Choose what to show*:
  1. pick a genome (suite references are marked *Everyone*);
  2. tick its tracks (Select all / none);
  3. say where to start, by typing a region or clicking a chromosome. A name the genome doesn't
     have is caught before the browser opens.

  *Open in browser* then loads igv.js on that genome with those tracks. The view (place, tracks,
  their order) is saved as you go, so `/genome/?genome=<id>` reopens where you left off, and
  *Change genome or tracks* goes back with your choices kept. The page uses the suite's header and
  look (`static/trellis/`, from Trellis Tools' `design/sync.sh`).
- **The viewer reads indexed files by range request** (`serving.py`):
  - `/api/genomes/<id>/reference.fa` and `.fai`: the sequence, indexed with `samtools faidx`. A
    FASTA with uneven line lengths is re-wrapped once so it can be indexed; the sequence is
    unchanged.
  - `/api/tracks/<id>/coverage.bw`: a BigWig track's file.
  - `/api/tracks/<id>/features.<bed|gff3>.gz` and `.tbi`: a sorted, bgzipped, tabix-indexed copy
    of a BED/GFF, so igv.js draws GFF gene models (exons, UTRs, CDS joined into transcripts) and
    fetches only the window on screen.

  These files are made on upload, or on first request for data uploaded before. Every route checks
  the same visibility as the rest of the API. The image gains Debian's `tabix` (bgzip, tabix).
- `GET /api/genomes` also returns `suite`, `chromosome_count` and `length_bp` for each genome,
  for the landing page.
- **Add a genome / Add a track** on the landing page. A FASTA (gzipped or not) becomes a genome
  of your lab; an administrator can also make it a suite reference. A GFF, BED, BigWig, BAM, wig
  or bedGraph becomes a track on the chosen genome. Uploads show progress and report renamed
  chromosomes, and a conversion is followed until it is ready. A gzipped FASTA is now unpacked on
  upload.
- **Many tracks at once.** *+ Add tracks* takes any number of files, chosen or dropped.
  - Each file gets a row with a name to edit (starting from the filename), its format and size,
    and × to leave it out. A file that isn't a track format is marked and skipped.
  - The files go one after another, each with its own progress bar and result: *Added* (with any
    renamed chromosomes), *Added, converting*, or *Refused* with the reason and, where it fits
    another genome, *Add it to … instead*.
  - A summary line counts them, and the added tracks arrive ticked. Those added cleanly are listed
    as one line under the tracks, however many there were.
- **Track settings: several tracks at once.** *Track settings* in the browser's toolbar lists the
  open tracks with checkboxes (Select all / none). What is set applies to every ticked track:
  - **height**: Short, Medium, Tall, or any number of pixels;
  - **colour**: nine swatches, any colour, or back to the original colours;
  - **gene display**: expanded, squished or collapsed (shown when a feature track is ticked);
  - **coverage scale**: each its own, *shared* (one scale across the ticked tracks so their heights
    compare), or a fixed min–max (shown when a coverage track is ticked);
  - **remove from view**.

  The settings are kept in the saved view, including changes made from a track's own ⚙ menu, so
  a genome reopens looking the same. igv.js's own multi-select toggle, which drops height and
  colour for mixed track types and kept nothing, is hidden in favour of this.
- **Figures to the Lab Notebook.** *Figure* in the browser's toolbar previews the view, as igv.js
  draws it, and offers:
  - **Save to Lab Notebook**: a new page under today in your notebook, through the notebook's own
    API as you. The page holds a line naming the genome, position, tracks and any selection, your
    note, the figure as a picture (a PNG at twice the screen size, since the notebook shows only
    raster images on its pages), and the SVG attached for download. *Open the page →* goes there.
  - **Download SVG**: the figure straight to your computer.
  - A notebook that cannot be reached leaves the SVG downloadable.
- **A track goes to the genome it was made for.** A refused upload lists the member's other
  genomes it does fit (`fits_genomes`), offered as *Add it to … instead*. A BED/GFF that only
  partly fits the chosen genome, while another genome has all of it, is refused and sent there
  rather than accepted with most of its features hidden.
- **Select any sequence, and send it on.**
  - **Shift-drag** across the tracks selects any stretch. **Select region** in the selection bar
    makes a plain drag select instead of moving along the genome. Coordinates can also be typed
    (`Chr1:3,631-5,899`). Double-clicking a feature still selects exactly that feature, without
    zooming, in either mode. A region drawn with igv's own region tool is selected too.
  - The selection is shown as a green band across every track, the bar turns green and reads
    *Selected*, and a band with a running length follows the drag. Escape leaves Select mode, then
    clears the selection.
  - The selection bar offers:
    - *Open in Sequence Editor*, up to 30 kb, with the annotations of the tracks on screen on the
      selection's own coordinates. Longer selections disable it and say why.
    - *BLAST at NCBI* and *Copy sequence*, up to 2 Mb.
    - *Download FASTA*, at any length.
  - A minus-strand feature is sent as it reads: reverse-complemented.
  - Served by `GET /api/genomes/<id>/sequence` (JSON with annotations, up to 2 Mb) and
    `GET /api/genomes/<id>/sequence.fasta` (streamed from `samtools faidx`, any length).
- The suite's tool switcher (`static/trellis/modules.js`) in the header, as in every other tool.
- **Help and About pages** (`/genome/help`, `/genome/about`, `pages.py`), linked from the header
  like every other tool's.
  - Help covers: golden rules, choosing what to show, moving around, the four ways to select,
    what each sending option takes, adding genomes and tracks, chromosome names, why a track is
    refused, who sees what, and what to do when something looks wrong.
  - About says what the tool holds, how tracks are checked, how it works with the other tools,
    where its data is kept (and that BLAST is the one thing that leaves the server), and credits
    igv.js with the paper to cite.
- **Guides in the tool:**
  - a short *New to the Genome Browser?* intro on the landing page;
  - a **Controls** panel in the browser listing every mouse and keyboard action, shown the first
    time and then one click away;
  - a **?** beside each landing-page step and in the upload dialogs, linking to the help;
  - tooltips on the selection buttons saying what each takes and does.

  The intro and the panel stay dismissed once dismissed, remembered in that browser only.
- **Chromosome names are checked on every track upload, and fixed when only the style differs.** A
  BED, GFF, BigWig, BAM or wig whose names differ from the reference genome's only in case, a
  `chr`/`chrom`/`chromosome` prefix, leading zeros or the spelling of an organelle (`MT`, `Pt`,
  `mitochondria`, `chloroplast`…) has them rewritten to the genome's. TAIR10's FASTA (`1`,
  `mitochondria`) and Araport11's GFF3 (`Chr1`, `ChrM`) now go together. A name is renamed only when
  it means one chromosome and no other. BigWigs get a new header, BAMs are reheadered with
  `samtools reheader` before conversion, and text files are rewritten line by line.
- **The upload says what it did.** `chrom_check` (`renamed`, `not_in_genome`, `messages`) comes back
  on the upload response and is kept on the track, returned by `GET /api/tracks/<id>` and its
  `/status`. A refusal now also lists the genome's chromosome names as `genome_chroms`.
- A gzipped BED/GFF is unpacked on upload. Before, it was stored packed and could never be read.
- **Suite reference genomes.** An administrator can upload a genome with `suite_reference=true`.
  It then belongs to the suite rather than a lab (`owner_lab` `suite`, visibility `suite`), every
  signed-in member sees it, and the administrator needs no lab of their own. Tracks an
  administrator adds to it are suite-wide too; anyone else's track on it stays their lab's.
- **Feature tracks are indexed.** Each stored BED/GFF gets a SQLite index beside it
  (`<file>.idx.sqlite`), built on upload, or on first read for a track stored before. A region read
  is now a lookup instead of a parse of the whole file: for Araport11's 289 MB GFF3, a few
  milliseconds instead of re-reading 930,000 features on every pan and zoom.

### Changed
- **BED and GFF are now checked against the genome.** Before, they were stored without any check. A
  file whose names are all missing from the genome is refused; one with only some missing is
  accepted, and `chrom_check` names the rest. A feature past the end of its chromosome is refused.
- `tracks` gains a `chrom_check` column, added to an existing database on start.

### Removed
- `POST /api/figure/notebook` and `notebook.py`. It posted figures to a Lab Notebook endpoint
  (`/api/entries`) that never existed, so it never worked; *Figure* in the browser saves to the
  notebook through the notebook's own API instead.

### Fixed
- **A line of a track file could be split in two.** The readers split on `str.splitlines`, which
  also breaks on form feeds, the ASCII separators and Unicode line separators. Araport11's note
  on AT5G64510 has a 0x1E in it, and its tail was read as a feature on a chromosome "c". Only a
  newline ends a line now.

## [0.2.0] - 2026-11-06

### Added
- **The interactive browser actually runs.** igv.js 3.8.10 is vendored into `static/`, which the
  page has asked for since 0.1.0 (and since suite 2.33.0 before that) while the file itself was
  never added: every visitor got the "igv.min.js did not load" note instead of a viewer. It is
  served from `/genome/static/`, never a CDN, so an isolated network is still enough.
- `NOTICE` records igv.js and the libraries its bundle carries (DOMPurify, vanilla-picker,
  regenerator-runtime, canvas2svg, portions of GMOD two-bit.js), and `static/igv.js-LICENSE.txt`
  carries the MIT notice that has to travel with the code. `NOTICE` also gives the igv.js paper to
  cite, since views from this tool are meant for figures.
- Tests over the page and the vendored bundle: the file is there and is the bundle, the licence is
  beside it, the page boots it from this tool's own mount with nothing fetched from the internet,
  and the asset route serves it while still refusing to climb out of `static/`. None of this was
  covered, which is how the missing file went unnoticed for five releases.

### Changed
- The note shown when the viewer cannot start says what it now means — the file is missing from the
  image — rather than telling the reader it has yet to be vendored.

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
