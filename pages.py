"""
The Genome Browser's Help and About pages, in the suite's format.

Every tool in Trellis Tools has the same two pages, linked from its header:
*Help* (how to use it, task by task) and *About* (what it is, how it fits with
the other tools, where its data lives). These are this tool's. They are static
HTML with the suite's header and look (``static/trellis/``) and the help layout
in ``static/browser.css``; the only values filled in are the mount, the asset
prefix and the version.
"""

from html import escape


def _page(title: str, sub: str, static: str, prefix: str, body: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title} · Genome Browser</title>
  <link rel="icon" href="{static}/trellis/genome.svg" type="image/svg+xml">
  <link rel="stylesheet" href="{static}/trellis/trellis.css">
  <link rel="stylesheet" href="{static}/browser.css">
</head>
<body class="tt">
  <header class="tt-bar">
    <a class="tt-home" href="/" title="All of Trellis Tools"><img src="{static}/trellis/trellis.svg" alt=""><span>Trellis Tools</span></a>
    <a class="tt-module" href="{prefix}/"><img src="{static}/trellis/genome.svg" alt="">Genome Browser</a><span class="tt-sub">{sub}</span>
    <nav class="tt-links">
      <a href="{prefix}/help">Help</a>
      <a href="{prefix}/about">About</a>
    </nav>
  </header>
  <main class="help">
{body}
  </main>
  <script src="{static}/trellis/modules.js" defer></script>
</body>
</html>
"""


def _glyph(static: str, name: str) -> str:
    return f'<img class="tt-glyph" src="{static}/trellis/glyphs/{name}.svg" alt="">'


def help_html(static_prefix: str, prefix: str) -> str:
    s, p = escape(static_prefix), escape(prefix)
    body = f"""
    <h1>How to use the Genome Browser</h1>
    <p class="lead">Look at a reference genome with your lab's data laid over it: gene annotations, peaks,
       coverage from sequencing. Select any stretch of sequence and take it to the Sequence Editor, to
       BLAST, or away as FASTA. This page covers each step.
       <a href="{p}/about">What the Genome Browser is →</a></p>

    <nav class="help-toc">
      <a href="#rules">Golden rules</a>
      <a href="#start">Choosing what to show</a>
      <a href="#move">Moving around</a>
      <a href="#select">Selecting sequence</a>
      <a href="#send">Sending sequence on</a>
      <a href="#figure">Saving a figure</a>
      <a href="#circos">Circos plots</a>
      <a href="#analyse">Analyse: profiles, heatmaps, comparisons</a>
      <a href="#genome">Adding a genome</a>
      <a href="#track">Adding tracks</a>
      <a href="#names">Chromosome names</a>
      <a href="#refused">When a track is refused</a>
      <a href="#sees">Who sees what, and sharing</a>
      <a href="#trouble">When something looks wrong</a>
    </nav>

    <section id="rules">
      <h2>Golden rules</h2>
      <ol class="rules">
        <li><b>A track belongs to the genome it was made for.</b> Reads aligned to TAIR10 go on TAIR10, not
            on another assembly of the same plant. The browser checks, and offers the right genome when it
            can.</li>
        <li><b>Look before you add a genome.</b> If the assembly is already in the list, use it: shared
            genomes (marked <i>Everyone</i>) are there for the whole suite.</li>
        <li><b>Name tracks so a colleague knows what they are:</b> <i>RNA-seq, root, 6 h ABA, rep 2</i>,
            not <i>sample3</i>.</li>
        <li><b>Upload coverage if you have it.</b> A BigWig is ready at once; a BAM is converted on the
            server first, which takes a while for a large one.</li>
        <li><b>Keep the header on sequence you take away.</b> It says which assembly and coordinates it came
            from, which is what makes it findable again.</li>
      </ol>
    </section>

    <section id="start">
      <h2>Choosing what to show</h2>
      <p>The Genome Browser opens on <b>Choose what to show</b>:</p>
      <ol class="steps">
        <li><b>Genome</b>: pick a reference. <i>Everyone</i> marks a suite genome, shared with all labs;
            <i>Your lab</i> marks one of your own lab's or chair group's.</li>
        <li><b>Tracks</b>: tick the data to lay over it. They are drawn in this order, top to bottom. The
            sequence and the ruler are always there. With many tracks, <b>search</b> above the list: type
            words from the name or the kind (<i>rnaseq rep2</i>, <i>bigwig</i>, <i>peaks</i>), then
            <b>Tick shown</b> to tick everything found. A match on another of your genomes is listed below,
            with <b>Go to that genome</b>.</li>
        <li><b>Start at</b>: a chromosome (click one) or a region such as <code>Chr1:3,600-9,000</code>.</li>
        <li><b>Open in browser</b>.</li>
      </ol>
      <p>The browser remembers where you were and which tracks you had, for each genome, so next time it
         opens where you left off. <b>Change genome or tracks</b> takes you back to these choices.</p>
    </section>

    <section id="move">
      <h2>Moving around</h2>
      <table class="help-keys">
        <tr><th>Move along the genome</th><td>Drag the tracks left or right.</td></tr>
        <tr><th>Zoom in or out</th><td>The <b>+</b> and <b>&minus;</b> buttons or the slider in the browser's
            toolbar. Double-clicking empty space zooms in on that spot.</td></tr>
        <tr><th>Go somewhere</th><td>Type a region (<code>Chr2:9,545,000-9,549,000</code>) or a chromosome
            in the box in the browser's toolbar, or pick a chromosome from the list beside it.</td></tr>
        <tr><th>See a feature's details</th><td>Click it once: its name, type, position and the rest of its
            annotation appear in a box.</td></tr>
        <tr><th>Change how tracks look</th><td><b>Bulk track settings</b>, at the left of the bar above the
            browser: tick the tracks to
            change (one, several, or <i>Select all</i>) and set their <b>height</b>, <b>colour</b> (a swatch, the colour picker, or a hex code such as
            <code>#1f5fa8</code>), how genes
            are drawn (<i>expanded</i>, every transcript on its own row; <i>squished</i>; <i>collapsed</i>) and,
            for coverage, the <b>scale</b>: each its own, <i>shared</i> so their heights compare, or a fixed
            range. All ticked tracks change together, and the settings are kept with your view. One track's
            gear (<b>⚙</b>) at its right does the same for that track alone: for a hex code there, choose
            <b>Set track color</b>, then <b>More Colors</b>, click in the box under the colour field (it shows
            the colour as <code>#rrggbb</code>), type the code and press <b>Enter</b> or <b>Ok</b>.</td></tr>
        <tr><th>Reorder or remove a track</th><td>Drag it by the grey handle at its right edge; remove it
            from its gear menu.</td></tr>
      </table>
      <p class="help-note">Gene annotations are drawn within 2 Mb. Zoomed out further than that, the track
         says so: zoom in to see the genes.</p>
    </section>

    <section id="select">
      <h2>Selecting sequence</h2>
      <p>What is selected is shown as a <b>green band</b> down every track, and the bar above the browser
         turns green and reads <b>Selected</b>, with its position and length. There are four ways to
         select:</p>
      <div class="help-cards">
        <div class="help-card">
          <h3>Double-click a feature</h3>
          <p>Selects exactly that gene, transcript or other feature: its full extent, on its own strand.
             (On a feature, double-clicking selects rather than zooms.)</p>
        </div>
        <div class="help-card">
          <h3>Shift + drag</h3>
          <p>Hold <kbd>Shift</kbd> and drag across the tracks to select any stretch. A band follows the
             mouse and shows its length as you go.</p>
        </div>
        <div class="help-card">
          <h3>Select region</h3>
          <p>Turn on <b>Select region</b> in the bar and a plain drag selects instead of moving along.
             <kbd>Esc</kbd> turns it off again.</p>
        </div>
        <div class="help-card">
          <h3>Type it</h3>
          <p>Type a region into the box in the bar (<code>Chr1:3,631-5,899</code>) and press <b>Select</b>.
             The browser goes there too.</p>
        </div>
      </div>
      <p><b>Clear</b> in the bar, <kbd>Esc</kbd>, or a click on the tracks outside the selected stretch removes the selection. A region drawn with the
         browser's own <b>ROI</b> tool in its toolbar becomes the selection as well.</p>
    </section>

    <section id="send">
      <h2>Sending sequence on</h2>
      <table class="help-keys">
        <tr><th>Open in Sequence Editor</th><td>Up to <b>30 kb</b>. Opens it as a new sequence in the
            Sequence Editor, in a new tab, with the annotations of the gene tracks on screen placed on it:
            exons, CDS, UTRs and so on. It is saved to your sequences there, named after the feature.</td></tr>
        <tr><th>BLAST at NCBI</th><td>Up to 2 Mb. Opens NCBI's nucleotide BLAST in a new tab with the
            sequence filled in; choose the database there and run it. <b>This sends the sequence to NCBI</b>,
            outside the lab's server.</td></tr>
        <tr><th>Copy sequence</th><td>Up to 2 Mb, to paste anywhere.</td></tr>
        <tr><th>Download FASTA</th><td><b>Any length</b>, even a whole chromosome. Anything too long for the
            others goes this way.</td></tr>
      </table>
      <p>A feature on the minus strand is sent <b>as it reads</b>: reverse-complemented, with its
         annotations turned with it. A dragged or typed region is the plus strand. The FASTA header names
         the feature, the assembly and the coordinates, for example
         <code>&gt;AT2G22480.1 TAIR10 Chr2:9545283-9548640</code>; <code>(-)</code> marks the minus strand.</p>
    </section>

    <section id="figure">
      <h2>Saving a figure</h2>
      <p><b>Save to notebook</b>, at the left of the bar above the browser, saves exactly what is on screen (the tracks, the
         position, and any selection band) as a figure:</p>
      <ol class="steps">
        <li>Get the view right first: position, tracks, their heights and colours (each track's ⚙).</li>
        <li>Click <b>Save to notebook</b>. A preview shows what will be saved.</li>
        <li>Give the page a <b>title</b> (it starts with the gene, the assembly and the position) and, if
            you like, a <b>note</b> saying what it shows.</li>
        <li><b>Save to Lab Notebook</b> makes a new page <b>under today</b> in your notebook, with the
            picture, a line saying which genome, position and tracks it shows, your note, and the figure as
            an <b>SVG</b> file. <b>Open the page →</b> takes you there.</li>
      </ol>
      <p><b>Every figure saved to the notebook gets a page under it: How to make it again.</b> It lists the inputs, each a
         link that downloads the track under the name the code uses, the software to install, and the code, with its files
         attached, so someone outside the suite can make the same figure:</p>
      <ul>
        <li><b>A view of the browser</b>: an igv.js page (<code>figure.html</code>) and the view (<code>session.json</code>:
            the place, the tracks, their colours, heights and scales, the selection), the commands to index the
            reference and annotation files, and how to open it with a local web server and save it as SVG.</li>
        <li><b>A Circos plot</b>: Circos's own configuration and data, zipped, and how to install and run Circos.</li>
        <li><b>A profile, heatmap or replicate check</b>: the deepTools commands exactly as they were run, as
            <code>run.sh</code>, with the regions used as BED files; the inputs are named as they download.</li>
        <li><b>A bar chart or dot plot</b>: a Python script that draws it again from the downloaded tracks.</li>
      </ul>
      <p><b>Download SVG</b> saves the figure straight to your computer instead. SVG is a vector format:
         it stays sharp at any size and its text and lines can be edited in Inkscape or Illustrator, which
         is what a journal figure needs. In the notebook, the SVG is attached for download, and a picture
         of it shows on the page.</p>
    </section>

    <section id="circos">
      <h2>Circos plots</h2>
      <p><b>Circos plot</b>, beside <b>Save to notebook</b>, draws the genome round a circle with every track you are
         showing as a ring, using <a href="https://circos.ca" target="_blank" rel="noopener">Circos</a>. Show the tracks
         you want, in the order you want, in the colours you want, then open it.</p>
      <ol class="steps">
        <li>Give it a <b>title</b>, or leave the genome and what is drawn as the title.</li>
        <li>For each track, choose how its ring is drawn: <b>histogram</b>, <b>line</b> or <b>heatmap</b> (a band
            shaded from white to the track's colour: compact, so many tracks fit), and how tall: <b>thin</b>,
            <b>normal</b> or <b>tall</b>.</li>
        <li>Choose what to draw: <b>the whole genome</b> (leaving out small chromosomes, under 1% of the genome, or
            not), <b>chosen chromosomes</b>, or <b>a region</b> such as <code>Chr1:1,000,000-5,000,000</code> (filled in
            with your selection if you have one, else the stretch in view; a region must be at least 1 kb).</li>
        <li>Choose the <b>detail</b> (coarse, medium, fine) or give an exact step size in base pairs.</li>
        <li>Optionally: let coverage rings <b>share one scale</b>, put coverage on a <b>log scale</b>,
            <b>highlight the selection</b> you made in the browser, and <b>label genes</b> by name.</li>
        <li><b>Draw</b>. A big genome with fine detail takes up to a minute. Your choices are kept for the next
            plot of this genome.</li>
      </ol>
      <h3>How it is drawn</h3>
      <ul>
        <li><b>The outer ring is the genome</b>: each chromosome as a band in proportion to its length, named, with
            ticks along it in megabases (kilobases for a small genome).</li>
        <li><b>The genome is cut into equal steps</b> (bins) all the way round: about 300 for coarse, 800 for medium
            and 2,000 for fine, the size rounded to 1, 2 or 5 times a power of ten. For Arabidopsis at medium that
            is 200 kb. The plot says which size it used.</li>
        <li><b>Each track is one ring</b>, inside the genome, in the order the tracks are shown (the top track
            outermost), in its own colour, as a histogram: one bar per step.
          <ul>
            <li><b>Coverage</b> (BigWig): the mean signal over the step. A step with no data has no bar.</li>
            <li><b>A gene annotation</b> (GFF): the number of genes that start in the step. Only features of type
                <i>gene</i> are counted, so each gene counts once, not again for its mRNAs and exons; a GFF with no
                genes has every feature counted.</li>
            <li><b>Other feature tracks</b> (BED): the number of features that start in the step.</li>
          </ul></li>
        <li><b>Every ring has its own scale</b>, from nothing at its inner edge to its own highest value at its
            outer edge, with a faint line halfway. So a ring shows where along the genome its track is high or low;
            heights in two rings are not comparable, unless coverage rings <b>share one scale</b>: then they all run
            to the highest value among them. A <b>log scale</b> draws log10(1 + the mean), so a few high peaks do not
            flatten everything else.</li>
        <li><b>Links</b> (BEDPE tracks) are drawn as arcs across the middle of the circle, inside the rings, in
            the track's colour, from one region to the other. Only links with both ends in what is drawn are shown;
            of very many, the 25,000 with the highest score (the file's eighth column) are drawn, and the plot says
            so.</li>
        <li><b>Gene labels</b> are written just inside the chromosomes, with a short line to each gene; names close
            together are spread out. They are looked up in the gene annotation tracks shown, by the gene's name
            there (for Araport11, its locus, such as <code>AT1G01010</code>); a name not found is listed with the
            plot.</li>
        <li><b>The selection</b>, if you tick <b>Highlight the selection</b> (it is off to start with), is shaded
            pale yellow through all the rings.</li>
        <li><b>A legend</b> in the corner names the rings, outermost first, and what each shows.</li>
        <li>Only tracks you may see are drawn, at most 12 (links included); a track still being converted is not. A track's chromosome
            names are matched to the genome's as on upload (<code>chr1</code> is <code>Chr1</code>).</li>
      </ul>
      <h3>Afterwards</h3>
      <p>The plot opens like a figure: save it to your Lab Notebook (a page under today, with the picture, the SVG
         and a note of what it shows), or download it as SVG. <b>Download Circos files</b> gives Circos's own files:
         <code>circos.conf</code> (the rings, colours and ticks), <code>karyotype.txt</code> (the chromosomes),
         <code>data/track-1.txt</code> and so on (one line per step: chromosome, start, end, value) and a
         <code>README.txt</code>. Change them (a ring's colour or height, a scale, a heatmap instead of a histogram)
         and run <code>circos -conf circos.conf</code> on a computer with Circos installed. If a Circos plot goes
         into a paper, cite Circos (see About).</p>
    </section>

    <section id="analyse">
      <h2>Analyse: profiles, heatmaps, comparisons</h2>
      <p><b>Analyse</b>, beside <b>Circos plot</b>, draws what the tracks say across the whole genome, not just the
         stretch on screen. Choose what to draw; the coverage tracks you are showing are ticked to start with. Each result
         opens like a figure: save it to your Lab Notebook, download it as SVG, and download its data.</p>
      <div class="help-cards">
        <div class="help-card">
          <h3>Profile around features, and heatmap</h3>
          <p>The average signal of each coverage track across a set of features: every gene (or transposon, or other
             type from a gene annotation, or the regions of a BED track), scaled from start to end, with 1 to 10 kb either
             side, or lined up on their starts or ends. Up to four feature types, each its own panel. With
             <b>Heatmap</b>, every feature is a row too, sorted by signal or grouped into clusters, so you see whether a
             change touches every gene a little or some genes a lot. Limit it to a list of genes with <b>Only these</b>.</p>
        </div>
        <div class="help-card">
          <h3>Difference between two tracks</h3>
          <p>A new coverage track of track A against track B: A minus B, log2(A / B), A / B, or their mean, in steps of
             50 bp to 1 kb. It is kept for your lab like an upload, shown at once, and can be drawn in Circos plots
             and profiles like any track. For a ratio, the pseudocount is added to both so that zero does not divide.</p>
        </div>
        <div class="help-card">
          <h3>Replicate check</h3>
          <p>The chosen tracks compared along the whole genome in bins (5 to 100 kb): a heatmap of their correlations
             (Pearson, or Spearman on ranks), clustered, with the values written in, and a principal component plot.
             Replicates should sit together and mutants apart; one that does not is worth a look before trusting it.</p>
        </div>
        <div class="help-card">
          <h3>Signal per feature</h3>
          <p>The mean of each coverage track over every feature of one type, as a table (CSV) with each feature's name
             and position: one row per gene, say. With two tracks, their difference too, and the fifty features that differ
             most are shown straight away.</p>
        </div>
        <div class="help-card">
          <h3>Tracks over the selection</h3>
          <p>A bar chart of each coverage track's mean over the region you selected, in the tracks' colours.</p>
        </div>
        <div class="help-card">
          <h3>Synteny dot plot</h3>
          <p>A links (BEDPE) track with the genome along both sides: each link a short line from one region to the
             other, drawn both ways round, so blocks that kept their order show as diagonals.</p>
        </div>
      </div>
      <p>The first four run on the server and can take a while: seconds for a replicate check, a few minutes for a
         profile of every gene with several tracks. The panel says how long it has been working and whether it is waiting
         for another job; you can close it and keep browsing, and up to three can wait at once. Results are kept for seven
         days. Only tracks you may see can be used. They are drawn with
         <a href="https://deeptools.readthedocs.io" target="_blank" rel="noopener">deepTools</a>; cite it if a figure goes
         into a paper (see About).</p>
    </section>

    <section id="genome">
      <h2>Adding a genome</h2>
      <ol class="steps">
        <li>On <b>Choose what to show</b>, click <b>+ Add a genome</b>.</li>
        <li>Give it a <b>name</b> people will recognise (<i>Arabidopsis thaliana Col-0 (TAIR10)</i>) and the
            <b>assembly</b> identifier (<i>TAIR10</i>).</li>
        <li>Choose the <b>FASTA</b> file: <code>.fa</code>, <code>.fasta</code>, <code>.fna</code> or
            <code>.fsa</code>, gzipped or not, up to 4 GB. The chromosome names are the first word of each
            <code>&gt;</code> line.</li>
        <li><b>Upload genome</b>. It belongs to your lab and is shared with your chair group.</li>
      </ol>
      <p>An administrator can tick <b>Share with everyone in the suite</b> to make it a <i>suite
         reference</i>: owned by no lab and seen by every member. Use this for the standard assemblies
         everyone works with.</p>
    </section>

    <section id="track">
      <h2>Adding tracks</h2>
      <ol class="steps">
        <li>Choose the genome first, then click <b>+ Add tracks</b>.</li>
        <li><b>Drop the files</b> on the box, or <b>choose files</b>: one, or a whole experiment's worth at
            once. Each gets a row.</li>
        <li>Each row's <b>name</b> starts as the file's name; change it to say what the track is
            (<i>RNA-seq, root, 6 h ABA, rep 2</i>). <b>×</b> leaves a file out. A file the browser
            doesn't take is marked and skipped.</li>
        <li><b>Upload</b>. The files go one after another, each with its own progress and result:
            <i>Added</i>, <i>Added, converting</i>, or <i>Refused</i> with the reason. Keep the window
            open until they are through.</li>
        <li>The added tracks appear in the list, ticked, ready to open.</li>
      </ol>
      <p>Name an <b>Experiment</b> above the files and every file added then is put in it. The track list groups
         tracks under their experiment (<b>Sort by</b> also offers name, or newest first); <b>Tick these</b> beside an
         experiment ticks all of its tracks. A track's <b>Details</b> changes its experiment, tags it with the
         Seed Bank lines, seed lots and plants it comes from, and lists the Lab Notebook pages that link to it.
         Only the track's own lab changes its experiment and tags.</p>
      <p><b>Published data.</b> A track from a paper says so: in its <b>Details</b> (or, for every file added at once,
         under <i>From a paper or GEO?</i> when adding them), give the <b>paper</b> as you would cite it, its
         <b>DOI</b> (<code>10.1016/j.cell.2012.10.054</code>, or its doi.org address) and its <b>GEO accessions</b>
         (series <code>GSE…</code> and samples <code>GSM…</code>). Everyone who sees the track sees them, linked to
         the paper and to GEO; the track list shows the GEO series, and the search finds tracks by any of them.</p>
      <p><b>Removing a track.</b> <b>Details</b>, <b>Remove track…</b> takes it away for good, for everyone: its
         file, tags and sharing. Only its own lab (or an administrator) can, and it cannot be undone, as the
         suite's backups do not keep track files. Notebook pages that link to it then say it is gone.</p>
      <h3>Marking a sequence: guides, primers, anything unlabelled</h3>
      <p><b>+ Mark a sequence</b> looks for a stretch of DNA everywhere on the genome, on both strands, and keeps
         what it finds as a track of its own. It is the way to show things the gene annotation does not name:
         CRISPR guides, primer binding sites, an enzyme site, a motif.</p>
      <ol class="steps">
        <li><b>Name the track</b> &mdash; <i>CRISPR guides, HSP101</i> &mdash; and type the <b>sequences</b>, one per
            line. Put a name before a sequence to label it (<code>sgRNA1&nbsp;&nbsp;GACGTTAGCACCTGTAGGAT</code>);
            a bare sequence takes the track's name. A FASTA <code>&gt;name</code> header works too.</li>
        <li><b>IUPAC codes</b> are allowed, so a degenerate primer or an enzyme site goes in as written:
            <code>GGTCTCNNNNN</code>, <code>GGWCC</code>.</li>
        <li><b>Find and add the track.</b> Each place found becomes a feature, on the strand it reads on, named
            after its sequence and numbered when there is more than one (<i>sgRNA1.1</i>, <i>sgRNA1.2</i>).</li>
      </ol>
      <p>From then on it is an ordinary track: tick it to show it, share it, export it, put it in an experiment
         or remove it. A sequence found nowhere on the genome makes no track and says so &mdash; usually it is
         from another assembly, or it has a vector end on it.</p>
      <p>The search is <b>exact</b>: it finds the sequence as written, so a guide's off-targets with a base or two
         changed are not marked. A sequence must be at least 8 bases. One that is all over the genome (an enzyme
         site) stops at 500 places, and the others typed with it are still searched in full.</p>
      <p>A refused file doesn't stop the rest. If it fits another of your genomes, its row offers
         <b>Add it to … instead</b>.</p>
      <table class="help-keys">
        <tr><th>GFF, GFF3</th><td>Gene annotations and other features, drawn as gene models.</td></tr>
        <tr><th>BED</th><td>Intervals: peaks, regions, primers.</td></tr>
        <tr><th>BEDPE</th><td>Links between two regions: synteny blocks, duplications, structural variants. One link a
            line: <code>chrom1 start1 end1 chrom2 start2 end2</code>, then optionally a name and a score (stronger
            links are kept first when a Circos plot has too many). Shown as arcs, and in a Circos plot as arcs across
            the circle. Both ends must be on this genome. At most 100 MB. Synteny tools (MCScanX, JCVI) list matching
            genes; turn their blocks into regions first, for example with JCVI's <code>jcvi.compara.synteny</code>.</td></tr>
        <tr><th>BigWig</th><td>Coverage or signal, ready at once.</td></tr>
        <tr><th>BAM, wig, bedGraph</th><td>Converted to coverage on the server. It shows as
            <i>converting</i> under the list, then joins the list when ready; you can leave the page
            meanwhile.</td></tr>
      </table>
      <p>Any of them may be gzipped. A track belongs to your lab; an administrator's track on a suite
         genome is shared with everyone.</p>
    </section>

    <section id="names">
      <h2>Chromosome names</h2>
      <p>The same chromosome is written many ways: <code>1</code>, <code>chr1</code>, <code>Chr1</code>,
         <code>chromosome_1</code>; the mitochondrion as <code>M</code>, <code>MT</code>, <code>ChrM</code> or
         <code>mitochondria</code>; the chloroplast as <code>C</code>, <code>Pt</code> or <code>ChrC</code>.
         When a track writes them differently from its genome, it is <b>renamed to match</b> on upload, and
         the browser tells you what it renamed. It only renames a name that can mean one chromosome and no
         other.</p>
    </section>

    <section id="refused">
      <h2>When a track is refused</h2>
      <p>A track is refused, and nothing is kept, when it does not fit the genome:</p>
      <ul>
        <li><b>None of its chromosome names are in the genome.</b> It was made for another assembly. If it
            fits one of your other genomes, the message offers <b>Add it to … instead</b>.</li>
        <li><b>A position runs past the end of a chromosome.</b> Same chromosome names, different assembly
            (an older release, say).</li>
        <li><b>A coverage file names a chromosome the genome doesn't have.</b> BigWig, BAM and wig must
            match completely.</li>
      </ul>
      <p>A GFF or BED with only <i>some</i> chromosomes missing (often unplaced scaffolds) is accepted, and
         the message lists the ones that won't be drawn. If another of your genomes has all of them, it is
         sent there instead.</p>
    </section>

    <section id="sees">
      <h2>Who sees what, and sharing</h2>
      <ul>
        <li>A <b>genome</b> is seen by its lab's chair group, or just its lab if the lab is in none.
            A suite reference is seen by everyone.</li>
        <li>A <b>track</b> is seen by its own lab. An administrator's track on a suite genome is seen by
            everyone.</li>
        <li>Your <b>saved view</b> of each genome (where you were, which tracks) is yours alone.</li>
      </ul>
      <p>To share, a member of the lab that owns it opens <b>Share</b>: <b>Share genome</b> beside the genome
         list, or <b>Share</b> at the right of a track. It lists who sees it and why, and shares it with
         another <b>lab</b> or a whole <b>chair group</b>; <b>Remove</b> takes a share back. It can be changed
         at any time, as often as you like: add and remove one at a time. Even a genome's own chair group,
         which sees it from the start, can be removed, so only your lab sees it, and added back. Anyone
         else sees the list but cannot change it.</p>
      <p>A track is only any use to someone who can see its genome, so <b>sharing a track with a lab that
         cannot see its genome shares the genome with that lab too</b>. That share is marked <i>through a
         shared track</i> and goes again when no track on the genome is shared with the lab any more. If the
         genome's own lab shares it on purpose as well, it stays.</p>
    </section>

    <section id="trouble">
      <h2>When something looks wrong</h2>
      <table class="help-keys">
        <tr><th>A track I added isn't in the list</th><td>It may be on another genome (choose that one), still
            converting (see under the list), or it was refused (the message said why).</td></tr>
        <tr><th>A gene track is empty</th><td>Zoom in to 2 Mb or less, or check you are on a chromosome it
            covers.</td></tr>
        <tr><th>Double-click zoomed instead of selecting</th><td>Double-click on the feature itself (its box
            or line), not the space beside it.</td></tr>
        <tr><th>Open in Sequence Editor is greyed out</th><td>The selection is over 30 kb. Select less, or
            use <b>Download FASTA</b>.</td></tr>
        <tr><th>Dragging selects instead of moving</th><td><b>Select region</b> is on: click it, or press
            <kbd>Esc</kbd>.</td></tr>
      </table>
    </section>
"""
    return _page('Help', 'Help', s, p, body)


def about_html(static_prefix: str, prefix: str, version: str) -> str:
    s, p, v = escape(static_prefix), escape(prefix), escape(version)
    body = f"""
    <h1>Genome Browser {v}</h1>
    <p class="lead">Reference genomes and the lab's data on them, side by side: gene annotations, peaks and
       sequencing coverage, viewed in an IGV-style browser, with any stretch of sequence a click away from
       the Sequence Editor, BLAST or a FASTA file.
       <a href="{p}/help">How to use it →</a></p>

    <nav class="help-toc">
      <a href="#holds">What it holds</a>
      <a href="#check">Every track is checked</a>
      <a href="#sees">Who sees what, and sharing</a>
      <a href="#together">With the other tools</a>
      <a href="#kept">Where the data is kept</a>
      <a href="#credits">Credits and licence</a>
    </nav>

    <section id="holds">
      <h2>What it holds</h2>
      <div class="help-cards">
        <div class="help-card">
          <h3>{_glyph(s, 'library')} Reference genomes</h3>
          <p>An assembly uploaded as FASTA, by a lab for its chair group, or by an administrator for the
             whole suite. It comes with Arabidopsis thaliana (TAIR10, with the Araport11 annotation) for
             everyone.</p>
        </div>
        <div class="help-card">
          <h3>{_glyph(s, 'line')} Tracks</h3>
          <p>Gene annotations and intervals (GFF, BED), coverage (BigWig, or BAM and wig converted on the
             server) and links between regions such as synteny (BEDPE), each on the genome it was made for.</p>
        </div>
        <div class="help-card">
          <h3>{_glyph(s, 'cut')} Sequence to take away</h3>
          <p>Select a feature or any stretch and open it in the Sequence Editor with its annotations, BLAST
             it, or download it as FASTA.</p>
        </div>
        <div class="help-card">
          <h3>{_glyph(s, 'ratio')} Circos plots</h3>
          <p>The whole genome, chosen chromosomes or a region round a circle, each track a ring (histogram, line or
             heatmap), links such as synteny as arcs across it, genes labelled and the selection highlighted. Drawn
             by Circos on the server; saved to the notebook, as SVG, or as Circos's own files to change.</p>
        </div>
        <div class="help-card">
          <h3>{_glyph(s, 'line')} Analyse</h3>
          <p>Metagene profiles and heatmaps around genes or transposons, difference tracks (mutant minus wild type),
             replicate checks (correlation and PCA), signal per gene as a table, bar charts of a selection and synteny
             dot plots, drawn with deepTools.</p>
        </div>
        <div class="help-card">
          <h3>{_glyph(s, 'history')} Your place, kept</h3>
          <p>Each person's view of each genome (where they were, which tracks) is remembered between
             visits.</p>
        </div>
      </div>
    </section>

    <section id="check">
      <h2>Every track is checked</h2>
      <p>A track drawn on the wrong assembly looks plausible and is wrong everywhere, so every upload is
         checked against its genome before it is kept:</p>
      <ul>
        <li>its chromosome names, with names written in another style (<code>chr1</code> for
            <code>Chr1</code>, <code>MT</code> for <code>ChrM</code>) renamed to the genome's;</li>
        <li>its positions, none of which may run past the end of a chromosome;</li>
        <li>and, when it doesn't fit, which of your other genomes it does fit.</li>
      </ul>
    </section>

    <section id="sees">
      <h2>Who sees what</h2>
      <p>As in the rest of the suite: a genome is its lab's and its chair group's, a track is its lab's,
         and suite references, with the tracks an administrator adds to them, are everyone's.</p>
    </section>

    <section id="together">
      <h2>With the other tools</h2>
      <ul>
        <li><b>Sequence Editor</b>: a selection opens there as a new sequence, annotated from the gene
            tracks on screen.</li>
        <li><b>Lab Notebook</b>: <b>Save to notebook</b> saves the view as a page under today, with the picture,
            what it shows, and the SVG to download.</li>
        <li><b>Seed Bank</b>: a track's <b>Details</b> tags it with the lines, seed lots and plants it comes from,
            found by code or line name; each tag opens its record in the Seed Bank, and the record lists the tracks
            and genomes tagged with it (those you can open here).</li>
        <li><b>Lab Notebook</b> pages link to a track: type <b>@</b> in a page and choose <i>A genome track</i>. The
            link opens the track's genome with it on, and the track's <b>Details</b> lists the pages.</li>
      </ul>
    </section>

    <section id="kept">
      <h2>Where the data is kept</h2>
      <p>On your own server, in the Genome Browser's own storage, backed up each night with the rest of the
         suite. Nothing is fetched from the internet to show it. The one thing that leaves the server is a
         sequence you choose to BLAST, which goes to NCBI.</p>
    </section>

    <section id="credits">
      <h2>Credits and licence</h2>
      <p>The browser view is <b>igv.js</b>, from the Broad Institute and the University of California, under
         the MIT licence. If a view from it goes into a paper, cite: Robinson JT, Thorvaldsdottir H, Turner D,
         Mesirov JP. <i>igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer
         (IGV).</i> Bioinformatics (2023).
         <a href="https://doi.org/10.1093/bioinformatics/btac830">doi:10.1093/bioinformatics/btac830</a></p>
      <p>The genome it comes with, <b>Arabidopsis thaliana TAIR10</b> with the <b>Araport11</b> genes and
         transposons, is from <a href="https://www.arabidopsis.org" rel="noopener">The Arabidopsis Information
         Resource (TAIR)</a>, unchanged, under the
         <a href="https://creativecommons.org/licenses/by/4.0/" rel="noopener">Creative Commons Attribution 4.0</a>
         licence. If it goes into a paper, cite: Cheng C-Y <i>et al.</i> <i>Araport11: a complete reannotation of
         the Arabidopsis thaliana reference genome.</i> The Plant Journal 89: 789&ndash;804 (2017).
         <a href="https://doi.org/10.1111/tpj.13415">doi:10.1111/tpj.13415</a>; and Lamesch P <i>et al.</i>
         <i>The Arabidopsis Information Resource (TAIR): improved gene annotation and new tools.</i> Nucleic
         Acids Research 40: D1202&ndash;D1210 (2012).
         <a href="https://doi.org/10.1093/nar/gkr1090">doi:10.1093/nar/gkr1090</a></p>
      <p>Circos plots are drawn by <b>Circos</b>, by Martin Krzywinski, under the GNU General Public License. If a
         Circos plot goes into a paper, cite: Krzywinski M <i>et al.</i> <i>Circos: an information aesthetic for
         comparative genomics.</i> Genome Research 19: 1639&ndash;1645 (2009).
         <a href="https://doi.org/10.1101/gr.092759.109">doi:10.1101/gr.092759.109</a></p>
      <p>Profiles, heatmaps, difference tracks and replicate checks are drawn by <b>deepTools</b> (installed from
         Debian's package, under its own open-source licence). If one goes into a paper, cite: Ram&iacute;rez F <i>et al.</i> <i>deepTools2: a next
         generation web server for deep-sequencing data analysis.</i> Nucleic Acids Research 44: W160&ndash;W165 (2016).
         <a href="https://doi.org/10.1093/nar/gkw257">doi:10.1093/nar/gkw257</a></p>
      <p>The Genome Browser is free software under the GNU General Public License v3.0.</p>
    </section>
"""
    return _page('About', 'About', s, p, body)
