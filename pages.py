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
      <a href="#genome">Adding a genome</a>
      <a href="#track">Adding tracks</a>
      <a href="#names">Chromosome names</a>
      <a href="#refused">When a track is refused</a>
      <a href="#sees">Who sees what</a>
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
            sequence and the ruler are always there.</li>
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
        <tr><th>Change how tracks look</th><td><b>Track settings</b> in the toolbar: tick the tracks to
            change (one, several, or <i>Select all</i>) and set their <b>height</b>, <b>colour</b>, how genes
            are drawn (<i>expanded</i>, every transcript on its own row; <i>squished</i>; <i>collapsed</i>) and,
            for coverage, the <b>scale</b>: each its own, <i>shared</i> so their heights compare, or a fixed
            range. All ticked tracks change together, and the settings are kept with your view. One track's
            gear (<b>⚙</b>) at its right does the same for that track alone.</td></tr>
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
      <p><b>Clear</b> in the bar, or <kbd>Esc</kbd>, removes the selection. A region drawn with the
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
      <p><b>Figure</b> in the browser's toolbar saves exactly what is on screen (the tracks, the
         position, and any selection band) as a figure:</p>
      <ol class="steps">
        <li>Get the view right first: position, tracks, their heights and colours (each track's ⚙).</li>
        <li>Click <b>Figure</b>. A preview shows what will be saved.</li>
        <li>Give the page a <b>title</b> (it starts with the gene, the assembly and the position) and, if
            you like, a <b>note</b> saying what it shows.</li>
        <li><b>Save to Lab Notebook</b> makes a new page <b>under today</b> in your notebook, with the
            picture, a line saying which genome, position and tracks it shows, your note, and the figure as
            an <b>SVG</b> file. <b>Open the page →</b> takes you there.</li>
      </ol>
      <p><b>Download SVG</b> saves the figure straight to your computer instead. SVG is a vector format:
         it stays sharp at any size and its text and lines can be edited in Inkscape or Illustrator, which
         is what a journal figure needs. In the notebook, the SVG is attached for download, and a picture
         of it shows on the page.</p>
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
      <p>A refused file doesn't stop the rest. If it fits another of your genomes, its row offers
         <b>Add it to … instead</b>.</p>
      <table class="help-keys">
        <tr><th>GFF, GFF3</th><td>Gene annotations and other features, drawn as gene models.</td></tr>
        <tr><th>BED</th><td>Intervals: peaks, regions, primers.</td></tr>
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
      <h2>Who sees what</h2>
      <ul>
        <li>A <b>genome</b> is seen by its lab's chair group, or just its lab if the lab is in none.
            A suite reference is seen by everyone.</li>
        <li>A <b>track</b> is seen by its own lab. An administrator's track on a suite genome is seen by
            everyone.</li>
        <li>Your <b>saved view</b> of each genome (where you were, which tracks) is yours alone.</li>
      </ul>
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
      <a href="#sees">Who sees what</a>
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
             whole suite.</p>
        </div>
        <div class="help-card">
          <h3>{_glyph(s, 'line')} Tracks</h3>
          <p>Gene annotations and intervals (GFF, BED) and coverage (BigWig, or BAM and wig converted on the
             server), each on the genome it was made for.</p>
        </div>
        <div class="help-card">
          <h3>{_glyph(s, 'cut')} Sequence to take away</h3>
          <p>Select a feature or any stretch and open it in the Sequence Editor with its annotations, BLAST
             it, or download it as FASTA.</p>
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
        <li><b>Lab Notebook</b>: <b>Figure</b> saves the view as a page under today, with the picture,
            what it shows, and the SVG to download.</li>
        <li><b>Seed Bank</b> and plasmids: tracks and genomes can be tagged with the plasmid or seed line
            they come from.</li>
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
      <p>The Genome Browser is free software under the GNU General Public License v3.0.</p>
    </section>
"""
    return _page('About', 'About', s, p, body)
