"""Circos plots: the files made from a genome and its tracks, who may draw what, and (with Circos installed) a plot."""

import base64
import io
import uuid
import zipfile

import pytest

import app as app_module
import circos_plot
import db
import identity
import labs_reader

GFF = ('##gff-version 3\n'
       'Chr1\tx\tgene\t100\t900\t.\t+\t.\tID=g1\n'
       'Chr1\tx\tgene\t1500\t1800\t.\t+\t.\tID=g2\n'
       'Chr1\tx\tmRNA\t1500\t1800\t.\t+\t.\tID=m2;Parent=g2\n'
       'Chr2\tx\tgene\t50\t60\t.\t-\t.\tID=g3\n')


def test_colours_as_the_page_has_them():
    assert circos_plot.rgb('#1f77b4', 'gff') == (31, 119, 180)
    assert circos_plot.rgb('#abc', 'gff') == (170, 187, 204)
    assert circos_plot.rgb('rgb(214, 39, 40)', 'bigwig') == (214, 39, 40)
    assert circos_plot.rgb('nonsense; inject', 'bigwig') == circos_plot.DEFAULT_COLOURS['bigwig']


def test_bins_and_chromosomes():
    assert [circos_plot.nice(v) for v in (0.5, 3, 7, 140_000, 260_000)] == [1, 5, 10, 200_000, 500_000]
    sizes = {'Chr1': 30_000_000, 'Chr2': 20_000_000, 'ChrM': 300_000}
    assert [c for c, _ in circos_plot.choose_chromosomes(sizes, None, True)] == ['Chr1', 'Chr2']
    assert [c for c, _ in circos_plot.choose_chromosomes(sizes, None, False)] == ['Chr1', 'Chr2', 'ChrM']
    assert circos_plot.choose_chromosomes(sizes, 'ChrM', True) == [('ChrM', 300_000)]
    with pytest.raises(circos_plot.CircosError):
        circos_plot.choose_chromosomes(sizes, 'Chr9', True)


def test_genes_per_bin_from_a_gff(tmp_path):
    path = tmp_path / 'genes.gff'
    path.write_text(GFF)
    rows, what = circos_plot.feature_bins(str(path), 'gff', [('Chr1', 0, 2000), ('Chr2', 0, 1000)], 1000)
    assert what == 'genes'                                    # the mRNA is not counted
    assert rows == [('Chr1', 0, 1000, 1.0), ('Chr1', 1000, 2000, 1.0), ('Chr2', 0, 1000, 1.0)]


def test_the_circos_files(tmp_path):
    chroms = [('Chr1', 2000), ('Chr 2', 1000)]
    assert circos_plot.karyotype(chroms).splitlines() == ['chr - Chr1 Chr1 0 2000 lgrey', 'chr - Chr_2 Chr_2 0 1000 vlgrey']
    assert circos_plot.data_file([('Chr 2', 0, 1000, 2.5)]) == 'Chr_2 0 999 2.5\n'
    conf = circos_plot.conf(chroms, [{'file': 'data/track-1.txt', 'rgb': (1, 2, 3)},
                                     {'file': 'data/track-2.txt', 'rgb': (4, 5, 6)}])
    assert conf.count('<plot>') == 2 and 'fill_color  = 1,2,3' in conf and 'karyotype = karyotype.txt' in conf
    svg = circos_plot.legend_svg('<svg></svg>', 'A & B', [{'name': '<x>', 'what': 'genes', 'rgb': (1, 2, 3)}])
    assert '&lt;x&gt;' in svg and 'A &amp; B' in svg and svg.endswith('</svg>')


@pytest.fixture
def world(tmp_path, monkeypatch):
    data_root = tmp_path / 'genome-data'

    class Config(app_module.Config):
        DATA_ROOT = data_root

    application = app_module.create_app(Config())
    application.config['TESTING'] = True
    genome, track, other = uuid.uuid4().hex, uuid.uuid4().hex, uuid.uuid4().hex
    (data_root / 'tracks').mkdir(parents=True, exist_ok=True)
    gff = data_root / 'tracks' / 'genes.gff'
    gff.write_text(GFF)
    conn = db.connect(data_root / 'genome.db')
    conn.execute('INSERT INTO reference_genomes (id, name, assembly_id, owner_lab, genome_visibility, created_by) '
                 "VALUES (?, 'Plant', 'P1', 'labA', 'owner_lab_only', 'alice')", (genome,))
    conn.executemany('INSERT INTO chrom_sizes (genome_id, chrom, length_bp) VALUES (?, ?, ?)',
                     [(genome, 'Chr1', 2000), (genome, 'Chr2', 1000)])
    conn.executemany('INSERT INTO tracks (id, genome_id, name, type, owner_lab, visibility, track_status, file_path, '
                     "created_by) VALUES (?, ?, ?, 'gff', ?, 'lab_only', 'ready', ?, 'alice')",
                     [(track, genome, 'Genes', 'labA', str(gff)), (other, genome, 'Theirs', 'labB', str(gff))])
    conn.commit()
    conn.close()
    members = {'alice': 'labA', 'bob': 'labB'}
    monkeypatch.setattr(identity, 'is_admin', lambda u: False)
    monkeypatch.setattr(labs_reader, 'labs_of', lambda u: [
        {'lab_id': members[u], 'lab': members[u], 'role': 'member'}] if u in members else [])
    monkeypatch.setattr(labs_reader, 'chair_group_of', lambda lab: None)
    monkeypatch.setattr(labs_reader, 'labs_in_chair_group', lambda g: [])
    monkeypatch.setattr(labs_reader, 'demo_lab_ids', lambda: set())
    client = application.test_client()

    def draw(user, tracks, **extra):
        return client.post(f'/api/genomes/{genome}/circos', headers={'X-Suite-User': user},
                           json={'tracks': [{'id': t, 'color': '#336699'} for t in tracks], **extra})
    return {'draw': draw, 'track': track, 'other': other}


def test_only_what_they_may_see_is_drawn(world, monkeypatch):
    seen = {}
    monkeypatch.setattr(circos_plot, 'make', lambda genome, sizes, tracks, **kw: seen.update(tracks=tracks, **kw) or {'ok': 1})
    assert world['draw']('bob', [world['track']]).status_code == 404            # not their genome
    assert world['draw']('alice', [world['track'], world['other']], chrom='Chr2', resolution='fine').status_code == 200
    assert [t['name'] for t in seen['tracks']] == ['Genes']                      # another lab's track is left out
    assert seen['tracks'][0]['color'] == '#336699' and seen['chrom'] == 'Chr2' and seen['resolution'] == 'fine'


@pytest.mark.skipif(not circos_plot.available(), reason='Circos is not installed here (it is in the image)')
def test_a_real_circos_plot(world):
    res = world['draw']('alice', [world['track']], skip_small=False)
    assert res.status_code == 200, res.get_json()
    out = res.get_json()
    assert out['svg'].lstrip().startswith('<?xml') or '<svg' in out['svg'][:500]
    assert 'trellis-legend' in out['svg'] and out['rings'][0]['what'] == 'genes per step'
    names = zipfile.ZipFile(io.BytesIO(base64.b64decode(out['bundle']))).namelist()
    assert {'circos/circos.conf', 'circos/karyotype.txt', 'circos/data/track-1.txt', 'circos/README.txt'} <= set(names)


def test_what_to_draw():
    sizes = {'Chr1': 30_000_000, 'Chr2': 20_000_000, 'ChrM': 300_000}
    assert circos_plot.choose_spans(sizes, ['Chr2', 'Chr1']) == [('Chr1', 0, 30_000_000), ('Chr2', 0, 20_000_000)]
    assert circos_plot.choose_spans(sizes, region='chr1:1,000,001-5,000,000') == [('Chr1', 1_000_000, 5_000_000)]
    for bad in ('Chr9:1-5000', 'Chr1', 'Chr1:1-10'):
        with pytest.raises(circos_plot.CircosError):
            circos_plot.choose_spans(sizes, region=bad)
    spans = [('Chr1', 1_000_000, 5_000_000)]
    assert circos_plot.step_size(spans, 'medium', None) == 5000
    assert circos_plot.step_size(spans, 'medium', 100_000) == 100_000
    assert circos_plot.step_size(spans, 'medium', 10) == 800          # never more than 5,000 steps
    assert circos_plot.where_text(spans, sizes, False) == 'Chr1:1,000,001-5,000,000'


def test_links_both_ends_drawn_strongest_first(tmp_path):
    path = tmp_path / 'links.bedpe'
    path.write_text('#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\n'
                    'Chr1\t100\t200\tChr2\t300\t400\tweak\t1\n'
                    'Chr1\t500\t600\tChr2\t700\t800\tstrong\t9\n'
                    'Chr1\t900\t950\tChrM\t10\t20\toff\t50\n'
                    'Chr1\tx\t950\tChr2\t10\t20\tbad\t50\n')
    pairs, total = circos_plot.links(str(path), [('Chr1', 0, 1000), ('Chr2', 0, 1000)], limit=1)
    assert total == 2 and pairs == [(('Chr1', 500, 600), ('Chr2', 700, 800))]
    assert circos_plot.links_file(pairs) == 'Chr1 500 599 Chr2 700 799\n'


def test_gene_labels_by_name(tmp_path):
    path = tmp_path / 'genes.gff'
    path.write_text(GFF.replace('ID=g1', 'ID=g1;Name=NAC001'))
    found, missing = circos_plot.gene_labels([(str(path), 'gff')], ['nac001', 'nothere'], [('Chr1', 0, 2000)])
    assert found == [('Chr1', 99, 900, 'NAC001')] and missing == ['nothere']


def test_a_plot_with_every_option_is_configured(tmp_path):
    rings = [{'file': 'data/track-1.txt', 'rgb': (10, 20, 30), 'style': 'heatmap', 'height': 1},
             {'file': 'data/track-2.txt', 'rgb': (1, 2, 3), 'style': 'line', 'height': 3, 'max': 7}]
    text = circos_plot.conf([('Chr1', 1_000_000, 5_000_000)], rings, [{'file': 'data/links-3.txt', 'rgb': (5, 6, 7)}],
                            'data/labels.txt', 'data/highlight.txt', {'Chr1': 30_000_000})
    assert 'chromosomes = Chr1:1000.000000-5000.000000' in text and 'chromosomes_display_default = no' in text
    assert 'type        = heatmap' in text and 'ring1_8 = 10,20,30' in text
    assert 'type        = line' in text and 'max         = 7' in text
    assert '<links>' in text and 'color         = 5,6,7,0.5' in text
    assert 'type          = text' in text and 'type       = highlight' in text


@pytest.mark.skipif(not circos_plot.available(), reason='Circos is not installed here (it is in the image)')
def test_a_real_circos_plot_with_everything(world, tmp_path):
    res = world['draw']('alice', [world['track']], region='Chr1:1-1500', step=200, labels='g1, nope',
                        highlight={'chrom': 'Chr1', 'start': 100, 'end': 900}, title='Mine', log_scale=True)
    assert res.status_code == 200, res.get_json()
    out = res.get_json()
    assert out['where'] == 'Chr1:1-1,500' and out['bin_size'] == 200 and 'Mine' in out['svg']
    assert any('nope' in n for n in out['notes'])
