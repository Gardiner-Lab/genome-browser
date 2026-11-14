"""The search box finds a feature by its name in the genome's own annotation, not at igv.org."""
from pathlib import Path

import feature_index
import regions
from test_analysis import world  # noqa: F401  (the same genome, tracks and members)


def test_a_gene_is_found_by_its_name_in_the_genomes_own_annotation(world):
    w = world
    path = f"/api/genomes/{w['genome']}/feature-search"
    found = w['call']('GET', path + '?name=GENE2').get_json()
    assert found[0] == {'chromosome': 'Chr1', 'start': 6000, 'end': 9000, 'name': 'GENE2', 'type': 'gene',
                        'track': 'Genes'}
    assert w['call']('GET', path + '?name=te1').get_json()[0]['type'] == 'transposable_element'   # any case
    assert w['call']('GET', path + '?name=NOPE').get_json() == []
    assert w['call']('GET', path + '?name=').get_json() == []
    # bob's lab cannot see alice's lab's genome: nothing about it, not even that the gene is there.
    assert w['call']('GET', path + '?name=GENE2', 'bob').status_code == 404


def test_genes_come_before_their_parts_and_old_indexes_get_the_name_index(tmp_path):
    gff = tmp_path / 'a.gff'
    gff.write_text('##gff-version 3\n'
                   'Chr3\tA\texon\t100\t200\t.\t-\t.\tID=x;Name=AT3G01835\n'
                   'Chr3\tA\tgene\t100\t900\t.\t-\t.\tID=AT3G01835;Name=AT3G01835\n'
                   'Chr3\tA\tlnc_RNA\t100\t900\t.\t-\t.\tID=AT3G01835.1;Name=AT3G01835.1;Parent=AT3G01835\n')
    parse = regions.feature_parser('gff')
    idx = feature_index.build(gff, parse)
    import sqlite3
    conn = sqlite3.connect(str(idx))
    conn.execute('DROP INDEX features_by_name')       # as an index built before searching existed
    conn.commit()
    conn.close()
    found = feature_index.named(gff, parse, 'at3g01835')
    assert [(f[4], f[1], f[2]) for f in found] == [('gene', 99, 900), ('exon', 99, 200)]
    conn = sqlite3.connect(str(idx))
    assert conn.execute("SELECT 1 FROM sqlite_master WHERE name = 'features_by_name'").fetchone()
    conn.close()
    assert feature_index.named(tmp_path / 'missing.gff', parse, 'x') == []


def test_the_browser_asks_this_server_not_igv_org():
    js = (Path(__file__).with_name('static') / 'browser.js').read_text(encoding='utf-8')
    assert "/feature-search?name=$FEATURE$" in js and 'search: {' in js


def test_a_gene_is_found_by_its_symbol_alias_or_full_name(tmp_path):
    gff = tmp_path / 'b.gff'
    gff.write_text('##gff-version 3\n'
                   'Chr5\tA\tgene\t3173382\t3179448\t.\t-\t.\tID=AT5G10140;Name=AT5G10140;'
                   'full_name=FLOWERING LOCUS C;symbol=FLC;Alias=FLF,AGL25\n'
                   'Chr5\tA\tmRNA\t3173382\t3179448\t.\t-\t.\tID=AT5G10140.1;Name=AT5G10140.1;Parent=AT5G10140;symbol=FLC\n'
                   'Chr1\tA\tgene\t100\t200\t.\t+\t.\tID=AT1G00001;Name=AT1G00001;symbol=AT1G00001\n')
    parse = regions.feature_parser('gff')
    for query in ('FLC', 'flc', 'FLOWERING LOCUS C', 'FLF', 'agl25', 'AT5G10140'):
        found = feature_index.named(gff, parse, query)
        assert [(f[0], f[1], f[2], f[3], f[4]) for f in found] == [('Chr5', 3173381, 3179448, 'AT5G10140', 'gene')], query
    assert feature_index.named(gff, parse, 'AT1G00001') == [('Chr1', 99, 200, 'AT1G00001', 'gene', 1)]   # once
    assert feature_index.named(gff, parse, 'FLOWERING') == []
    feature_index.prepare_search(gff, parse, 'gff')         # again: already made, nothing to do
