"""The feature index answers region reads exactly as a full scan of the file would."""

import os
import random

import feature_index
import regions


def _brute(features, chrom, start, end):
    out = [{'chrom': c, 'start': s, 'end': e, 'name': n}
           for c, s, e, n in features
           if c == chrom and s < end and e > start]
    return sorted(out, key=lambda f: (f['start'], f['end']))


def test_index_matches_a_full_scan(tmp_path):
    rng = random.Random(7)
    lines = ['##gff-version 3\n']
    for i in range(3000):
        chrom = rng.choice(['Chr1', 'Chr2', 'ChrM'])
        start = rng.randint(1, 100_000)
        # Mostly short features, a few very long ones (genes spanning others).
        length = rng.randint(1, 50_000) if i % 97 == 0 else rng.randint(1, 900)
        lines.append(f'{chrom}\tsrc\tgene\t{start}\t{start + length}\t.\t+\t.\t'
                     f'ID=f{i};Name=feat{i}\n')
    path = tmp_path / 'genes.gff3'
    path.write_text(''.join(lines), encoding='utf-8')
    features = regions._scan_gff_features(''.join(lines))
    parse = regions.feature_parser('gff')

    for _ in range(200):
        chrom = rng.choice(['Chr1', 'Chr2', 'ChrM', 'ChrC'])
        start = rng.randint(0, 110_000)
        end = start + rng.randint(1, 20_000)
        assert feature_index.overlapping(path, parse, chrom, start, end) == \
            _brute(features, chrom, start, end)


def test_index_is_rebuilt_when_the_file_changes(tmp_path):
    path = tmp_path / 'f.bed'
    path.write_text('Chr1\t10\t20\ta\n', encoding='utf-8')
    parse = regions.feature_parser('bed')
    assert [f['name'] for f in
            feature_index.overlapping(path, parse, 'Chr1', 0, 100)] == ['a']
    path.write_text('Chr1\t30\t40\tb\n', encoding='utf-8')
    later = os.path.getmtime(feature_index.index_path(path)) + 5
    os.utime(path, (later, later))
    assert [f['name'] for f in
            feature_index.overlapping(path, parse, 'Chr1', 0, 100)] == ['b']


def test_crlf_files_give_clean_names(tmp_path):
    path = tmp_path / 'f.gff'
    path.write_bytes(b'Chr1\t.\tgene\t1\t9\t.\t+\t.\tID=g1;Name=ABC\r\n')
    found = feature_index.overlapping(path, regions.feature_parser('gff'),
                                      'Chr1', 0, 10)
    assert found == [{'chrom': 'Chr1', 'start': 0, 'end': 9, 'name': 'ABC'}]
