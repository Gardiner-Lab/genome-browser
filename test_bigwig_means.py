"""Means over a region count the bases a coverage track has no data for as zero."""
import numpy as np
import pytest

pyBigWig = pytest.importorskip('pyBigWig')

import bigwig_means


@pytest.fixture
def sparse(tmp_path):
    """10 over chr1:0-2000 and 5 over chr1:50000-60000; nothing else. chr2 has no reads at all."""
    path = tmp_path / 'sparse.bw'
    bw = pyBigWig.open(str(path), 'w')
    bw.addHeader([('chr1', 100_000), ('chr2', 30_000)])
    bw.addEntries(['chr1', 'chr1'], [0, 50_000], ends=[2000, 60_000], values=[10.0, 5.0])
    bw.close()
    return path


def test_mean_counts_missing_bases_as_zero(sparse):
    bw = pyBigWig.open(str(sparse))
    try:
        assert bw.stats('chr1', 0, 4000, type='mean')[0] == pytest.approx(10.0)   # what pyBigWig says
        assert bigwig_means.mean(bw, 'chr1', 0, 4000) == pytest.approx(5.0)      # what it is
        assert bigwig_means.mean(bw, 'chr1', 20_000, 24_000) == 0.0               # no reads: zero, not None
        assert bigwig_means.mean(bw, 'chr2', 0, 1000) == 0.0
        assert bigwig_means.mean(bw, 'chr1', 99_000, 200_000) == 0.0             # clipped to the chromosome
        assert bigwig_means.mean(bw, 'chrX', 0, 1000) is None                    # not in the track at all
    finally:
        bw.close()


def test_bin_means_with_a_shorter_last_bin(sparse):
    bw = pyBigWig.open(str(sparse))
    try:
        assert bigwig_means.bin_means(bw, 'chr1', 0, 100_000, 10_000) == pytest.approx(
            [2.0, 0, 0, 0, 0, 5.0, 0, 0, 0, 0])
        assert bigwig_means.bin_means(bw, 'chr1', 0, 3000, 2500) == pytest.approx([8.0, 0.0])
        assert bigwig_means.bin_means(bw, 'chrX', 0, 3000, 1000) == []
    finally:
        bw.close()


def test_features_table(sparse, tmp_path):
    bed = tmp_path / 'regions.bed'
    bed.write_text('chr1\t0\t4000\tGENE1\t0\t+\nchr1\t20000\t24000\tEMPTY\t0\t-\nchrX\t0\t10\tELSEWHERE\t0\t+\n')
    out = tmp_path / 'raw.tab'
    bigwig_means.main(['features', str(bed), str(sparse), '--labels', 'A', '-o', str(out)])
    assert out.read_text().splitlines() == ["#'chr'\t'start'\t'end'\t'A'", 'chr1\t0\t4000\t5',
                                            'chr1\t20000\t24000\t0', 'chrX\t0\t10\tnan']


def test_bins_keep_the_bins_where_only_one_track_has_reads(sparse, tmp_path):
    other = tmp_path / 'other.bw'
    bw = pyBigWig.open(str(other), 'w')
    bw.addHeader([('chr1', 100_000)])
    bw.addEntries(['chr1'], [0], ends=[1000], values=[4.0])
    bw.close()
    npz, raw = tmp_path / 's.npz', tmp_path / 'bins.tab'
    bigwig_means.main(['bins', str(sparse), str(other), '--labels', 'A', 'B', '-bs', '10000',
                       '-o', str(npz), '--outRawCounts', str(raw)])
    z = np.load(npz)
    assert list(z['labels']) == ['A', 'B']
    assert z['matrix'].shape == (10, 2)          # chr1 only: chr2 is not in the other track
    assert not np.isnan(z['matrix']).any()
    assert z['matrix'][0] == pytest.approx([2.0, 0.4])
    assert z['matrix'][5] == pytest.approx([5.0, 0.0])   # deepTools left this NaN, and plotCorrelation dropped it
    assert raw.read_text().splitlines()[6] == 'chr1\t50000\t60000\t5\t0'
