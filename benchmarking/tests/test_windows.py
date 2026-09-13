"""Hamming distance scored in fixed windows, the block-length-independent variant.

The block-wise Hamming distance of ``evaluate_phasing`` counts, per intersection
block, the variants that must move to the other haplotype.  One switch in the
middle of a block therefore costs about half the block, so a caller that builds
longer blocks scores worse at the same local error density.  These tests pin down
that the windowed statistic does not have that property.
"""
import pytest

from evaluate_phasing import (WINDOW_SECTIONS, evaluate_phasing, evaluate_windowed_hamming,
                              get_overlapping_sites, load_phasing, window_columns)


def _snp(chrom, pos, gt, ps=1):
    return (chrom, pos, 'ACGT'[pos % 4], 'TGCA'[pos % 4], gt, ps)


def _load(write_vcf, target_records, truth_records):
    target = load_phasing(write_vcf('t.vcf', target_records))
    truth = load_phasing(write_vcf('g.vcf', truth_records))
    return get_overlapping_sites(target, truth), target, truth


def _windowed(write_vcf, target_records, truth_records, window, unit='bp'):
    pairs, target, truth = _load(write_vcf, target_records, truth_records)
    return evaluate_windowed_hamming(pairs, target, truth, window, unit)


# One variant per kb throughout, so a 10 kb window holds exactly ten of them.
def _spread(n, ps_every=None, flip_from=None, start=1000, step=1000):
    """Target records: ``n`` variants a fixed distance apart.

    ``ps_every`` starts a new phase set every that many variants; ``flip_from``
    inverts the genotype from that index on, which plants one switch error.
    """
    out = []
    for k in range(n):
        ps = 1 if ps_every is None else 1 + k // ps_every
        gt = '1|0' if flip_from is not None and k >= flip_from else '0|1'
        out.append(_snp('chr1', start + k * step, gt, ps))
    return out


def _truth(n, start=1000, step=1000):
    return [_snp('chr1', start + k * step, '0|1', None) for k in range(n)]


# --------------------------------------------------------------------------
# the property the metric exists for
# --------------------------------------------------------------------------

def test_block_length_changes_the_blockwise_rate_but_not_the_windowed_one(write_vcf):
    """Same errors, same positions, different block lengths.

    Sixty variants a kb apart with one switch planted at the 33rd.  Phased as one
    block the block-wise Hamming distance charges the whole tail to it; phased as
    three 20 kb blocks it charges only the block the switch falls in.  The 10 kb
    windows are the same genomic windows either way, so they charge the same.
    """
    truth = _truth(60)
    long_block = _spread(60, flip_from=33)
    short_blocks = _spread(60, ps_every=20, flip_from=33)

    long_pairs, long_t, g = _load(write_vcf, long_block, truth)
    short_pairs, short_t, _ = _load(write_vcf, short_blocks, truth)

    long_within = evaluate_phasing(long_pairs, long_t, g)
    short_within = evaluate_phasing(short_pairs, short_t, g)
    assert long_within.hamming == 27            # min(33 agreeing, 27 disagreeing)
    assert short_within.hamming == 7            # only the middle block is mixed
    assert long_within.hamming / long_within.covered_variants == pytest.approx(0.45)
    assert short_within.hamming / short_within.covered_variants == pytest.approx(7 / 60)

    long_win = evaluate_windowed_hamming(long_pairs, long_t, g, 10000)
    short_win = evaluate_windowed_hamming(short_pairs, short_t, g, 10000)
    assert long_win.hamming == short_win.hamming == 3      # the one mixed window
    assert long_win.covered_variants == short_win.covered_variants == 60
    assert long_win.hamming_rate == short_win.hamming_rate == pytest.approx(0.05)


def test_one_switch_costs_half_the_block_but_only_half_a_window(write_vcf):
    """The mechanism itself: a lone switch scales with what it sits in."""
    truth = _truth(100)
    target = _spread(100, flip_from=62)
    pairs, t, g = _load(write_vcf, target, truth)

    assert evaluate_phasing(pairs, t, g).hamming == 38      # the whole tail of the block
    # only the one window the switch falls in is charged, and only for its short side
    assert evaluate_windowed_hamming(pairs, t, g, 10000).hamming == 2
    assert evaluate_windowed_hamming(pairs, t, g, 25000).hamming == 12


def test_windowed_hamming_never_exceeds_the_blockwise_one(write_vcf):
    """Windows only ever subdivide the blocks, and min() is subadditive."""
    truth = _truth(80)
    target = _spread(80, flip_from=17)
    pairs, t, g = _load(write_vcf, target, truth)
    blockwise = evaluate_phasing(pairs, t, g).hamming
    for window in (5000, 10000, 20000, 40000, 80000, 1000000):
        assert evaluate_windowed_hamming(pairs, t, g, window).hamming <= blockwise


def test_a_window_larger_than_the_blocks_reproduces_the_blockwise_statistic(write_vcf):
    """Why the window has to be chosen below the block N50."""
    truth = _truth(40)
    target = _spread(40, ps_every=20, flip_from=25)
    pairs, t, g = _load(write_vcf, target, truth)
    within = evaluate_phasing(pairs, t, g)
    huge = evaluate_windowed_hamming(pairs, t, g, 10 ** 9)
    assert (huge.hamming, huge.covered_variants) == (within.hamming, within.covered_variants)
    assert huge.windows == within.intersection_blocks


# --------------------------------------------------------------------------
# window construction
# --------------------------------------------------------------------------

def test_windows_do_not_cross_a_block_boundary(write_vcf):
    """Orientation between two blocks is undefined, so the window is cut there."""
    # ten variants inside one 10 kb window; the second half is a separate,
    # inverted phase set
    target = [_snp('chr1', 1000 + k * 100, '0|1' if k < 5 else '1|0', 1 if k < 5 else 2)
              for k in range(10)]
    truth = [_snp('chr1', 1000 + k * 100, '0|1', None) for k in range(10)]
    result = _windowed(write_vcf, target, truth, 10000)
    assert result.windows == 2               # not one window of ten
    assert result.hamming == 0               # each block is internally correct


def test_windows_are_anchored_on_absolute_coordinates(write_vcf):
    """Two call sets are scored on the same genomic windows, whatever they phase.

    The second call set leaves the first two variants unphased, which would shift
    a block-anchored tiling by two variants; with absolute anchoring the windows
    it shares with the first call set are unchanged.
    """
    truth = _truth(20)
    full = _spread(20)
    partial = [_snp('chr1', 1000 + k * 1000, '0/1' if k < 2 else '0|1', None if k < 2 else 1)
               for k in range(20)]
    a = _windowed(write_vcf, full, truth, 10000)
    b = _windowed(write_vcf, partial, truth, 10000)
    assert a.windows == 2 and a.covered_variants == 20
    assert b.windows == 2 and b.covered_variants == 18   # same two windows, two fewer variants


def test_groups_with_a_single_variant_are_skipped(write_vcf):
    """As single-variant blocks are: they can never disagree, so they only dilute."""
    truth = _truth(11)
    target = _spread(11)
    result = _windowed(write_vcf, target, truth, 10000)
    assert result.windows == 1               # the 11th variant is alone in its window
    assert result.covered_variants == 10


def test_variant_windows_cut_into_runs_of_fixed_length(write_vcf):
    truth = _truth(25)
    target = _spread(25)
    result = _windowed(write_vcf, target, truth, 10, unit='variants')
    assert (result.windows, result.covered_variants) == (3, 25)   # 10 + 10 + 5


def test_variant_windows_ignore_the_spacing_bp_windows_depend_on(write_vcf):
    """The two units answer different questions; this is the one they differ on."""
    positions = (1000, 2000, 3000, 4000, 501000, 502000)     # two clusters far apart
    truth = [_snp('chr1', p, '0|1', None) for p in positions]
    target = [_snp('chr1', p, '0|1', 1) for p in positions]
    pairs, t, g = _load(write_vcf, target, truth)
    assert evaluate_windowed_hamming(pairs, t, g, 10000).windows == 2      # two clusters
    assert evaluate_windowed_hamming(pairs, t, g, 6, unit='variants').windows == 1


def test_an_unknown_unit_is_rejected(write_vcf):
    pairs, t, g = _load(write_vcf, _spread(4), _truth(4))
    with pytest.raises(ValueError):
        evaluate_windowed_hamming(pairs, t, g, 10000, unit='kb')


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def test_window_columns_shape_is_the_same_with_and_without_a_result(write_vcf):
    result = _windowed(write_vcf, _spread(20, flip_from=15), _truth(20), 10000)
    filled = window_columns('winbp', result)
    empty = window_columns('winbp', None)
    assert filled.keys() == empty.keys()
    assert all(v == '' for v in empty.values())
    assert filled['winbp_window'] == 10000
    assert filled['winbp_hamming_rate'] == result.hamming / result.covered_variants


def test_both_window_sections_reach_the_summary_tsv(run_cli, write_vcf, tmp_path):
    target = write_vcf('t.vcf', _spread(30, flip_from=22))
    truth = write_vcf('g.vcf', _truth(30))
    out = tmp_path / 'summary.tsv'
    stdout = run_cli('--vcf', target, '--gt_vcf', truth, '--hamming_window', '10000',
                     '--hamming_window_variants', '10', '--summary_tsv', out).stdout
    header, values = (line.split('\t') for line in out.read_text().rstrip('\n').split('\n'))
    row = dict(zip(header, values))
    for prefix, _, heading in WINDOW_SECTIONS:
        assert heading in stdout
        assert row[f'{prefix}_window'] != ''
        assert row[f'{prefix}_hamming'] != ''
    assert row['winbp_window'] == '10000' and row['winvar_window'] == '10'
    assert row['within_hamming'] != ''       # the existing statistics are untouched


def test_the_window_tests_can_be_turned_off_without_moving_the_header(run_cli, write_vcf,
                                                                     tmp_path):
    target = write_vcf('t.vcf', _spread(30, flip_from=22))
    truth = write_vcf('g.vcf', _truth(30))
    on, off = tmp_path / 'on.tsv', tmp_path / 'off.tsv'
    run_cli('--vcf', target, '--gt_vcf', truth, '--summary_tsv', on)
    stdout = run_cli('--vcf', target, '--gt_vcf', truth, '--hamming_window', '0',
                     '--hamming_window_variants', '0', '--summary_tsv', off).stdout
    assert on.read_text().split('\n')[0] == off.read_text().split('\n')[0]
    row = dict(zip(off.read_text().split('\n')[0].split('\t'),
                   off.read_text().split('\n')[1].split('\t')))
    assert row['winbp_window'] == '' and row['winvar_hamming_rate'] == ''
    for _, _, heading in WINDOW_SECTIONS:
        assert heading not in stdout
