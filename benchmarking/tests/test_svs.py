"""evaluate_svs: how each non-SNP sits relative to its flanking SNPs."""
import numpy as np

from evaluate_phasing import (allele_length, evaluate_svs, get_overlapping_sites,
                              load_phasing, nearest_snp_indices)


def _snp(chrom, pos, gt, ps=1):
    return (chrom, pos, 'ACGT'[pos % 4], 'TGCA'[pos % 4], gt, ps)


def _ins(chrom, pos, gt, ps=1, length=3):
    return (chrom, pos, 'A', 'A' + 'T' * length, gt, ps)


def _del(chrom, pos, gt, ps=1, length=3):
    return (chrom, pos, 'A' + 'T' * length, 'A', gt, ps)


def _evaluate(write_vcf, target_records, truth_records, **kwargs):
    target = load_phasing(write_vcf('t.vcf', target_records, contigs=('chr1', 'chr2')))
    truth = load_phasing(write_vcf('g.vcf', truth_records, contigs=('chr1', 'chr2')))
    pairs = get_overlapping_sites(target, truth)
    return evaluate_svs(pairs, target, truth, **kwargs)


def _check_identities(res):
    assert res.correct + res.flipped + res.ambiguous + res.one_sided + res.no_anchor == res.svs
    assert res.connections == 2 * (res.correct + res.flipped + res.ambiguous) + res.one_sided


def test_allele_length():
    class V:
        ref = np.array(['AT'], dtype=object)
        alt = np.empty(1, dtype=object)
    V.alt[0] = ('A', 'ATTTT')
    assert allele_length(V, 0) == 5


def test_nearest_snp_indices_skips_non_snps():
    class T:
        variant_type = np.array(['SNP', 'INS', 'DEL', 'SNP'], dtype=object)
    chain = [(i, i) for i in range(4)]
    left, right = nearest_snp_indices(chain, T)
    assert left == [None, 0, 0, 0]
    assert right == [3, 3, 3, None]


def test_sv_agreeing_with_both_anchors_is_correct(write_vcf):
    t = [_snp('chr1', 10, '0|1'), _ins('chr1', 20, '0|1'), _snp('chr1', 30, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _ins('chr1', 20, '0|1', None),
         _snp('chr1', 30, '0|1', None)]
    res = _evaluate(write_vcf, t, g)
    assert (res.total, res.svs, res.connections, res.errors) == (1, 1, 2, 0)
    assert (res.correct, res.flipped, res.ambiguous) == (1, 0, 0)
    _check_identities(res)


def test_sv_out_of_phase_with_both_anchors_is_flipped(write_vcf):
    t = [_snp('chr1', 10, '0|1'), _ins('chr1', 20, '1|0'), _snp('chr1', 30, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _ins('chr1', 20, '0|1', None),
         _snp('chr1', 30, '0|1', None)]
    res = _evaluate(write_vcf, t, g)
    assert (res.connections, res.errors) == (2, 2)     # both sides wrong
    assert (res.correct, res.flipped, res.ambiguous) == (0, 1, 0)
    _check_identities(res)


def test_switch_between_the_anchors_is_ambiguous(write_vcf):
    """The anchors disagree with each other, so the SV cannot be blamed."""
    t = [_snp('chr1', 10, '0|1'), _ins('chr1', 20, '0|1'), _snp('chr1', 30, '1|0')]
    g = [_snp('chr1', 10, '0|1', None), _ins('chr1', 20, '0|1', None),
         _snp('chr1', 30, '0|1', None)]
    res = _evaluate(write_vcf, t, g)
    assert (res.connections, res.errors) == (2, 1)     # exactly one side wrong either way
    assert (res.correct, res.flipped, res.ambiguous) == (0, 0, 1)
    _check_identities(res)


def test_anchor_in_another_phase_set_leaves_one_side(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _ins('chr1', 20, '0|1', 2), _snp('chr1', 30, '0|1', 2)]
    g = [_snp('chr1', 10, '0|1', None), _ins('chr1', 20, '0|1', None),
         _snp('chr1', 30, '0|1', None)]
    res = _evaluate(write_vcf, t, g)
    assert (res.one_sided, res.connections, res.errors) == (1, 1, 0)
    assert res.no_anchor == 0
    _check_identities(res)


def test_sv_with_no_snp_anchor_at_all(write_vcf):
    t = [_ins('chr1', 10, '0|1'), _del('chr1', 20, '0|1')]
    g = [_ins('chr1', 10, '0|1', None), _del('chr1', 20, '0|1', None)]
    res = _evaluate(write_vcf, t, g)
    assert (res.svs, res.no_anchor, res.connections) == (2, 2, 0)
    _check_identities(res)


def test_adjacent_non_snps_share_the_same_snp_anchors(write_vcf):
    """A run of non-SNPs is anchored on the SNPs either side, not on each other."""
    t = [_snp('chr1', 10, '0|1'), _ins('chr1', 20, '0|1'),
         _del('chr1', 30, '1|0'), _snp('chr1', 40, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _ins('chr1', 20, '0|1', None),
         _del('chr1', 30, '0|1', None), _snp('chr1', 40, '0|1', None)]
    res = _evaluate(write_vcf, t, g)
    assert res.svs == 2 and res.connections == 4
    assert (res.correct, res.flipped) == (1, 1)        # judged independently
    assert res.errors == 2
    _check_identities(res)


def test_total_counts_unphased_non_snps(write_vcf):
    """`total` is a property of the call set, not of how much got phased."""
    t = [_snp('chr1', 10, '0|1'), _ins('chr1', 20, '0|1'),
         _del('chr1', 30, '0/1'), _snp('chr1', 40, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _ins('chr1', 20, '0|1', None),
         _del('chr1', 30, '0|1', None), _snp('chr1', 40, '0|1', None)]
    res = _evaluate(write_vcf, t, g)
    assert res.total == 2        # both non-SNPs
    assert res.svs == 1          # only the phased one is scorable


def test_min_sv_length_filters_by_longest_allele(write_vcf):
    t = [_snp('chr1', 10, '0|1'), _ins('chr1', 20, '0|1', length=2),
         _snp('chr1', 30, '0|1'), _ins('chr1', 40, '0|1', length=60), _snp('chr1', 50, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _ins('chr1', 20, '0|1', None, length=2),
         _snp('chr1', 30, '0|1', None), _ins('chr1', 40, '0|1', None, length=60),
         _snp('chr1', 50, '0|1', None)]
    assert _evaluate(write_vcf, t, g).svs == 2
    filtered = _evaluate(write_vcf, t, g, min_sv_length=50)
    assert (filtered.total, filtered.svs) == (1, 1)


def test_multi_and_other_types_are_counted(write_vcf):
    """Neither branch ever executed on the real call set."""
    t = [_snp('chr1', 10, '0|1'),
         ('chr1', 20, 'A', 'G,T', '1|2', 1),        # Multi
         ('chr1', 30, 'AT', 'GC', '0|1', 1),        # Other: equal length, not 1 bp
         _snp('chr1', 40, '0|1')]
    g = [_snp('chr1', 10, '0|1', None),
         ('chr1', 20, 'A', 'G,T', '1|2', None),
         ('chr1', 30, 'AT', 'GC', '0|1', None),
         _snp('chr1', 40, '0|1', None)]
    res = _evaluate(write_vcf, t, g)
    assert res.by_type['Multi'] == 1 and res.by_type['Other'] == 1
    assert res.svs == 2
    _check_identities(res)


def test_anchors_are_not_taken_across_a_chromosome_boundary(write_vcf):
    t = [_snp('chr1', 10, '0|1'), _ins('chr1', 20, '0|1'),
         _snp('chr2', 30, '0|1'), _ins('chr2', 40, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _ins('chr1', 20, '0|1', None),
         _snp('chr2', 30, '0|1', None), _ins('chr2', 40, '0|1', None)]
    res = _evaluate(write_vcf, t, g)
    # each SV has only its own chromosome's SNP to the left, none to the right
    assert res.svs == 2 and res.one_sided == 2 and res.connections == 2
    _check_identities(res)
