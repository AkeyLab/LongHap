"""evaluate_junctions: the connections between phase blocks."""
from evaluate_phasing import (evaluate_junctions, evaluate_phasing, get_overlapping_sites,
                              load_phasing)


def _snp(chrom, pos, gt, ps=1):
    return (chrom, pos, 'ACGT'[pos % 4], 'TGCA'[pos % 4], gt, ps)


def _load(write_vcf, target_records, truth_records):
    target = load_phasing(write_vcf('t.vcf', target_records, contigs=('chr1', 'chr2')))
    truth = load_phasing(write_vcf('g.vcf', truth_records, contigs=('chr1', 'chr2')))
    return get_overlapping_sites(target, truth), target, truth


def test_no_junction_inside_a_single_block(write_vcf):
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)]
    pairs, target, truth = _load(write_vcf, t, g)
    res = evaluate_junctions(pairs, target, truth)
    assert (res.junctions, res.errors) == (0, 0)


def test_junction_scored_at_a_block_boundary(write_vcf):
    """Two blocks, joined in the correct relative orientation."""
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '0|1', 2), _snp('chr1', 40, '0|1', 2)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    pairs, target, truth = _load(write_vcf, t, g)
    res = evaluate_junctions(pairs, target, truth)
    assert (res.junctions, res.errors) == (1, 0)


def test_junction_error_when_the_second_block_is_inverted(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '1|0', 2), _snp('chr1', 40, '1|0', 2)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    pairs, target, truth = _load(write_vcf, t, g)
    res = evaluate_junctions(pairs, target, truth)
    assert (res.junctions, res.errors) == (1, 1)
    assert res.positions == [('chr1', 20, 30)]


def test_junctions_and_within_block_pairs_partition_the_chain(write_vcf):
    """The identity relied on when reporting the two statistics side by side."""
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '1|0', 1),
         _snp('chr1', 30, '0|1', 2), _snp('chr1', 40, '0|1', 2),
         _snp('chr1', 50, '1|0', 3), _snp('chr1', 60, '0|1', 3)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40, 50, 60)]
    pairs, target, truth = _load(write_vcf, t, g)

    within = evaluate_phasing(pairs, target, truth)
    junctions = evaluate_junctions(pairs, target, truth)
    chromosome_wide = evaluate_phasing(pairs, target, truth, ignore_phase_blocks=True)

    assert within.assessed_pairs + junctions.junctions == chromosome_wide.assessed_pairs
    assert within.switches + junctions.errors == chromosome_wide.switches
    within_iv = set(within.switch_positions)
    assert within_iv | set(junctions.positions) == set(chromosome_wide.switch_positions)
    assert within_iv & set(junctions.positions) == set()


def test_no_junction_across_a_chromosome_boundary(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr2', 30, '1|0', 2), _snp('chr2', 40, '1|0', 2)]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 20, '0|1', None),
         _snp('chr2', 30, '0|1', None), _snp('chr2', 40, '0|1', None)]
    pairs, target, truth = _load(write_vcf, t, g)
    res = evaluate_junctions(pairs, target, truth)
    assert (res.junctions, res.errors) == (0, 0)


def test_a_singleton_target_block_is_counted_separately(write_vcf):
    """A call set that strands a variant at its old phase set fragments the chain.

    Site 30 sits in a phase set of its own, so what would be one junction between
    two blocks is scanned as two -- and because the stranded site keeps the
    orientation it had before the merge, one of the two is charged whichever way
    the merge went.  The pairs are still scored, so the partition holds; they are
    counted a second time so that the contamination is visible.
    """
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '0|1', 2),                   # stranded at the old PS
         _snp('chr1', 40, '1|0', 1), _snp('chr1', 50, '1|0', 1)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40, 50)]
    pairs, target, truth = _load(write_vcf, t, g)
    res = evaluate_junctions(pairs, target, truth)
    assert res.junctions == 2                          # 20|30 and 30|40
    assert res.at_singleton == 2                       # both touch the stranded site
    assert (res.errors, res.errors_at_singleton) == (1, 1)


def test_junctions_between_intact_blocks_are_not_counted_as_singleton(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '1|0', 2), _snp('chr1', 40, '1|0', 2)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    pairs, target, truth = _load(write_vcf, t, g)
    res = evaluate_junctions(pairs, target, truth)
    assert (res.junctions, res.errors) == (1, 1)
    assert (res.at_singleton, res.errors_at_singleton) == (0, 0)


def test_a_stranded_variant_breaks_the_partition_identity(write_vcf):
    """Which is why the singleton junctions are counted separately.

    Site 30 sits in a phase set of its own *inside* the span of block 1, so the
    within-block pair 20|40 steps over it while the junction scan charges both
    20|30 and 30|40 for the same stretch.  The two statistics then overlap and no
    longer add up to the chromosome-wide one -- the identity
    ``test_junctions_and_within_block_pairs_partition_the_chain`` relies on holds
    only for call sets whose blocks do not interleave positionally.  ``at_singleton``
    is exactly the count of the junctions responsible.
    """
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '0|1', 2),
         _snp('chr1', 40, '1|0', 1), _snp('chr1', 50, '1|0', 1)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40, 50)]
    pairs, target, truth = _load(write_vcf, t, g)
    within = evaluate_phasing(pairs, target, truth)
    junctions = evaluate_junctions(pairs, target, truth)
    chromosome_wide = evaluate_phasing(pairs, target, truth, ignore_phase_blocks=True)

    overcount = within.assessed_pairs + junctions.junctions - chromosome_wide.assessed_pairs
    assert overcount == 1                      # the pair 20|40 counted on both sides
    assert junctions.at_singleton == 2
    # the split is a breakdown of the reported figures, never an exclusion
    assert junctions.at_singleton <= junctions.junctions
    assert junctions.errors_at_singleton <= junctions.errors
