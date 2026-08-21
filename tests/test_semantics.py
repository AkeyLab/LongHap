"""Behaviour of the phasing model itself.

These tests are about answers, not crashes: does longhap recover the haplotype
it was given, does it break a block rather than guess, and does it stay put when
the input is perturbed in ways that should not matter.
"""
import numpy as np
import pytest

from synthetic import ReadSpec, alternating_snvs, ref_snv, simulate_read

POSITIONS = list(range(400, 3600, 200))
PATTERN = '0110100101101001'[:len(POSITIONS)]


def hap_string(lh):
    return ''.join('.' if h == -1 else str(h) for h in lh.haplotypes[0])


def truth_string(locus):
    return ''.join(str(v.hap[0]) for v in locus.variants)


def agrees_up_to_flip(got, want):
    flip = {'0': '1', '1': '0'}
    pairs = [(g, w) for g, w in zip(got, want) if g != '.']
    if not pairs:
        return True
    return all(g == w for g, w in pairs) or all(g == flip[w] for g, w in pairs)


def count_switches(got, want):
    """Switch errors: places where the relative phase flips between neighbours."""
    calls = [(g, w) for g, w in zip(got, want) if g != '.']
    rel = [int(g) ^ int(w) for g, w in calls]
    return sum(a != b for a, b in zip(rel, rel[1:]))


# --------------------------------------------------------------------------- #
# stability
# --------------------------------------------------------------------------- #
def test_phasing_is_deterministic(make_locus, run_longhap):
    locus = make_locus(lambda seq: alternating_snvs(seq, POSITIONS, PATTERN))
    a = run_longhap(locus, write=False)
    b = run_longhap(locus, write=False)
    assert np.array_equal(a.haplotypes, b.haplotypes)
    assert a.block_ends == b.block_ends


def test_swapping_the_truth_haplotypes_gives_the_mirror_answer(make_locus, run_longhap):
    """Relabelling which haplotype is which must not change the block structure."""
    flipped = ''.join('1' if c == '0' else '0' for c in PATTERN)
    a = run_longhap(make_locus(lambda s: alternating_snvs(s, POSITIONS, PATTERN)),
                    write=False)
    b = run_longhap(make_locus(lambda s: alternating_snvs(s, POSITIONS, flipped)),
                    write=False)
    assert a.block_ends == b.block_ends
    assert (a.haplotypes[0] == -1).tolist() == (b.haplotypes[0] == -1).tolist()


def test_extra_coverage_does_not_change_the_answer(make_locus, run_longhap):
    thin = [ReadSpec(s, s + 400, i % 2) for i, s in enumerate(range(0, 3600, 200))]
    thick = [ReadSpec(s, s + 400, i % 2) for i, s in enumerate(range(0, 3600, 50))]
    a = run_longhap(make_locus(lambda s: alternating_snvs(s, POSITIONS, PATTERN),
                               read_specs=thin), write=False)
    b = run_longhap(make_locus(lambda s: alternating_snvs(s, POSITIONS, PATTERN),
                               read_specs=thick), write=False)
    assert agrees_up_to_flip(hap_string(a), truth_string(
        make_locus(lambda s: alternating_snvs(s, POSITIONS, PATTERN))))
    assert count_switches(hap_string(b), PATTERN) == 0


# --------------------------------------------------------------------------- #
# noise tolerance
# --------------------------------------------------------------------------- #
def test_a_minority_of_wrong_reads_does_not_cause_a_switch(make_locus, run_longhap,
                                                           write_bam):
    """Two of twenty reads carry the opposite haplotype's alleles at one site."""
    holder = {}

    def build(seq):
        holder['seq'] = seq
        return alternating_snvs(seq, POSITIONS, PATTERN)

    locus = make_locus(build, read_specs=[])
    specs = [ReadSpec(s, s + 400, i % 2, name=f'r{i}')
             for i, s in enumerate(range(0, 3600, 100))]
    # flip two reads outright: they now look like the other haplotype everywhere
    specs[7].hap = 1 - specs[7].hap
    specs[13].hap = 1 - specs[13].hap
    locus.bam_path = write_bam(holder['seq'], specs, locus.variants, name='noisy.bam')

    lh = run_longhap(locus, write=False)
    assert count_switches(hap_string(lh), truth_string(locus)) == 0, \
        f'switch error: {hap_string(lh)} vs {truth_string(locus)}'


def test_a_site_seen_on_only_one_allele_is_not_confidently_phased(make_locus,
                                                                  run_longhap,
                                                                  write_bam):
    """Every read shows REF at one site, so there is no evidence to phase it."""
    holder = {}

    def build(seq):
        holder['seq'] = seq
        return alternating_snvs(seq, POSITIONS, PATTERN)

    locus = make_locus(build, read_specs=[])
    # all reads rendered from haplotype 0 across the middle site only
    specs = [ReadSpec(s, s + 400, i % 2, name=f'r{i}')
             for i, s in enumerate(range(0, 3600, 100))]
    variants = list(locus.variants)
    monomorphic = variants[8]
    monomorphic.hap = (monomorphic.hap[0], monomorphic.hap[0])
    locus.bam_path = write_bam(holder['seq'], specs, variants, name='mono.bam')

    lh = run_longhap(locus, write=False)
    assert lh.allele_coverage[:, 8].min() == 0, lh.allele_coverage[:, 8]


# --------------------------------------------------------------------------- #
# block breaks are preferred over guesses
# --------------------------------------------------------------------------- #
def test_an_unlinked_junction_breaks_rather_than_guesses(make_locus, run_longhap):
    positions = [400, 600, 800, 3000, 3200, 3400]
    reads = ([ReadSpec(300, 900, i % 2, name=f'l{i}') for i in range(12)] +
             [ReadSpec(2900, 3500, i % 2, name=f'r{i}') for i in range(12)])
    locus = make_locus(lambda s: alternating_snvs(s, positions, '010110'),
                       read_specs=reads, length=3600)
    lh = run_longhap(locus, write=False)
    assert len(lh.block_ends) >= 2, f'block_ends={lh.block_ends}'
    # within each block the phase must still be right
    got, want = hap_string(lh), truth_string(locus)
    assert agrees_up_to_flip(got[:3], want[:3]), (got[:3], want[:3])
    assert agrees_up_to_flip(got[3:], want[3:]), (got[3:], want[3:])


def test_phase_sets_change_across_a_block_boundary(make_locus, run_longhap,
                                                   read_phased_vcf):
    positions = [400, 600, 800, 3000, 3200, 3400]
    reads = ([ReadSpec(300, 900, i % 2, name=f'l{i}') for i in range(12)] +
             [ReadSpec(2900, 3500, i % 2, name=f'r{i}') for i in range(12)])
    locus = make_locus(lambda s: alternating_snvs(s, positions, '010110'),
                       read_specs=reads, length=3600)
    lh = run_longhap(locus)
    records = read_phased_vcf(lh.output_vcf)
    left = {ps for pos, gt, ps in records if pos <= 800 and '|' in gt}
    right = {ps for pos, gt, ps in records if pos >= 3000 and '|' in gt}
    assert left and right
    assert left != right, f'both blocks share phase set {left}'


# --------------------------------------------------------------------------- #
# the transition-index convention
# --------------------------------------------------------------------------- #
def test_transition_columns_are_indexed_by_left_variant(make_locus, run_longhap):
    """Column i of transition_matrix describes the step out of variant i.

    Nothing in the code states this, and several places depend on it, so it is
    pinned here: with contiguous phaseable variants, an informative junction
    between variant i and i+1 must show up in column i.
    """
    positions = [400, 600, 800]
    reads = [ReadSpec(300, 900, i % 2) for i in range(12)]
    locus = make_locus(lambda s: alternating_snvs(s, positions, '010'),
                       read_specs=reads, length=1200)
    lh = run_longhap(locus, phase=False)
    assert lh.transition_matrix.shape[2] == 2
    for i in range(2):
        col = lh.transition_matrix[:, :, i]
        assert col.max() > 0.5, f'column {i} is uninformative: {col}'


def test_haplotypes_array_is_indexed_by_variant_not_by_phaseable(make_locus,
                                                                 run_longhap):
    """haplotypes has one column per input variant, including unphaseable ones."""
    locus = make_locus(lambda seq: alternating_snvs(seq, POSITIONS, PATTERN))
    lh = run_longhap(locus, write=False)
    assert lh.haplotypes.shape == (2, lh.num_variants)
    assert lh.delta.shape == (2, lh.num_variants)


# --------------------------------------------------------------------------- #
# read haplotagging agrees with the variant phasing
# --------------------------------------------------------------------------- #
def test_haplotagged_reads_match_the_haplotype_they_came_from(make_locus,
                                                              run_longhap,
                                                              tmp_path):
    import pysam
    locus = make_locus(lambda seq: alternating_snvs(seq, POSITIONS, PATTERN))
    lh = run_longhap(locus, output_bam=str(tmp_path / 'tagged.bam'),
                     output_read_assignments=str(tmp_path / 'assign.tsv'))

    # reads were named r<i> with haplotype i % 2 by make_locus's default tiling
    with pysam.AlignmentFile(lh.output_bam, check_sq=False) as bam:
        tagged = {r.query_name: r.get_tag('HP') for r in bam if r.has_tag('HP')}
    assert tagged, 'nothing was tagged'

    groups = {0: set(), 1: set()}
    for name, hp in tagged.items():
        groups[int(name[4:]) % 2].add(hp)
    # each true haplotype must land consistently on one tag, and the two differ
    assert all(len(v) == 1 for v in groups.values()), groups
    assert groups[0] != groups[1], groups
