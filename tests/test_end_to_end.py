"""The whole pipeline on small synthetic loci with a known answer.

A locus is built from an explicit truth haplotype, reads are rendered from it,
and longhap is asked to recover it.  Phasing is only defined up to a global
flip, so assertions compare the *pattern* of the recovered haplotype against the
truth and its complement.
"""
import numpy as np
import pytest

from synthetic import ReadSpec, alternating_snvs, ref_deletion, ref_insertion, ref_snv


def phase_string(locus, lh):
    """Recovered haplotype 0 as a '0'/'1' string, '.' where unphased."""
    return ''.join('.' if h == -1 else str(h) for h in lh.haplotypes[0])


def truth_string(locus):
    """Which allele haplotype 0 carries at each site, as a '0'/'1' string."""
    return ''.join(str(v.hap[0]) for v in locus.variants)


def matches_up_to_flip(got, want):
    """True when ``got`` equals ``want`` or its complement, ignoring '.' sites."""
    flipped = ''.join({'0': '1', '1': '0'}.get(c, c) for c in want)
    same = all(g == w for g, w in zip(got, want) if g != '.')
    comp = all(g == f for g, f in zip(got, flipped) if g != '.')
    return same or comp


# --------------------------------------------------------------------------- #
# the happy path
# --------------------------------------------------------------------------- #
def test_phases_a_simple_snv_locus(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    pattern = '0101100110110100'[:len(positions)]
    locus = make_locus(lambda seq: alternating_snvs(seq, positions, pattern))
    lh = run_longhap(locus)

    assert lh.num_variants == len(positions)
    got = phase_string(locus, lh)
    assert '.' not in got, f'left {got.count(".")} variants unphased: {got}'
    assert matches_up_to_flip(got, truth_string(locus)), \
        f'recovered {got}, truth {truth_string(locus)}'


def test_one_phase_block_when_reads_span_every_junction(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    locus = make_locus(lambda seq: alternating_snvs(seq, positions, '0' * len(positions)))
    lh = run_longhap(locus)
    # every junction is covered, so there should be exactly one block
    assert len(lh.block_ends) == 1, f'block_ends={lh.block_ends}'


def test_phased_vcf_round_trips(make_locus, run_longhap, read_phased_vcf):
    positions = list(range(400, 3600, 200))
    pattern = '0110100101101001'[:len(positions)]
    locus = make_locus(lambda seq: alternating_snvs(seq, positions, pattern))
    lh = run_longhap(locus)
    records = read_phased_vcf(lh.output_vcf)

    assert [pos for pos, _, _ in records] == positions
    assert all('|' in gt for _, gt, _ in records), records
    got = ''.join(gt[0] for _, gt, _ in records)
    assert matches_up_to_flip(got, truth_string(locus)), f'{got} vs {truth_string(locus)}'


def test_allele_coverage_counts_both_alleles(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    locus = make_locus(lambda seq: alternating_snvs(seq, positions))
    lh = run_longhap(locus, phase=False)
    assert (lh.allele_coverage > 0).all(), lh.allele_coverage


def test_transition_matrix_rows_are_normalised(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    locus = make_locus(lambda seq: alternating_snvs(seq, positions))
    lh = run_longhap(locus, phase=False)
    sums = lh.transition_matrix.sum(axis=1)
    assert np.allclose(sums, 1.0), sums


def test_haplotypes_are_complementary(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    locus = make_locus(lambda seq: alternating_snvs(seq, positions, '0110110010101001'[:16]))
    lh = run_longhap(locus)
    phased = lh.haplotypes[0] != -1
    assert (lh.haplotypes[0][phased] != lh.haplotypes[1][phased]).all()


# --------------------------------------------------------------------------- #
# indels
# --------------------------------------------------------------------------- #
def test_phases_a_locus_containing_an_insertion(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))

    def build(seq):
        variants = alternating_snvs(seq, positions, '0' * len(positions))
        variants[5] = ref_insertion(seq, positions[5], 'TTAGGC')
        return variants

    locus = make_locus(build)
    lh = run_longhap(locus)
    got = phase_string(locus, lh)
    assert got[5] != '.', f'insertion left unphased: {got}'
    assert matches_up_to_flip(got, truth_string(locus)), got


def test_phases_a_locus_containing_a_deletion(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))

    def build(seq):
        variants = alternating_snvs(seq, positions, '0' * len(positions))
        variants[5] = ref_deletion(seq, positions[5], length=6)
        return variants

    locus = make_locus(build)
    lh = run_longhap(locus)
    got = phase_string(locus, lh)
    assert got[5] != '.', f'deletion left unphased: {got}'
    assert matches_up_to_flip(got, truth_string(locus)), got


# --------------------------------------------------------------------------- #
# block structure
# --------------------------------------------------------------------------- #
def test_uncovered_junction_breaks_the_block(make_locus, run_longhap):
    """No read spans the gap between the two clusters, so phasing must break."""
    positions = [300, 400, 500, 2500, 2600, 2700]
    reads = ([ReadSpec(200, 600, i % 2) for i in range(20)] +
             [ReadSpec(2400, 2800, i % 2) for i in range(20)])
    locus = make_locus(lambda seq: alternating_snvs(seq, positions, '010101'),
                       read_specs=reads)
    lh = run_longhap(locus)
    assert len(lh.block_ends) >= 2, f'expected a break, got block_ends={lh.block_ends}'


def test_variant_records_are_all_written_back(make_locus, run_longhap, read_phased_vcf):
    """Every input record must appear in the output, phased or not."""
    positions = list(range(400, 3600, 200))
    genotypes = ['0/1'] * len(positions)
    genotypes[3] = '1/1'          # homozygous, longhap must pass it through
    genotypes[7] = './.'          # missing, likewise
    locus = make_locus(lambda seq: alternating_snvs(seq, positions), genotypes=genotypes)
    lh = run_longhap(locus)
    records = read_phased_vcf(lh.output_vcf)
    assert [pos for pos, _, _ in records] == positions


def test_haplotype_blocks_bed_is_written(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    locus = make_locus(lambda seq: alternating_snvs(seq, positions))
    lh = run_longhap(locus)
    rows = [line.split() for line in open(lh.output_blocks) if line.strip()]
    assert rows, 'no phase blocks written'
    for chrom, start, end in rows:
        assert chrom == locus.chrom
        assert int(start) <= int(end)


# --------------------------------------------------------------------------- #
# read filtering
# --------------------------------------------------------------------------- #
def test_low_mapq_reads_are_ignored(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    reads = [ReadSpec(s, s + 400, i % 2, mapq=0)
             for i, s in enumerate(range(0, 3600, 100))]
    locus = make_locus(lambda seq: alternating_snvs(seq, positions), read_specs=reads)
    lh = run_longhap(locus, phase=False)
    assert lh.allele_coverage.sum() == 0, 'MAPQ 0 reads were counted'


def test_secondary_alignments_are_ignored(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    good = [ReadSpec(s, s + 400, i % 2, name=f'r{i}')
            for i, s in enumerate(range(0, 3600, 100))]
    noise = [ReadSpec(s, s + 400, (i + 1) % 2, name=f'sec{i}', is_secondary=True)
             for i, s in enumerate(range(0, 3600, 100))]
    locus = make_locus(lambda seq: alternating_snvs(seq, positions),
                       read_specs=good + noise)
    lh = run_longhap(locus, phase=False)
    baseline = make_locus(lambda seq: alternating_snvs(seq, positions), read_specs=good)
    lh_ref = run_longhap(baseline, phase=False)
    assert np.array_equal(lh.allele_coverage, lh_ref.allele_coverage)
