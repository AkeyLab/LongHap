"""Degenerate inputs and boundary conditions.

Each test states the behaviour longhap *should* have.  Tests marked ``xfail``
pin a defect written up in ``longhap_claude_issues.md``; when the fix lands they
turn into unexpected passes and the marker comes off.
"""
import numpy as np
import pytest

from synthetic import ReadSpec, alternating_snvs, ref_snv


# --------------------------------------------------------------------------- #
# empty and near-empty loci
# --------------------------------------------------------------------------- #
def test_locus_with_no_heterozygous_variants(make_locus, run_longhap, read_phased_vcf):
    positions = [400, 800, 1200]
    locus = make_locus(lambda seq: alternating_snvs(seq, positions),
                       genotypes=['1/1'] * 3)
    lh = run_longhap(locus)
    assert lh.num_variants == 0
    # nothing is phaseable, but the outputs must still be produced so a
    # per-chromosome pipeline does not lose a file
    assert not lh.block_ends
    records = read_phased_vcf(lh.output_vcf)
    assert [pos for pos, _, _ in records] == positions
    assert all('|' not in gt for _, gt, _ in records)


def test_locus_with_a_single_heterozygous_variant(make_locus, run_longhap):
    locus = make_locus(lambda seq: [ref_snv(seq, 800)])
    lh = run_longhap(locus, write=False)
    assert lh.num_variants == 1
    # a lone variant has no junction, so it cannot be phased
    assert lh.haplotypes.shape[1] == 1


def test_locus_with_two_heterozygous_variants(make_locus, run_longhap):
    """Two variants on one read is the smallest phaseable locus."""
    locus = make_locus(lambda seq: alternating_snvs(seq, [800, 900], '01'))
    lh = run_longhap(locus, write=False)
    assert lh.num_variants == 2
    assert (lh.haplotypes[0] != -1).all()


def test_locus_with_no_reads(make_locus, run_longhap):
    positions = [400, 800, 1200]
    locus = make_locus(lambda seq: alternating_snvs(seq, positions), read_specs=[])
    lh = run_longhap(locus, write=False)
    assert lh.allele_coverage.sum() == 0
    assert (lh.haplotypes[0] == -1).all(), \
        f'nothing is phaseable without reads, got {lh.haplotypes[0]}'


def test_trailing_variant_with_no_junction_is_left_unphased(make_locus, run_longhap):
    """The last site sits far past read reach, so nothing connects it."""
    positions = [400, 500, 600, 3800]
    reads = [ReadSpec(300, 700, i % 2) for i in range(20)]
    locus = make_locus(lambda seq: alternating_snvs(seq, positions, '0101'),
                       read_specs=reads, length=4000)
    lh = run_longhap(locus, write=False)
    assert lh.haplotypes[0][3] == -1, \
        f'unsupported trailing variant reported as phased: {lh.haplotypes[0]}'


# --------------------------------------------------------------------------- #
# output plumbing
# --------------------------------------------------------------------------- #
def test_blocks_file_when_nothing_is_phaseable(make_locus, run_longhap, tmp_path):
    """Every transition stays uninformative, so the scan runs off the end."""
    positions = [400, 800, 1200, 1600]
    # reads never span two variants, so no junction is ever observed
    reads = [ReadSpec(p - 60, p + 60, i % 2) for i, p in enumerate(positions * 6)]
    locus = make_locus(lambda seq: alternating_snvs(seq, positions), read_specs=reads)
    lh = run_longhap(locus)
    assert open(lh.output_blocks).read() == ''


def test_haplotagging_without_a_read_assignment_file(make_locus, run_longhap, tmp_path):
    positions = list(range(400, 3600, 200))
    locus = make_locus(lambda seq: alternating_snvs(seq, positions))
    lh = run_longhap(locus, output_bam=str(tmp_path / 'tagged.bam'),
                     output_read_assignments=None)
    assert lh.output_bam


def test_haplotagging_with_a_read_assignment_file(make_locus, run_longhap, tmp_path):
    import pysam
    positions = list(range(400, 3600, 200))
    locus = make_locus(lambda seq: alternating_snvs(seq, positions))
    lh = run_longhap(locus, output_bam=str(tmp_path / 'tagged.bam'),
                     output_read_assignments=str(tmp_path / 'assignments.tsv'))
    with pysam.AlignmentFile(lh.output_bam, check_sq=False) as bam:
        tags = [r.get_tag('HP') for r in bam if r.has_tag('HP')]
    assert tags, 'no read was haplotagged'
    assert set(tags) == {1, 2}, f'reads landed on only one haplotype: {set(tags)}'


def test_phase_set_ids_are_stable_within_a_block(make_locus, run_longhap, read_phased_vcf):
    positions = list(range(400, 3600, 200))
    locus = make_locus(lambda seq: alternating_snvs(seq, positions))
    lh = run_longhap(locus)
    ps = {r[2] for r in read_phased_vcf(lh.output_vcf) if '|' in r[1]}
    assert len(ps) == 1, f'one block should mean one PS, got {ps}'


# --------------------------------------------------------------------------- #
# read-level edge cases
# --------------------------------------------------------------------------- #
def test_supplementary_alignments_are_ignored(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    good = [ReadSpec(s, s + 400, i % 2, name=f'r{i}')
            for i, s in enumerate(range(0, 3600, 100))]
    # same names, opposite haplotype, flagged supplementary
    supp = [ReadSpec(s, s + 400, (i + 1) % 2, name=f'r{i}', is_supplementary=True)
            for i, s in enumerate(range(0, 3600, 100))]

    baseline = make_locus(lambda seq: alternating_snvs(seq, positions), read_specs=good)
    expected = run_longhap(baseline, phase=False).allele_coverage

    locus = make_locus(lambda seq: alternating_snvs(seq, positions),
                       read_specs=good + supp)
    got = run_longhap(locus, phase=False).allele_coverage
    assert np.array_equal(got, expected), \
        f'supplementary alignments changed coverage:\n{got}\nvs\n{expected}'


def test_reverse_strand_reads_are_used(make_locus, run_longhap):
    positions = list(range(400, 3600, 200))
    reads = [ReadSpec(s, s + 400, i % 2, is_reverse=(i % 3 == 0))
             for i, s in enumerate(range(0, 3600, 100))]
    locus = make_locus(lambda seq: alternating_snvs(seq, positions), read_specs=reads)
    lh = run_longhap(locus)
    assert '.' not in ''.join('.' if h == -1 else 'x' for h in lh.haplotypes[0])


def test_variant_at_the_very_first_base_of_the_contig(make_locus, run_longhap):
    """Realignment windows reach upstream; position 1 has nothing upstream."""
    locus = make_locus(lambda seq: alternating_snvs(seq, [1, 200, 400], '010'),
                       read_specs=[ReadSpec(0, 500, i % 2) for i in range(20)])
    lh = run_longhap(locus, write=False)
    assert lh.num_variants == 3


# --------------------------------------------------------------------------- #
# variant filtering
# --------------------------------------------------------------------------- #
def test_max_allele_length_excludes_long_alleles(make_locus, run_longhap):
    from synthetic import ref_insertion

    def build(seq):
        variants = alternating_snvs(seq, [400, 800, 1200])
        variants[1] = ref_insertion(seq, 800, 'A' * 60)
        return variants

    locus = make_locus(build)
    lh = run_longhap(locus, phase=False, max_allele_length=20)
    assert lh.num_variants == 2, 'the 61 bp allele should have been filtered out'


def test_snvs_only_excludes_indels(make_locus, run_longhap):
    from synthetic import ref_insertion

    def build(seq):
        variants = alternating_snvs(seq, [400, 800, 1200])
        variants[1] = ref_insertion(seq, 800, 'TTAG')
        return variants

    locus = make_locus(build)
    lh = run_longhap(locus, phase=False, snvs_only=True)
    assert lh.num_variants == 2
    assert list(lh.variant_type) == ['SNP', 'SNP']


def test_multiallelic_sites_are_excluded_by_default(make_locus, run_longhap):
    from synthetic import Variant

    def build(seq):
        variants = alternating_snvs(seq, [400, 800, 1200])
        base = seq[799]
        others = [b for b in 'ACGT' if b != base][:2]
        variants[1] = Variant(800, base, others, hap=(1, 2))
        return variants

    locus = make_locus(build, genotypes=['0/1', '1/2', '0/1'])
    lh = run_longhap(locus, phase=False)
    assert lh.num_variants == 2

    lh_multi = run_longhap(locus, phase=False, multiallelics=True)
    assert lh_multi.num_variants == 3
    assert 'MULTI' in lh_multi.variant_type


def test_lv_vcf_split_variants_stay_in_register(make_locus, run_longhap, read_phased_vcf):
    """A vg/pggb VCF marks split complex variants with '_' in the ID."""
    lv_header = '##INFO=<ID=LV,Number=1,Type=Integer,Description="Level">\n'
    positions = list(range(400, 3600, 200))

    def build(seq):
        variants = alternating_snvs(seq, positions, '0101010101010101'[:len(positions)])
        for n, v in enumerate(variants):
            v.vid = f'v{n}_1' if n == 5 else f'v{n}'
        return variants

    # dropping the split variant widens that junction to 400 bp, so the reads
    # have to be long enough to span it or the block legitimately breaks
    locus = make_locus(build, infos={n: 'LV=0' for n in range(len(positions))},
                       extra_header=lv_header, read_length=600)
    lh = run_longhap(locus)
    records = read_phased_vcf(lh.output_vcf)
    phased = [(pos, gt) for pos, gt, _ in records if '|' in gt]
    truth = {v.pos: str(v.hap[0]) for v in locus.variants}
    # a register slip shows up as the phase of one site being attributed to another
    flip = {'0': '1', '1': '0'}
    same = all(gt[0] == truth[pos] for pos, gt in phased)
    comp = all(gt[0] == flip[truth[pos]] for pos, gt in phased)
    assert same or comp, f'phases do not line up with the truth: {phased}'
