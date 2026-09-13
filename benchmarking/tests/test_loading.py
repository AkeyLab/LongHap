"""load_phasing and load_annotations."""
import numpy as np
import pytest

from evaluate_phasing import classify, load_annotations, load_phasing


def test_classify():
    assert classify('A', ('G',)) == 'SNP'
    assert classify('AT', ('A',)) == 'DEL'
    assert classify('A', ('AT',)) == 'INS'
    assert classify('AT', ('GC',)) == 'Other'      # equal length, not 1 bp
    assert classify('A', ('G', 'T')) == 'Multi'


def test_heterozygous_only(write_vcf):
    path = write_vcf('a.vcf', [
        ('chr1', 10, 'A', 'G', '0|1', 1),
        ('chr1', 20, 'C', 'T', '1|1', 1),      # hom alt, excluded
        ('chr1', 30, 'G', 'A', '0/0', None),   # hom ref, excluded
        ('chr1', 40, 'T', 'C', '1|0', 1),
    ])
    vs = load_phasing(path)
    assert list(vs.position) == [9, 39]         # 0-based
    assert vs.n_total == 4                      # every record, het or not


def test_missing_allele_is_heterozygous_but_incomplete(write_vcf):
    """whatshap treats a genotype with a missing allele as not homozygous."""
    path = write_vcf('a.vcf', [
        ('chr1', 10, 'A', 'G', '.|1', 1),
        ('chr1', 20, 'C', 'T', '0|.', 1),
        ('chr1', 30, 'G', 'A', '.|.', 1),
        ('chr1', 40, 'T', 'C', '0|1', 1),
    ])
    vs = load_phasing(path)
    assert len(vs) == 4
    assert list(vs.complete) == [False, False, False, True]
    assert vs.allele_a[0] is None and vs.allele_a[3] == 'T'


def test_duplicate_positions_dropped_like_whatshap(write_vcf):
    records = [
        ('chr1', 10, 'A', 'G', '0|1', 1),
        ('chr1', 10, 'A', 'T', '0|1', 1),      # same position, dropped
        ('chr1', 20, 'C', 'T', '0|1', 1),
    ]
    path = write_vcf('a.vcf', records)
    vs = load_phasing(path)
    assert len(vs) == 2 and vs.n_duplicate == 1 and vs.n_total == 2

    kept = load_phasing(path, keep_duplicate_positions=True)
    assert len(kept) == 3 and kept.n_duplicate == 0 and kept.n_total == 3


def test_duplicate_position_check_is_per_chromosome(write_vcf):
    """The same coordinate on two chromosomes is not a duplicate."""
    path = write_vcf('a.vcf', [
        ('chr1', 10, 'A', 'G', '0|1', 1),
        ('chr2', 10, 'A', 'G', '0|1', 1),
    ], contigs=('chr1', 'chr2'))
    assert len(load_phasing(path)) == 2


def test_phase_set_and_phased_flag(write_vcf):
    path = write_vcf('a.vcf', [
        ('chr1', 10, 'A', 'G', '0|1', 7),
        ('chr1', 20, 'C', 'T', '0/1', 7),      # unphased genotype
        ('chr1', 30, 'G', 'A', '0|1', None),   # phased, no PS
    ])
    vs = load_phasing(path)
    assert list(vs.phased) == [True, False, True]
    assert list(vs.phase_block) == [7, 7, -1]  # absent PS gets the sentinel


def test_only_snvs(write_vcf):
    path = write_vcf('a.vcf', [
        ('chr1', 10, 'A', 'G', '0|1', 1),
        ('chr1', 20, 'CA', 'C', '0|1', 1),
        ('chr1', 30, 'G', 'GT', '0|1', 1),
    ])
    assert len(load_phasing(path)) == 3
    assert list(load_phasing(path, only_snvs=True).variant_type) == ['SNP']


def test_multiallelic_cap(write_vcf):
    """whatshap refuses sites with 16 or more alleles."""
    many = ','.join('ACGT'[i % 4] * (i + 2) for i in range(16))
    path = write_vcf('a.vcf', [
        ('chr1', 10, 'A', many, '0|1', 1),
        ('chr1', 20, 'C', 'T,G', '1|2', 1),
    ])
    vs = load_phasing(path)
    assert list(vs.position) == [19]
    assert list(vs.variant_type) == ['Multi']
    assert (vs.allele_a[0], vs.allele_b[0]) == ('T', 'G')   # 1|2 picks the two ALTs


def test_sample_selection(write_vcf):
    path = write_vcf('a.vcf', [
        ('chr1', 10, 'A', 'G', ('0|1', '1|1'), (1, 1)),
        ('chr1', 20, 'C', 'T', ('1|1', '0|1'), (1, 1)),
    ], samples=('S1', 'S2'))
    assert list(load_phasing(path, sample='S1').position) == [9]
    assert list(load_phasing(path, sample='S2').position) == [19]


def test_missing_sample_is_an_error(write_vcf):
    path = write_vcf('a.vcf', [('chr1', 10, 'A', 'G', '0|1', 1)])
    with pytest.raises(SystemExit):
        load_phasing(path, sample='nope')


def test_allele_arrays_are_object_dtype(write_vcf):
    """Regression: a fixed-width unicode array over a long allele needs tens of GB."""
    long_alt = 'A' * 50000
    path = write_vcf('a.vcf', [
        ('chr1', 10, 'A', long_alt, '0|1', 1),
        ('chr1', 20, 'C', 'T', '0|1', 1),
    ])
    vs = load_phasing(path)
    for arr in (vs.ref, vs.alt, vs.allele_a, vs.allele_b):
        assert arr.dtype == object
    assert vs.ref.nbytes < 1000


def test_alt_array_stays_one_dimensional(write_vcf):
    """Regression: np.array over a list of 1-tuples builds a 2-D array."""
    path = write_vcf('a.vcf', [
        ('chr1', 10, 'A', 'G', '0|1', 1),
        ('chr1', 20, 'C', 'T', '0|1', 1),
    ])
    vs = load_phasing(path)
    assert vs.alt.ndim == 1
    assert vs.alt[0] == ('G',)


def test_empty_and_all_homozygous(write_vcf):
    assert len(load_phasing(write_vcf('e.vcf', []))) == 0
    hom = write_vcf('h.vcf', [('chr1', 10, 'A', 'G', '1|1', 1)])
    vs = load_phasing(hom)
    assert len(vs) == 0 and vs.n_total == 1


def test_load_annotations_filters_and_sorts(write_bed):
    path = write_bed('e.bed', [
        ('chr2', 50, 60, 'G2'),
        ('chr1', 30, 40, 'G1b'),
        ('chr1', 10, 20, 'G1a'),
        ('#comment', 0, 0, 'x'),
    ])
    exons = load_annotations(path)
    assert set(exons) == {'chr1', 'chr2'}
    assert exons['chr1'] == [(10, 20, 'G1a'), (30, 40, 'G1b')]

    only1 = load_annotations(path, {'chr1'})
    assert set(only1) == {'chr1'}


def test_load_annotations_gzip(write_bed):
    path = write_bed('e.bed.gz', [('chr1', 10, 20, 'G1')])
    assert load_annotations(path)['chr1'] == [(10, 20, 'G1')]


def test_load_annotations_rejects_short_rows(write_bed):
    path = write_bed('bad.bed', [('chr1', 10, 20)])
    with pytest.raises(SystemExit):
        load_annotations(path)
