"""Properties that must hold for any input, checked over randomised datasets.

Two kinds of check live here.  First, naive reference implementations of the two
pieces of real algorithmic machinery -- the exon sweep and the anchor scan --
compared against the optimised versions on random inputs; these are oracles, not
restatements of the code.  Second, metamorphic properties: transformations of the
input that must leave the statistics unchanged, or change them in a known way.
A hand-written fixture can only assert one point; these assert a whole class.
"""
import random

import pytest

from evaluate_phasing import (evaluate_genes, evaluate_junctions, evaluate_phasing,
                              evaluate_svs, genes_by_variant, get_overlapping_sites,
                              load_annotations, load_phasing, nearest_snp_indices)

SEEDS = list(range(12))


# --------------------------------------------------------------------------
# reference implementations
# --------------------------------------------------------------------------

def naive_genes_by_variant(indices, target, chrom_exons):
    """O(variants x exons): test every variant against every exon."""
    out = {}
    for i in indices:
        start = int(target.position[i])
        end = start + len(target.ref[i])
        for exon_start, exon_end, gene in chrom_exons:
            if start < exon_end and end > exon_start:
                out.setdefault(gene, [])
                if i not in out[gene]:
                    out[gene].append(i)
    return out


def naive_nearest_snp_indices(chain, target):
    """Scan outward from every position looking for a SNP."""
    left, right = [], []
    for k in range(len(chain)):
        prev = [j for j in range(k) if target.variant_type[chain[j][0]] == 'SNP']
        nxt = [j for j in range(k + 1, len(chain))
               if target.variant_type[chain[j][0]] == 'SNP']
        left.append(prev[-1] if prev else None)
        right.append(nxt[0] if nxt else None)
    return left, right


# --------------------------------------------------------------------------
# random dataset generation
# --------------------------------------------------------------------------

GENOTYPES = ['0|1', '1|0']


def _records(rng, chroms=('chr1',), n=40, unphased_rate=0.15, indel_rate=0.3):
    """A sorted, position-unique set of records with random types and phase sets."""
    records = []
    for chrom in chroms:
        pos = sorted(rng.sample(range(10, 10000), n))
        block = 1
        for p in pos:
            if rng.random() < 0.12:
                block += 1
            if rng.random() < indel_rate:
                ref, alt = ('A', 'A' + 'T' * rng.randint(1, 4)) if rng.random() < 0.5 \
                    else ('A' + 'T' * rng.randint(1, 4), 'A')
            else:
                ref, alt = rng.sample('ACGT', 2)
            gt = rng.choice(GENOTYPES)
            if rng.random() < unphased_rate:
                gt = gt.replace('|', '/')
            records.append((chrom, p, ref, alt, gt, block))
    return records


def _truth_of(records, rng, disagree_rate=0.1):
    """Truth carrying the same variants, chromosome-wide, with planted flips."""
    out = []
    for chrom, p, ref, alt, gt, _ in records:
        base = gt.replace('/', '|')
        if rng.random() < disagree_rate:
            base = '1|0' if base == '0|1' else '0|1'
        out.append((chrom, p, ref, alt, base, None))
    return out


def _invert(records):
    """Swap the two haplotypes of every record."""
    flipped = []
    for chrom, p, ref, alt, gt, ps in records:
        sep = '|' if '|' in gt else '/'
        a, b = gt.split(sep)
        flipped.append((chrom, p, ref, alt, f'{b}{sep}{a}', ps))
    return flipped


def _shift(records, offset):
    return [(c, p + offset, r, a, gt, ps) for c, p, r, a, gt, ps in records]


def _all_stats(write_vcf, target_records, truth_records, name='a'):
    target = load_phasing(write_vcf(f'{name}t.vcf', target_records,
                                    contigs=('chr1', 'chr2')))
    truth = load_phasing(write_vcf(f'{name}g.vcf', truth_records, contigs=('chr1', 'chr2')))
    pairs = get_overlapping_sites(target, truth)
    within = evaluate_phasing(pairs, target, truth)
    junctions = evaluate_junctions(pairs, target, truth)
    chrwide = evaluate_phasing(pairs, target, truth, ignore_phase_blocks=True)
    svs = evaluate_svs(pairs, target, truth)
    return dict(
        assessed=within.assessed_pairs, switches=within.switches, hamming=within.hamming,
        flips=within.switch_flips, blocks=within.intersection_blocks,
        junctions=junctions.junctions, junction_errors=junctions.errors,
        chr_assessed=chrwide.assessed_pairs, chr_switches=chrwide.switches,
        sv_total=svs.total, sv_evaluated=svs.svs, sv_errors=svs.errors,
        sv_correct=svs.correct, sv_flipped=svs.flipped, sv_ambiguous=svs.ambiguous,
    )


# --------------------------------------------------------------------------
# oracle comparisons
# --------------------------------------------------------------------------

@pytest.mark.parametrize('seed', SEEDS)
def test_exon_sweep_matches_a_naive_scan(seed, write_vcf):
    rng = random.Random(seed)
    records = _records(rng, n=30)
    target = load_phasing(write_vcf('t.vcf', records))
    exons = sorted((s, s + rng.randint(1, 400), f'G{rng.randint(0, 6)}')
                   for s in (rng.randint(0, 10000) for _ in range(40)))
    indices = list(range(len(target)))

    fast = {g: sorted(v) for g, v in genes_by_variant(indices, target, exons).items()}
    slow = {g: sorted(v) for g, v in naive_genes_by_variant(indices, target, exons).items()}
    assert fast == slow


@pytest.mark.parametrize('seed', SEEDS)
def test_anchor_scan_matches_a_naive_scan(seed, write_vcf):
    rng = random.Random(seed)
    target = load_phasing(write_vcf('t.vcf', _records(rng, n=30)))
    chain = [(i, i) for i in range(len(target))]
    assert nearest_snp_indices(chain, target) == naive_nearest_snp_indices(chain, target)


# --------------------------------------------------------------------------
# metamorphic properties
# --------------------------------------------------------------------------

@pytest.mark.parametrize('seed', SEEDS)
def test_within_and_junctions_partition_the_chromosome_wide_pairs(seed, write_vcf):
    rng = random.Random(seed)
    records = _records(rng)
    stats = _all_stats(write_vcf, records, _truth_of(records, rng))
    assert stats['assessed'] + stats['junctions'] == stats['chr_assessed']
    assert stats['switches'] + stats['junction_errors'] == stats['chr_switches']


@pytest.mark.parametrize('seed', SEEDS)
def test_inverting_the_target_changes_nothing(seed, write_vcf):
    """Absolute orientation is arbitrary; only relative orientation is measured."""
    rng = random.Random(seed)
    records = _records(rng)
    truth = _truth_of(records, rng)
    assert _all_stats(write_vcf, records, truth, 'a') == \
        _all_stats(write_vcf, _invert(records), truth, 'b')


@pytest.mark.parametrize('seed', SEEDS)
def test_inverting_the_truth_changes_nothing(seed, write_vcf):
    rng = random.Random(seed)
    records = _records(rng)
    truth = _truth_of(records, rng)
    assert _all_stats(write_vcf, records, truth, 'a') == \
        _all_stats(write_vcf, records, _invert(truth), 'b')


@pytest.mark.parametrize('seed', SEEDS)
def test_shifting_all_coordinates_changes_nothing(seed, write_vcf):
    rng = random.Random(seed)
    records = _records(rng)
    truth = _truth_of(records, rng)
    assert _all_stats(write_vcf, records, truth, 'a') == \
        _all_stats(write_vcf, _shift(records, 1_000_000), _shift(truth, 1_000_000), 'b')


@pytest.mark.parametrize('seed', SEEDS[:6])
def test_duplicating_the_data_onto_a_second_chromosome_doubles_everything(seed, write_vcf):
    """Nothing may leak across chromosomes, so every count doubles exactly."""
    rng = random.Random(seed)
    records = _records(rng)
    truth = _truth_of(records, rng)
    one = _all_stats(write_vcf, records, truth, 'a')

    def onto_chr2(rs):
        return [('chr2', p, r, a, gt, ps) for _, p, r, a, gt, ps in rs]
    two = _all_stats(write_vcf, records + onto_chr2(records),
                     truth + onto_chr2(truth), 'b')

    for key, value in one.items():
        if key == 'flips':
            assert two[key] == (value[0] * 2, value[1] * 2), key
        else:
            assert two[key] == value * 2, key


@pytest.mark.parametrize('seed', SEEDS)
def test_gene_pairs_are_a_subset_of_the_within_block_pairs(seed, write_vcf, write_bed):
    """Cross-check between two independently written code paths.

    The two deliberately pair different things: the gene metric walks *every*
    heterozygous site, so a pair touching an unphased site is unresolved, while
    the within-block metric walks only the scorable chain and pairs across such a
    site.  A gene pair that does get scored is therefore always also a
    within-block pair -- both its ends are scorable and nothing sits between them
    -- so the gene statistics can never exceed the within-block ones.
    """
    rng = random.Random(seed)
    records = _records(rng)
    target = load_phasing(write_vcf('t.vcf', records))
    truth = load_phasing(write_vcf('g.vcf', _truth_of(records, rng)))
    pairs = get_overlapping_sites(target, truth)
    exons = load_annotations(write_bed('e.bed', [('chr1', 0, 100_000_000, 'ALL')]))

    genes = evaluate_genes(pairs, target, truth, exons)
    within = evaluate_phasing(pairs, target, truth)
    assert genes.connections <= within.assessed_pairs
    assert genes.errors <= within.switches
    assert set(genes.positions) <= set(within.switch_positions)
    assert genes.connections + genes.unresolved == genes.sites - genes.genes


@pytest.mark.parametrize('seed', SEEDS)
def test_gene_and_switch_metrics_coincide_when_every_site_is_phased(seed, write_vcf,
                                                                   write_bed):
    """With nothing unphased the two pairings are identical, so the counts must be."""
    rng = random.Random(seed)
    records = _records(rng, unphased_rate=0.0)
    target = load_phasing(write_vcf('t.vcf', records))
    truth = load_phasing(write_vcf('g.vcf', _truth_of(records, rng)))
    pairs = get_overlapping_sites(target, truth)
    exons = load_annotations(write_bed('e.bed', [('chr1', 0, 100_000_000, 'ALL')]))

    genes = evaluate_genes(pairs, target, truth, exons)
    within = evaluate_phasing(pairs, target, truth)
    assert genes.connections == within.assessed_pairs
    assert genes.errors == within.switches
    assert sorted(genes.positions) == sorted(within.switch_positions)


@pytest.mark.parametrize('seed', SEEDS)
def test_sv_identities_hold_on_random_data(seed, write_vcf):
    rng = random.Random(seed)
    records = _records(rng)
    target = load_phasing(write_vcf('t.vcf', records))
    truth = load_phasing(write_vcf('g.vcf', _truth_of(records, rng)))
    res = evaluate_svs(get_overlapping_sites(target, truth), target, truth)

    assert res.correct + res.flipped + res.ambiguous + res.one_sided + res.no_anchor == res.svs
    assert res.connections == 2 * (res.correct + res.flipped + res.ambiguous) + res.one_sided
    assert res.svs <= res.total
    assert 2 * res.flipped + res.ambiguous <= res.errors <= res.connections
    assert sum(res.by_type.values()) == res.svs


@pytest.mark.parametrize('seed', SEEDS)
def test_min_sv_length_is_monotonic(seed, write_vcf):
    rng = random.Random(seed)
    records = _records(rng, indel_rate=0.5)
    target = load_phasing(write_vcf('t.vcf', records))
    truth = load_phasing(write_vcf('g.vcf', _truth_of(records, rng)))
    pairs = get_overlapping_sites(target, truth)
    counts = [evaluate_svs(pairs, target, truth, min_sv_length=n).total
              for n in (0, 2, 3, 5, 50)]
    assert counts == sorted(counts, reverse=True)


@pytest.mark.parametrize('seed', SEEDS[:6])
def test_evaluation_is_deterministic(seed, write_vcf):
    rng = random.Random(seed)
    records = _records(rng)
    truth = _truth_of(records, rng)
    assert _all_stats(write_vcf, records, truth, 'a') == \
        _all_stats(write_vcf, records, truth, 'b')
