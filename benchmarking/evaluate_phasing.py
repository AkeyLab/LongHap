#!/usr/bin/env python3
"""Evaluate variant phasing against a ground truth phasing.

The variant counting and switch-error statistics are defined to reproduce
``whatshap compare`` exactly (see the notes on each option below), so results
from this script can be checked against it directly.
"""
import argparse
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from cyvcf2 import VCF
import numpy as np

# whatshap refuses to model sites with this many alleles or more (whatshap.core)
MAX_GENOTYPE_ALLELES = 16

COUNT_WIDTH = 9
LABEL_WIDTH = 40


@dataclass
class VariantSet:
    """Heterozygous variants of a single sample, in VCF order.

    ``position`` is 0-based . ``allele_a``/``allele_b``
    are the actual allele *sequences* carried by haplotype 0 and haplotype 1,
    and are ``None`` where the genotype has a missing allele.
    """
    chromosome: np.ndarray
    position: np.ndarray
    ref: np.ndarray
    alt: np.ndarray
    genotype: np.ndarray      # (n, 2) allele indices, -1 for a missing allele
    phased: np.ndarray
    phase_block: np.ndarray
    variant_type: np.ndarray
    allele_a: np.ndarray
    allele_b: np.ndarray
    complete: np.ndarray      # genotype has no missing allele
    n_records: int = 0        # records read from the file
    n_duplicate: int = 0      # records dropped because a record shared their position
    n_total: int = 0          # records kept (het and non-het), i.e. whatshap's "all"
    all_keys: Optional[set] = None    # identity of every kept record, het or not

    def __len__(self):
        return len(self.position)

    def keys(self):
        """whatshap's variant identity: (chromosome, 0-based position, REF, ALT tuple)."""
        return list(zip(self.chromosome.tolist(), self.position.tolist(),
                        self.ref.tolist(), self.alt.tolist()))


def classify(ref, alt):
    if len(alt) > 1:
        return 'Multi'
    if len(ref) == len(alt[0]) == 1:
        return 'SNP'
    if len(ref) > len(alt[0]):
        return 'DEL'
    if len(ref) < len(alt[0]):
        return 'INS'
    return 'Other'


def load_phasing(vcf_f, sample=None, only_snvs=False, keep_duplicate_positions=False):
    """Load heterozygous variants of one sample from a VCF.

    A variant is heterozygous whenever its genotype is not homozygous; as in
    whatshap, a genotype with a missing allele (``.|1``, ``0|.``, ``.|.``) is
    *not* homozygous and therefore counts as heterozygous.  Such variants can
    never be compared, so they are excluded from the comparison later on, but
    they are part of the reported variant counts.

    Unless ``keep_duplicate_positions`` is set, a record whose position equals
    the previous record's is dropped, which is what whatshap does.  Splitting
    multi-allelic records (``bcftools norm -m -any``) creates many such
    records, so keeping them makes the counts diverge from whatshap.
    """
    variant_calls = VCF(vcf_f)
    if sample is not None:
        variant_calls.set_samples([sample])
    if not variant_calls.samples:
        raise SystemExit(f'no sample to evaluate in {vcf_f}')

    chromosome, position, ref_l, alt_l = [], [], [], []
    genotypes, phased, phase_block = [], [], []
    variant_type, allele_a, allele_b, complete = [], [], [], []
    n_records = n_duplicate = n_total = 0
    all_keys = set()
    prev_key = None

    for variant in variant_calls:
        n_records += 1
        if not variant.ALT:
            continue
        ref = variant.REF.upper()
        alt = tuple(a.upper() for a in variant.ALT)
        if len(alt) >= MAX_GENOTYPE_ALLELES:
            continue
        if only_snvs and not (len(ref) == 1 and all(len(a) == 1 for a in alt)):
            continue

        # whatshap keys on the 0-based start position
        pos = variant.start
        key = (variant.CHROM, pos)
        if not keep_duplicate_positions and key == prev_key:
            n_duplicate += 1
            continue
        prev_key = key
        n_total += 1
        all_keys.add((variant.CHROM, pos, ref, alt))

        gt = list(variant.genotypes[0][:2])
        if len(gt) < 2:
            continue
        if gt[0] == gt[1] and gt[0] != -1:      # homozygous
            continue

        alleles = (ref,) + alt
        has_all = -1 not in gt
        chromosome.append(variant.CHROM)
        position.append(pos)
        ref_l.append(ref)
        alt_l.append(alt)
        genotypes.append(gt)
        phased.append(bool(variant.genotypes[0][2]))
        phase_block.append(int(variant.format('PS')[0][0]) if 'PS' in variant.FORMAT else -1)
        variant_type.append(classify(ref, alt))
        allele_a.append(alleles[gt[0]] if has_all else None)
        allele_b.append(alleles[gt[1]] if has_all else None)
        complete.append(has_all)

    # np.array() on a list of tuples builds a 2-D array, so the ALT tuples have
    # to be placed into an object array one by one to stay 1-D.
    alt_arr = np.empty(len(alt_l), dtype=object)
    for i, a in enumerate(alt_l):
        alt_arr[i] = a

    # NOTE: object dtype is required.  np.array() on a list of allele strings
    # builds a fixed-width unicode array sized by the *longest* allele, which
    # for a VCF containing structural variants (alleles of ~100 kb) needs tens
    # of gigabytes.
    return VariantSet(
        chromosome=np.array(chromosome, dtype=object),
        position=np.array(position, dtype=np.int64),
        ref=np.array(ref_l, dtype=object),
        alt=alt_arr,
        genotype=np.array(genotypes, dtype=np.int16).reshape(-1, 2),
        phased=np.array(phased, dtype=bool),
        phase_block=np.array(phase_block, dtype=np.int64),
        variant_type=np.array(variant_type, dtype=object),
        allele_a=np.array(allele_a, dtype=object),
        allele_b=np.array(allele_b, dtype=object),
        complete=np.array(complete, dtype=bool),
        n_records=n_records,
        n_duplicate=n_duplicate,
        n_total=n_total,
        all_keys=all_keys,
    )


def get_overlapping_sites(target, gt, match='strict'):
    """Index pairs of heterozygous variants present in both call sets.

    ``strict`` matches on (chromosome, 0-based position, REF, ALT), which is
    whatshap's variant identity.  ``alleles`` instead matches on the unordered
    pair of alleles carried by the sample, which additionally pairs up sites
    that the two files represent differently (an insertion ``C -> CA`` in one
    file against a deletion ``CA -> C`` in the other).  Those pairs are not
    the same variant, and since REF differs their allele indices mean
    different things, so ``strict`` is the default.
    """
    index: Dict[tuple, int] = {}
    if match == 'strict':
        for j in range(len(gt)):
            key = (gt.chromosome[j], gt.position[j], gt.ref[j], gt.alt[j])
            index.setdefault(key, j)
        keys = ((target.chromosome[i], target.position[i], target.ref[i], target.alt[i])
                for i in range(len(target)))
    elif match == 'alleles':
        for j in range(len(gt)):
            index.setdefault((gt.chromosome[j], gt.position[j]), j)
        keys = ((target.chromosome[i], target.position[i]) for i in range(len(target)))
    else:
        raise ValueError(f'unknown match mode {match!r}')

    pairs = []
    for i, key in enumerate(keys):
        j = index.get(key)
        if j is None:
            continue
        if match == 'alleles':
            a1, b1, a2, b2 = target.allele_a[i], target.allele_b[i], gt.allele_a[j], gt.allele_b[j]
            if None in (a1, b1, a2, b2):
                continue
            if not ((a1 == a2 and b1 == b2) or (a1 == b2 and b1 == a2)):
                continue
        pairs.append((i, j))
    return np.array(pairs, dtype=np.int64).reshape(-1, 2)


def switch_encoding(phasing: str) -> str:
    """'0001011' -> '001110': 1 wherever consecutive entries differ."""
    return ''.join('0' if phasing[i - 1] == phasing[i] else '1' for i in range(1, len(phasing)))


def hamming(s0, s1) -> int:
    return sum(c0 != c1 for c0, c1 in zip(s0, s1))


def compute_switch_flips(phasing0: str, phasing1: str) -> Tuple[int, int]:
    """Decompose disagreements into switches and flips (a flip = two adjacent switches)."""
    s0, s1 = switch_encoding(phasing0), switch_encoding(phasing1)
    switches = flips = 0
    in_a_row = 0
    for i, (p0, p1) in enumerate(zip(s0, s1)):
        if p0 != p1:
            in_a_row += 1
        if i + 1 == len(s0) or p0 == p1:
            flips += in_a_row // 2
            switches += in_a_row % 2
            in_a_row = 0
    return switches, flips


@dataclass
class ComparisonResult:
    intersection_blocks: int = 0
    covered_variants: int = 0
    assessed_pairs: int = 0
    switches: int = 0
    switch_flips: Tuple[int, int] = (0, 0)
    hamming: int = 0
    diff_genotypes: int = 0
    switch_positions: Optional[List[Tuple[str, int, int]]] = None
    blocks_target: int = 0
    covered_target: int = 0
    blocks_gt: int = 0
    covered_gt: int = 0

    @property
    def switch_rate(self):
        return self.switches / self.assessed_pairs if self.assessed_pairs else float('nan')

    @property
    def switchflip_rate(self):
        total = sum(self.switch_flips)
        return total / self.assessed_pairs if self.assessed_pairs else float('nan')


def group_blocks(overlapping_sites, target, gt, ignore_phase_blocks=False):
    """Assign every comparable variant pair to a phase block in each call set.

    Returns ``(blocks, blocks_target, blocks_gt)``.  ``blocks`` is keyed on the
    *pair* of block ids and holds the intersection blocks -- stretches that are
    contiguously phased in both files.  ``ignore_phase_blocks`` puts every
    phased variant of a chromosome into one block, which is what stripping
    ``FORMAT/PS`` before running whatshap compare does.
    """
    blocks: Dict[tuple, List[tuple]] = defaultdict(list)
    blocks_target: Dict[tuple, List[int]] = defaultdict(list)
    blocks_gt: Dict[tuple, List[int]] = defaultdict(list)

    def block_id(vs, k):
        """Block a variant belongs to, or None if it carries no usable phase."""
        if not (vs.phased[k] and vs.complete[k]):
            return None
        return (vs.chromosome[k], 0 if ignore_phase_blocks else vs.phase_block[k])

    for i, j in overlapping_sites:
        # As in whatshap, the per-file block statistics count a variant as soon
        # as *that* file phases it; only the intersection needs both files.
        b_target, b_gt = block_id(target, i), block_id(gt, j)
        if b_target is not None:
            blocks_target[b_target].append(i)
        if b_gt is not None:
            blocks_gt[b_gt].append(j)
        if b_target is not None and b_gt is not None:
            blocks[(b_target, b_gt)].append((i, j))
    return blocks, blocks_target, blocks_gt


def agreement(block, target, gt):
    """For each variant pair, whether both files put the same allele on haplotype 0.

    Compares allele sequences rather than allele indices, because the two files
    may order REF/ALT differently, in which case equal indices do not mean the
    same allele.
    """
    return [target.allele_a[i] == gt.allele_a[j] for i, j in block]


def interval(target, pair0, pair1):
    """BED interval spanning two variants, in whatshap's --switch-error-bed convention."""
    return (target.chromosome[pair0[0]],
            int(target.position[pair0[0]]) + 1,
            int(target.position[pair1[0]]) + 1)


@dataclass
class JunctionResult:
    junctions: int = 0
    errors: int = 0
    positions: Optional[List[Tuple[str, int, int]]] = None
    # only filled in when a baseline call set is supplied
    preexisting: int = 0                # pairs the baseline already connected
    no_truth_frame: int = 0             # pairs the truth splits, so unevaluable
    skipped_no_baseline_phase: int = 0  # variants the baseline does not phase
    baseline_blocks: int = 0            # blocks entering the scan
    target_blocks: int = 0              # target blocks those baseline blocks form
    singleton_baseline_blocks: int = 0  # of those, holding one scorable variant
    structural_joins: int = 0           # sum over target blocks of (baseline blocks - 1)

    @property
    def error_rate(self):
        return self.errors / self.junctions if self.junctions else float('nan')


def evaluate_phasing(overlapping_sites, target, gt, ignore_phase_blocks=False):
    """Count switch errors inside blocks that are contiguously phased in both call sets.

    Mirrors ``whatshap compare``: blocks holding a single variant are skipped,
    and a switch is a position where the two phasings stop agreeing (or stop
    disagreeing).  Junctions *between* blocks are not assessed -- for those see
    :func:`evaluate_junctions`.  With ``ignore_phase_blocks`` the blocks are
    concatenated into one chain per chromosome, so the junctions are assessed
    along with everything else.
    """
    blocks, blocks_target, blocks_gt = group_blocks(
        overlapping_sites, target, gt, ignore_phase_blocks=ignore_phase_blocks)

    result = ComparisonResult(switch_positions=[])
    for block in blocks.values():
        if len(block) < 2:
            continue
        result.intersection_blocks += 1
        result.covered_variants += len(block)
        result.assessed_pairs += len(block) - 1

        agree = agreement(block, target, gt)
        same_genotype = [{target.allele_a[i], target.allele_b[i]} ==
                         {gt.allele_a[j], gt.allele_b[j]} for i, j in block]
        result.diff_genotypes += len(block) - sum(same_genotype)

        h0 = ''.join('1' if a else '0' for a in agree)
        h1 = '1' * len(block)
        result.switches += hamming(switch_encoding(h0), switch_encoding(h1))
        switches, flips = compute_switch_flips(h0, h1)
        result.switch_flips = (result.switch_flips[0] + switches, result.switch_flips[1] + flips)
        result.hamming += min(sum(1 for a in agree if not a), sum(1 for a in agree if a))

        for k in range(len(block) - 1):
            if agree[k] != agree[k + 1]:
                result.switch_positions.append(interval(target, block[k], block[k + 1]))

    result.blocks_target = sum(1 for b in blocks_target.values() if len(b) > 1)
    result.covered_target = sum(len(b) for b in blocks_target.values() if len(b) > 1)
    result.blocks_gt = sum(1 for b in blocks_gt.values() if len(b) > 1)
    result.covered_gt = sum(len(b) for b in blocks_gt.values() if len(b) > 1)
    return result


def evaluate_junctions(overlapping_sites, target, gt):
    """Score only the connections *between* phase blocks.

    Walks the variants of each chromosome in order and assesses a pair whenever
    consecutive variants fall into different intersection blocks, i.e. exactly
    the pairs that :func:`evaluate_phasing` cannot see but the chromosome-wide
    comparison adds.  An error means the two blocks were joined in the wrong
    orientation relative to the truth.
    """
    blocks, _, _ = group_blocks(overlapping_sites, target, gt)
    block_of = {}
    for key, block in blocks.items():
        for pair in block:
            block_of[tuple(pair)] = key

    chains: Dict[str, List[tuple]] = defaultdict(list)
    for i, j in overlapping_sites:
        pair = (int(i), int(j))
        if pair in block_of:
            chains[target.chromosome[i]].append(pair)

    result = JunctionResult(positions=[])
    for chain in chains.values():
        agree = agreement(chain, target, gt)
        for k in range(len(chain) - 1):
            if block_of[chain[k]] == block_of[chain[k + 1]]:
                continue
            result.junctions += 1
            if agree[k] != agree[k + 1]:
                result.errors += 1
                result.positions.append(interval(target, chain[k], chain[k + 1]))
    return result


def baseline_block_map(target, baseline, match='strict'):
    """Map each target variant to the phase block it sat in before the new joins.

    Pairing goes through :func:`get_overlapping_sites`, so ``--match`` is
    honoured here exactly as it is against the truth.  Variants the baseline
    leaves unphased get no entry.
    """
    baseline_of = {}
    for i, j in get_overlapping_sites(target, baseline, match=match):
        if baseline.phased[j] and baseline.complete[j]:
            baseline_of[int(i)] = (baseline.chromosome[j], baseline.phase_block[j])
    return baseline_of


def evaluate_new_junctions(overlapping_sites, target, gt, baseline_of):
    """Score only the connections the target established but the baseline did not.

    Walks the comparable variants of each chromosome in position order, keeping
    only those the baseline also phases -- so a variant the baseline cannot
    place is stepped over and the nearest phased variants on either side are
    compared directly.  An adjacent pair is a *new connection* when the target
    puts both variants in one phase block, the baseline puts them in two, and
    the truth phases both within one block to supply a common frame.  It is an
    error when the two disagree in orientation relative to the truth, i.e. the
    baseline blocks were joined the wrong way round.

    Note that new connections are a subset of the pairs
    :func:`evaluate_phasing` assesses *within* blocks, not of the ones
    :func:`evaluate_junctions` finds between them: once the target has merged
    two baseline blocks, both variants carry the same target PS.
    """
    blocks, _, _ = group_blocks(overlapping_sites, target, gt)
    block_of = {}
    for key, block in blocks.items():
        for pair in block:
            block_of[tuple(pair)] = key

    result = JunctionResult(positions=[])
    chains: Dict[str, List[tuple]] = defaultdict(list)
    scored_target: Counter = Counter()
    scored_baseline: Counter = Counter()
    for i, j in overlapping_sites:
        pair = (int(i), int(j))
        if pair not in block_of:
            continue
        if pair[0] not in baseline_of:
            result.skipped_no_baseline_phase += 1
            continue
        # Both counts are taken over exactly the variants scanned below, so that
        # baseline blocks - target blocks equals the number of joins seen.
        scored_target[(target.chromosome[i], target.phase_block[i])] += 1
        scored_baseline[baseline_of[pair[0]]] += 1
        chains[target.chromosome[i]].append(pair)

    result.target_blocks = len(scored_target)
    result.baseline_blocks = len(scored_baseline)
    result.singleton_baseline_blocks = sum(1 for n in scored_baseline.values() if n == 1)

    # How many joins the target made in total, counted over every variant it and
    # the baseline both phase -- independent of whether the truth can score them,
    # so joins lost for want of a comparable truth variant stay visible.
    merged = defaultdict(set)
    for i, block in baseline_of.items():
        if target.phased[i] and target.complete[i]:
            merged[(target.chromosome[i], target.phase_block[i])].add(block)
    result.structural_joins = sum(len(v) - 1 for v in merged.values())

    for chain in chains.values():
        agree = agreement(chain, target, gt)
        for k in range(len(chain) - 1):
            (i0, j0), (i1, j1) = chain[k], chain[k + 1]
            if target.phase_block[i0] != target.phase_block[i1]:
                continue                        # the target claims no connection here
            if baseline_of[i0] == baseline_of[i1]:
                result.preexisting += 1         # the baseline already made it
                continue
            if gt.phase_block[j0] != gt.phase_block[j1]:
                result.no_truth_frame += 1      # no shared frame to judge against
                continue
            result.junctions += 1
            if agree[k] != agree[k + 1]:
                result.errors += 1
                result.positions.append(interval(target, chain[k], chain[k + 1]))
    return result


def print_stat(label, value=None):
    label = label.rjust(LABEL_WIDTH)
    if value is None:
        print(label)
    elif isinstance(value, int):
        print(label + ':', str(value).rjust(COUNT_WIDTH))
    else:
        print(label + ':', str(value).rjust(COUNT_WIDTH))


def percent(numerator, denominator):
    return '--' if denominator == 0 else f'{numerator * 100.0 / denominator:.2f}%'


def report_switches(result, header):
    print()
    print(f'{header}:'.rjust(LABEL_WIDTH), '-' * COUNT_WIDTH)
    print_stat('phased pairs of variants assessed', result.assessed_pairs)
    print_stat('switch errors', result.switches)
    print_stat('switch error rate', percent(result.switches, result.assessed_pairs))
    print_stat('switch/flip decomposition', '{}/{}'.format(*result.switch_flips))
    print_stat('switch/flip rate', percent(sum(result.switch_flips), result.assessed_pairs))
    print_stat('Block-wise Hamming distance', result.hamming)
    print_stat('Block-wise Hamming distance [%]', percent(result.hamming, result.covered_variants))
    print_stat('Different genotypes', result.diff_genotypes)
    print_stat('Different genotypes [%]', percent(result.diff_genotypes, result.covered_variants))


def report(target, gt, overlapping_sites, within, junctions, chromosome_wide,
           new_connections=None, names=('truth', 'query')):
    gt_name, target_name = names
    print('VARIANT COUNTS (heterozygous / all): ')
    for name, vs in ((gt_name, gt), (target_name, target)):
        print(f'{name}:'.rjust(LABEL_WIDTH), str(len(vs)).rjust(COUNT_WIDTH),
              '/', str(vs.n_total).rjust(COUNT_WIDTH))
    het_gt, het_target = set(gt.keys()), set(target.keys())
    for label, het, allk in (('UNION', het_gt | het_target, gt.all_keys | target.all_keys),
                             ('INTERSECTION', het_gt & het_target, gt.all_keys & target.all_keys)):
        print(f'{label}:'.rjust(LABEL_WIDTH), str(len(het)).rjust(COUNT_WIDTH),
              '/', str(len(allk)).rjust(COUNT_WIDTH))
    print_stat('common heterozygous variants', len(overlapping_sites))
    print_stat(f'non-singleton blocks in {gt_name}', within.blocks_gt)
    print_stat('--> covered variants', within.covered_gt)
    print_stat(f'non-singleton blocks in {target_name}', within.blocks_target)
    print_stat('--> covered variants', within.covered_target)
    print_stat('non-singleton intersection blocks', within.intersection_blocks)
    print_stat('--> covered variants', within.covered_variants)

    report_switches(within, 'WITHIN PHASE BLOCKS')

    print()
    print('BETWEEN PHASE BLOCKS:'.rjust(LABEL_WIDTH), '-' * COUNT_WIDTH)
    print_stat('block junctions assessed', junctions.junctions)
    print_stat('junction errors', junctions.errors)
    print_stat('junction error rate', percent(junctions.errors, junctions.junctions))

    if new_connections is not None:
        n = new_connections
        print()
        print('NEW CONNECTIONS (vs baseline):'.rjust(LABEL_WIDTH), '-' * COUNT_WIDTH)
        print_stat('new connections assessed', n.junctions)
        print_stat('connection errors', n.errors)
        print_stat('connection error rate', percent(n.errors, n.junctions))
        print_stat('baseline blocks scanned', n.baseline_blocks)
        print_stat('--> target blocks they form', n.target_blocks)
        print_stat('--> joins between them', n.baseline_blocks - n.target_blocks)
        print_stat('of those, holding one scorable variant', n.singleton_baseline_blocks)
        print_stat('joins the target made over all variants', n.structural_joins)
        print_stat('--> not assessed, no truth variant',
                   n.structural_joins - n.junctions - n.no_truth_frame)
        print_stat('--> not assessed, no truth frame', n.no_truth_frame)
        print_stat('variants skipped, unphased in baseline', n.skipped_no_baseline_phase)

    report_switches(chromosome_wide, 'CHROMOSOME-WIDE (blocks concatenated)')


def write_bed(path, intervals, annotation):
    with open(path, 'w') as out:
        for chrom, start, end in sorted(intervals):
            print(chrom, start, end, annotation, sep='\t', file=out)


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", required=True, help='VCF to evaluate')
    parser.add_argument("--gt_vcf", required=True, help='Ground truth VCF')
    parser.add_argument('--sample', help='Sample to evaluate (default: first sample in each VCF)')
    parser.add_argument('--only_snvs', action='store_true', default=False,
                        help='Ignore everything but SNVs (whatshap compare --only-snvs)')
    parser.add_argument('--keep_duplicate_positions', action='store_true', default=False,
                        help='Keep records sharing a position with the preceding record. '
                             'whatshap drops these, so this makes counts diverge from it.')
    parser.add_argument('--match', choices=('strict', 'alleles'), default='strict',
                        help="How to pair variants across the two files. 'strict' uses "
                             "(position, REF, ALT) like whatshap; 'alleles' pairs any two "
                             "sites carrying the same unordered allele pair.")
    parser.add_argument('--switch_error_bed',
                        help='Write within-block switch error intervals to this BED file')
    parser.add_argument('--junction_error_bed',
                        help='Write misoriented block junctions to this BED file')
    parser.add_argument('--new_connection_bed',
                        help='Write misoriented new connections to this BED file')
    parser.add_argument('--evaluate_sv', action='store_true', default=False,
                        help='Evaluate structural variants instead of SNPs')
    parser.add_argument('--baseline_vcf', required=False,
                        help='Phasing the target VCF was built on, e.g. the pre-methylation '
                             'call set. Given this, the error rate of the connections the '
                             'target newly established is reported.')
    parser.add_argument('--annotations', required=False,
                        help='Bed file with genome annotations, e.g., CDS, '
                             'for which to evaluate the phasing of compound heterozygous variants')
    args = parser.parse_args()

    load_kwargs = dict(sample=args.sample, only_snvs=args.only_snvs,
                       keep_duplicate_positions=args.keep_duplicate_positions)
    target = load_phasing(args.vcf, **load_kwargs)
    gt = load_phasing(args.gt_vcf, **load_kwargs)

    overlapping_sites = get_overlapping_sites(target, gt, match=args.match)
    # switches inside stretches that both files phase contiguously
    within = evaluate_phasing(overlapping_sites, target, gt)
    # orientation of the connections between those stretches
    junctions = evaluate_junctions(overlapping_sites, target, gt)
    # both together: every consecutive pair of comparable variants on a chromosome
    chromosome_wide = evaluate_phasing(overlapping_sites, target, gt, ignore_phase_blocks=True)

    new_connections = None
    if args.baseline_vcf:
        # only the joins the target added on top of the baseline
        baseline = load_phasing(args.baseline_vcf, **load_kwargs)
        baseline_of = baseline_block_map(target, baseline, match=args.match)
        new_connections = evaluate_new_junctions(overlapping_sites, target, gt, baseline_of)

    report(target, gt, overlapping_sites, within, junctions, chromosome_wide, new_connections)

    if args.switch_error_bed:
        write_bed(args.switch_error_bed, within.switch_positions, 'switch')
    if args.junction_error_bed:
        write_bed(args.junction_error_bed, junctions.positions, 'junction')
    if args.new_connection_bed and new_connections is not None:
        write_bed(args.new_connection_bed, new_connections.positions, 'new_connection')

if __name__ == '__main__':
    main(sys.argv[1:])
