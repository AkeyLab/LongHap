#!/usr/bin/env python3
"""Evaluate variant phasing against a ground truth phasing.

The variant counting and switch-error statistics are defined to reproduce
``whatshap compare`` exactly (see the notes on each option below), so results
from this script can be checked against it directly.  The statistics beyond it --
the block junctions, the new connections, the SV and gene tests, and the windowed
Hamming distance -- have no counterpart there.
"""
import argparse
import gzip
import heapq
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from cyvcf2 import VCF
import numpy as np

# whatshap refuses to model sites with this many alleles or more (whatshap.core)
MAX_GENOTYPE_ALLELES = 16

# The fixed-size windows the Hamming distance is additionally scored in, as
# (column prefix, unit, report heading).  See evaluate_windowed_hamming().
WINDOW_SECTIONS = (
    ('winbp', 'bp', 'HAMMING IN FIXED bp WINDOWS'),
    ('winvar', 'variants', 'HAMMING IN FIXED VARIANT WINDOWS'),
)

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
    allele_frequency: np.ndarray      # INFO/AF, NaN where the field is absent
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


def parse_allele_frequency(variant):
    """``INFO/AF`` as a float, NaN when the variant carries no AF at all.

    NaN is how "not in the reference panel" is represented, and it must never
    count as rare: ``np.nan < threshold`` is already False, so the comparison
    does the right thing on its own -- do not "simplify" it into something that
    treats missing as zero.

    ``bcftools norm -m -any`` subsets a Number=A field per split record, so the
    value is normally scalar; where a sequence still arrives the **largest** is
    taken, so a site counts as rare only when every ALT of it is rare.
    """
    af = variant.INFO.get('AF')
    if af is None:
        return float('nan')
    if isinstance(af, (tuple, list)):
        values = [v for v in af if v is not None]
        return max(values) if values else float('nan')
    return float(af)


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
    allele_frequency = []
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
        allele_frequency.append(parse_allele_frequency(variant))
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
        allele_frequency=np.array(allele_frequency, dtype=float),
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


def interval(target, i0, i1):
    """BED interval spanning two variants, in whatshap's --switch-error-bed convention."""
    return (target.chromosome[i0],
            int(target.position[i0]) + 1,
            int(target.position[i1]) + 1)


def scorable_chains(overlapping_sites, target, gt):
    """Comparable variant pairs per chromosome, in position order.

    A pair is comparable when both files phase it with a complete genotype and
    the truth can place it, i.e. it sits in one of the intersection blocks.
    """
    blocks, _, _ = group_blocks(overlapping_sites, target, gt)
    scorable = {tuple(pair) for block in blocks.values() for pair in block}
    chains: Dict[str, List[tuple]] = defaultdict(list)
    for i, j in overlapping_sites:
        pair = (int(i), int(j))
        if pair in scorable:
            chains[target.chromosome[i]].append(pair)
    return chains


def target_block_sizes(chains, target):
    """Scorable variants per target phase block, over the chains of one scan.

    A block holding a single one of them can carry no junction of its own, and a
    call set that strands such variants -- methphaser leaves the first variant of
    every block it merges behind, still carrying its old PS -- fragments the chain
    the scans walk.  Both scans count them, so that fragmentation shows up in the
    report instead of silently moving pairs between the statistics.
    """
    sizes: Counter = Counter()
    for chain in chains.values():
        for i, _ in chain:
            sizes[(target.chromosome[i], int(target.phase_block[i]))] += 1
    return sizes


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
    singleton_target_blocks: int = 0    # target blocks holding one scorable variant
    interleaved_blocks: int = 0         # target blocks whose baseline blocks alternate
    # the same joins scored against the baseline blocks' own boundary variants
    naive_junctions: int = 0
    naive_errors: int = 0
    naive_no_truth_frame: int = 0       # of those, unevaluable at the naive endpoints
    naive_differs: int = 0              # joins where the naive pair is not the scored pair
    # only filled in without a baseline call set
    at_singleton: int = 0               # junctions touching a singleton target block
    errors_at_singleton: int = 0        # of those, the ones that are errors

    @property
    def error_rate(self):
        return self.errors / self.junctions if self.junctions else float('nan')


@dataclass
class SvResult:
    total: int = 0          # non-SNPs in the call set, phased or not
    svs: int = 0            # of those, the ones the truth can score
    connections: int = 0
    errors: int = 0
    correct: int = 0        # anchors agree with each other, the SV agrees with both
    flipped: int = 0        # anchors agree with each other, the SV disagrees with both
    ambiguous: int = 0      # anchors disagree, so a real switch sits at the SV
    one_sided: int = 0      # only one of the two connections is evaluable
    no_anchor: int = 0      # neither side is evaluable
    anchor_is_target: int = 0   # connections whose anchor is itself a target
    by_type: Optional[Counter] = None
    positions: Optional[List[Tuple[str, int, int]]] = None

    @property
    def error_rate(self):
        return self.errors / self.connections if self.connections else float('nan')


@dataclass
class GeneResult:
    # The first three depend only on the call set and the annotation, never on how
    # well anything was phased, so they are identical across tools run on one VCF
    # and give the comparison a fixed denominator.
    genes: int = 0                 # genes carrying at least two het sites, phased or not
    sites: int = 0                 # het sites inside those genes, phased or not
    single_site_genes: int = 0     # one het site, so no connection to make
    sites_scorable: int = 0        # of `sites`, those the truth can score
    connections: int = 0           # consecutive pairs the phasing resolves
    errors: int = 0                # of those, pairs placed cis when truth says trans, or vice versa
    unresolved: int = 0            # consecutive pairs it cannot place, for any reason
    genes_correct: int = 0
    genes_with_error: int = 0
    genes_unresolved: int = 0      # no errors, but at least one pair it cannot place
    per_gene: Optional[List[tuple]] = None
    positions: Optional[List[Tuple[str, int, int]]] = None

    @property
    def error_rate(self):
        return self.errors / self.connections if self.connections else float('nan')


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
                result.switch_positions.append(interval(target, block[k][0], block[k + 1][0]))

    result.blocks_target = sum(1 for b in blocks_target.values() if len(b) > 1)
    result.covered_target = sum(len(b) for b in blocks_target.values() if len(b) > 1)
    result.blocks_gt = sum(1 for b in blocks_gt.values() if len(b) > 1)
    result.covered_gt = sum(len(b) for b in blocks_gt.values() if len(b) > 1)
    return result


@dataclass
class WindowResult:
    """Hamming distance scored in windows of one fixed size."""
    window: int = 0
    unit: str = 'bp'
    windows: int = 0
    covered_variants: int = 0
    hamming: int = 0

    @property
    def hamming_rate(self):
        return self.hamming / self.covered_variants if self.covered_variants else float('nan')


def evaluate_windowed_hamming(overlapping_sites, target, gt, window, unit='bp'):
    """Hamming distance restricted to windows of one fixed size.

    The block-wise Hamming distance of :func:`evaluate_phasing` grows with block
    length: a single switch in the middle of an N-variant block leaves ~N/4
    variants on the wrong haplotype, so a call set that builds longer blocks
    scores worse than one that builds shorter blocks at the same switch density.
    That makes ``within_hamming_rate`` incomparable between tools whose phase
    blocks differ in length.  Scoring every call set in windows of one size
    removes the confounding, and what is left is the local error density.

    ``unit='bp'`` groups on the absolute coordinate (``position // window``), so
    every call set is scored on the same genomic windows and the tiling does not
    depend on where a call set happens to start a block.  ``unit='variants'``
    cuts each block into runs of ``window`` consecutive comparable variants,
    which fixes the denominator instead of the span.

    Windows never cross an intersection block: relative orientation is undefined
    between two blocks, so a window overlapping a block boundary is scored as one
    group per block.  Groups holding fewer than two variants are skipped, as
    single-variant blocks are in :func:`evaluate_phasing`.

    The window has to be chosen well below the *smallest* block N50 being
    compared; once it exceeds the blocks it is clipped by them and the statistic
    decays back into the block-wise one it was meant to replace.
    """
    if unit not in ('bp', 'variants'):
        raise ValueError(f"unit must be 'bp' or 'variants', not {unit!r}")
    blocks, _, _ = group_blocks(overlapping_sites, target, gt)

    result = WindowResult(window=window, unit=unit)
    for block in blocks.values():
        if len(block) < 2:
            continue
        agree = agreement(block, target, gt)
        if unit == 'bp':
            groups = defaultdict(list)
            for k, (i, _) in enumerate(block):
                groups[int(target.position[i]) // window].append(agree[k])
            chunks = groups.values()
        else:
            chunks = (agree[k:k + window] for k in range(0, len(agree), window))
        for chunk in chunks:
            if len(chunk) < 2:
                continue
            result.windows += 1
            result.covered_variants += len(chunk)
            result.hamming += min(sum(1 for a in chunk if a), sum(1 for a in chunk if not a))
    return result


def evaluate_junctions(overlapping_sites, target, gt):
    """Score only the connections *between* phase blocks.

    Walks the variants of each chromosome in order and assesses a pair whenever
    consecutive variants fall into different intersection blocks, i.e. exactly
    the pairs that :func:`evaluate_phasing` cannot see but the chromosome-wide
    comparison adds.  An error means the two blocks were joined in the wrong
    orientation relative to the truth.

    Every such pair is scored, so that the within-block pairs and these together
    partition the chromosome-wide comparison.  Pairs touching a single-variant
    target block are additionally counted on their own: such a block is a call set
    fragmenting its own phasing, and it contributes two junctions where an intact
    block would contribute one.  See :func:`target_block_sizes`.
    """
    blocks, _, _ = group_blocks(overlapping_sites, target, gt)
    block_of = {}
    for key, block in blocks.items():
        for pair in block:
            block_of[tuple(pair)] = key

    chains = scorable_chains(overlapping_sites, target, gt)
    sizes = target_block_sizes(chains, target)

    result = JunctionResult(positions=[])
    for chain in chains.values():
        agree = agreement(chain, target, gt)
        for k in range(len(chain) - 1):
            if block_of[chain[k]] == block_of[chain[k + 1]]:
                continue
            result.junctions += 1
            error = agree[k] != agree[k + 1]
            if any(sizes[(target.chromosome[i], int(target.phase_block[i]))] == 1
                   for i, _ in (chain[k], chain[k + 1])):
                result.at_singleton += 1
                result.errors_at_singleton += error
            if error:
                result.errors += 1
                result.positions.append(interval(target, chain[k][0], chain[k + 1][0]))
    return result


def allele_length(vs, k):
    """Length of the longest allele of a variant, used to size indels and SVs."""
    return max(len(vs.ref[k]), max(len(a) for a in vs.alt[k]))


def nearest_snp_indices(chain, target):
    """For every position in the chain, the nearest SNP position to its left and right.

    Anchoring on SNPs keeps both reference points on variants a phaser places from
    direct read support, so a disagreement is attributable to the variant between them.
    """
    left, right = [None] * len(chain), [None] * len(chain)
    last = None
    for k, (i, _) in enumerate(chain):
        left[k] = last
        if target.variant_type[i] == 'SNP':
            last = k
    last = None
    for k in range(len(chain) - 1, -1, -1):
        right[k] = last
        if target.variant_type[chain[k][0]] == 'SNP':
            last = k
    return left, right


def evaluate_anchored(overlapping_sites, target, gt, is_target):
    """Score how well each variant of one class is placed relative to its flanking SNPs.

    For every scorable target the two connections ``anchor_a - target - anchor_b`` are
    tested against the truth, where the anchors are the nearest SNPs on either side.
    The two connections are judged independently: a connection needs only its own two
    variants to share a target phase set and a truth phase set, so a target at the start
    of a block is still scored against the anchor that follows it.

    ``is_target(i)`` picks the class under test -- non-SNPs for ``--evaluate_sv``, rare
    variants and rare non-SNPs for ``--rare_variants``.  Anchoring on SNPs keeps both
    reference points on variants a phaser places from direct read support whichever
    class is being scored, so the three sections share one reference frame and differ
    only in what sits between the anchors.

    Reads nothing but position, alleles, genotype, AF and PS, so it applies to any
    phasing tool that writes standard phase sets.
    """
    result = SvResult(by_type=Counter(), positions=[])
    # Every target in the call set, phased or not, so the denominator does not move
    # with how much of the call set a given tool managed to phase.
    for i in range(len(target)):
        if is_target(i):
            result.total += 1

    for chain in scorable_chains(overlapping_sites, target, gt).values():
        agree = agreement(chain, target, gt)
        left, right = nearest_snp_indices(chain, target)
        for k, (i, j) in enumerate(chain):
            if not is_target(i):
                continue
            result.svs += 1
            result.by_type[target.variant_type[i]] += 1

            sides = []
            for anchor in (left[k], right[k]):
                if anchor is None:
                    continue
                a, b = chain[anchor], chain[k]
                # each connection stands on its own two variants
                if target.phase_block[a[0]] != target.phase_block[b[0]]:
                    continue
                if gt.phase_block[a[1]] != gt.phase_block[b[1]]:
                    continue
                sides.append(anchor)

            if not sides:
                result.no_anchor += 1
                continue

            errors = 0
            for anchor in sides:
                result.connections += 1
                # an anchor drawn from the class under test cannot separate a
                # misplaced target from a misplaced anchor; counted, not excluded
                if is_target(chain[anchor][0]):
                    result.anchor_is_target += 1
                if agree[anchor] != agree[k]:
                    errors += 1
                    lo, hi = sorted((anchor, k))
                    result.positions.append(interval(target, chain[lo][0], chain[hi][0]))
            result.errors += errors

            if len(sides) == 1:
                result.one_sided += 1
            elif errors == 0:
                result.correct += 1
            elif errors == 2:
                result.flipped += 1
            else:
                # the two anchors disagree with each other, so a switch genuinely lies
                # between them and the SV's own placement cannot be singled out
                result.ambiguous += 1
    return result


def evaluate_svs(overlapping_sites, target, gt, min_sv_length=0):
    """Score how well each non-SNP is placed relative to its flanking SNPs."""
    def is_sv(i):
        return (target.variant_type[i] != 'SNP'
                and not (min_sv_length and allele_length(target, i) < min_sv_length))
    return evaluate_anchored(overlapping_sites, target, gt, is_sv)


def evaluate_rare(overlapping_sites, target, gt, max_af, min_sv_length=0, svs_only=False):
    """Score how well each rare variant is placed relative to its flanking SNPs.

    ``svs_only`` narrows the class to rare non-SNPs, honouring ``min_sv_length`` so the
    two SV sections stay comparable.  A variant with no ``INFO/AF`` is never a target --
    ``np.nan < max_af`` is False -- but stays available as an anchor, since not being in
    the reference panel says nothing about how well the phaser placed it.
    """
    def is_rare(i):
        if not target.allele_frequency[i] < max_af:
            return False
        if not svs_only:
            return True
        return (target.variant_type[i] != 'SNP'
                and not (min_sv_length and allele_length(target, i) < min_sv_length))
    return evaluate_anchored(overlapping_sites, target, gt, is_rare)


def load_annotations(path, chromosomes=None):
    """Read a 4-column BED (chrom, start, end, gene) into per-chromosome exon lists.

    One gene occupies many rows -- one per exon, and one per transcript isoform,
    so rows belonging to the same gene routinely overlap each other.  Rows on
    chromosomes absent from the VCFs are dropped, which keeps a genome-wide
    annotation cheap when only one chromosome is being evaluated.
    """
    exons: Dict[str, List[tuple]] = defaultdict(list)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                raise SystemExit(f'{path}:{line_no}: expected 4 columns '
                                 f'(chrom, start, end, gene), found {len(fields)}')
            chrom = fields[0]
            if chromosomes is not None and chrom not in chromosomes:
                continue
            exons[chrom].append((int(fields[1]), int(fields[2]), fields[3]))
    for chrom in exons:
        exons[chrom].sort()
    return exons


def genes_by_variant(indices, target, chrom_exons):
    """Indices of the variants overlapping each gene's exons, in position order.

    Sweeps variants and exons together, keeping the exons that still reach the
    current variant in a heap ordered by end.  A variant is recorded once per
    gene however many of that gene's exons it falls in, and a variant is
    compared by its full reference span so an indel reaching into an exon counts.
    """
    by_gene: Dict[str, List[int]] = defaultdict(list)
    e = 0
    active: List[tuple] = []
    for i in indices:
        start = int(target.position[i])
        end = start + len(target.ref[i])
        while e < len(chrom_exons) and chrom_exons[e][0] < end:
            heapq.heappush(active, (chrom_exons[e][1], chrom_exons[e][2]))
            e += 1
        while active and active[0][0] <= start:
            heapq.heappop(active)
        for gene in {gene for _, gene in active}:
            by_gene[gene].append(i)
    return by_gene


def evaluate_genes(overlapping_sites, target, gt, exons):
    """Score the phasing of heterozygous sites that share a gene.

    Whether two heterozygous sites in one gene are in cis or in trans is what
    decides a compound-heterozygous call, so for each gene the sites overlapping
    any of its exons are taken in position order and every consecutive pair is
    checked against the truth.

    Genes and sites are counted over *every* heterozygous site, phased or not.
    That keeps the denominator a property of the call set and the annotation
    alone, so two phasing tools run on the same VCF are compared over the same
    genes; how much each of them managed to phase shows up in the split between
    resolved and unresolved connections instead.  A pair is unresolved whenever
    it carries no cis/trans claim -- either side unphased or unmatched in the
    truth, or the two sides in different phase sets.
    """
    # truth counterpart of every site the comparison can score
    scorable = {}
    blocks, _, _ = group_blocks(overlapping_sites, target, gt)
    for block in blocks.values():
        for i, j in block:
            scorable[int(i)] = int(j)

    by_chrom: Dict[str, List[int]] = defaultdict(list)
    for i in range(len(target)):
        by_chrom[target.chromosome[i]].append(i)

    result = GeneResult(per_gene=[], positions=[])
    for chrom, indices in by_chrom.items():
        genes = genes_by_variant(indices, target, exons.get(chrom, ()))
        for gene, sites in sorted(genes.items()):
            sites = sorted(set(sites))
            if len(sites) < 2:
                result.single_site_genes += 1
                continue
            result.genes += 1
            result.sites += len(sites)
            result.sites_scorable += sum(1 for i in sites if i in scorable)

            connections = errors = unresolved = 0
            for m in range(len(sites) - 1):
                a, b = sites[m], sites[m + 1]
                if a not in scorable or b not in scorable:
                    unresolved += 1
                    continue
                ja, jb = scorable[a], scorable[b]
                if (target.phase_block[a] != target.phase_block[b] or
                        gt.phase_block[ja] != gt.phase_block[jb]):
                    unresolved += 1
                    continue
                connections += 1
                if (target.allele_a[a] == gt.allele_a[ja]) != (target.allele_a[b] == gt.allele_a[jb]):
                    errors += 1
                    result.positions.append(interval(target, a, b))

            result.connections += connections
            result.errors += errors
            result.unresolved += unresolved
            if errors:
                result.genes_with_error += 1
            elif unresolved:
                result.genes_unresolved += 1
            else:
                result.genes_correct += 1
            result.per_gene.append((chrom, gene, len(sites), connections, errors, unresolved))
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

    Walks each *target phase block* in position order, keeping only the variants
    the baseline also phases -- so a variant the baseline cannot place is stepped
    over and the nearest phased variants on either side are compared directly.  A
    consecutive pair within one target block is a *new connection* when the
    baseline puts the two variants in different blocks and the truth phases both
    within one block to supply a common frame.  It is an error when the two
    disagree in orientation relative to the truth, i.e. the baseline blocks were
    joined the wrong way round.

    Walking per target block rather than along the whole chromosome is what makes
    the scan robust to a call set that strands variants at the old PS: methphaser
    relabels and flips a merged block but leaves its first variant behind, so the
    two variants either side of the join are *not* adjacent on the chromosome.
    The stranded variant lands in a block of its own, contributes nothing, and the
    join is scored between the variants that actually carry the decision.  It must
    not be folded back into the surrounding block instead: it still holds the
    un-flipped genotype, so scoring against it would test the base phaser's
    arbitrary cross-block orientation rather than the join.

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
    chains: Dict[tuple, List[tuple]] = defaultdict(list)
    scored_baseline: Counter = Counter()
    # First and last scorable variant of each baseline block, taken over every
    # variant rather than over the per-target-block chains below.  A variant the
    # target stranded at its old PS is missing from those chains but present
    # here, which is what lets the naive statistic reach it.  overlapping_sites
    # is in VCF order, so first come first and last wins.
    first_of_baseline: Dict[tuple, tuple] = {}
    last_of_baseline: Dict[tuple, tuple] = {}
    for i, j in overlapping_sites:
        pair = (int(i), int(j))
        if pair not in block_of:
            continue
        if pair[0] not in baseline_of:
            result.skipped_no_baseline_phase += 1
            continue
        scored_baseline[baseline_of[pair[0]]] += 1
        first_of_baseline.setdefault(baseline_of[pair[0]], pair)
        last_of_baseline[baseline_of[pair[0]]] = pair
        chains[(target.chromosome[i], int(target.phase_block[i]))].append(pair)

    # Both counts are taken over exactly the variants scanned below.  Their
    # difference is the number of joins only when the target left no block
    # fragmented; singleton_target_blocks is exactly what accounts for the gap.
    result.target_blocks = len(chains)
    result.baseline_blocks = len(scored_baseline)
    result.singleton_baseline_blocks = sum(1 for n in scored_baseline.values() if n == 1)
    result.singleton_target_blocks = sum(1 for chain in chains.values() if len(chain) == 1)

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
        changes = 0
        for k in range(len(chain) - 1):
            (i0, j0), (i1, j1) = chain[k], chain[k + 1]
            if baseline_of[i0] == baseline_of[i1]:
                result.preexisting += 1         # the baseline already made it
                continue
            changes += 1
            if gt.phase_block[j0] != gt.phase_block[j1]:
                result.no_truth_frame += 1      # no shared frame to judge against
                continue
            result.junctions += 1
            if agree[k] != agree[k + 1]:
                result.errors += 1
                result.positions.append(interval(target, chain[k][0], chain[k + 1][0]))

            # The same join, scored where a reader would naively look for it: at
            # the two baseline blocks' own boundary variants.  For a target that
            # relabelled every variant those are the two just scored and the
            # counts cannot diverge; where one was stranded this reads it, which
            # is the whole point of reporting the two rates side by side.
            naive_a = last_of_baseline[baseline_of[i0]]
            naive_b = first_of_baseline[baseline_of[i1]]
            if (naive_a, naive_b) != (chain[k], chain[k + 1]):
                result.naive_differs += 1
            if gt.phase_block[naive_a[1]] != gt.phase_block[naive_b[1]]:
                result.naive_no_truth_frame += 1
            else:
                result.naive_junctions += 1
                agree_a = target.allele_a[naive_a[0]] == gt.allele_a[naive_a[1]]
                agree_b = target.allele_a[naive_b[0]] == gt.allele_a[naive_b[1]]
                if agree_a != agree_b:
                    result.naive_errors += 1
        # A block whose baseline blocks alternate along it would be scanned as
        # more joins than it made.  Never seen in practice; counted so that it
        # cannot pass unnoticed if it starts happening.
        if changes > len({baseline_of[i] for i, _ in chain}) - 1:
            result.interleaved_blocks += 1
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


def fraction(numerator, denominator):
    """Machine-readable counterpart of percent(): a float, or '' when undefined."""
    return '' if denominator == 0 else numerator / denominator


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


def report_windows(result, header):
    print()
    print(f'{header}:'.rjust(LABEL_WIDTH), '-' * COUNT_WIDTH)
    print_stat(f'window size [{result.unit}]', result.window)
    print_stat('windows scored', result.windows)
    print_stat('--> covered variants', result.covered_variants)
    print_stat('Windowed Hamming distance', result.hamming)
    print_stat('Windowed Hamming distance [%]',
               percent(result.hamming, result.covered_variants))


def report_anchored(result, header, noun):
    print()
    print(f'{header}:'.rjust(LABEL_WIDTH), '-' * COUNT_WIDTH)
    print_stat(f'{noun} in the call set', result.total)
    print_stat('--> evaluated', result.svs)
    print_stat('anchored connections', result.connections)
    print_stat('connection errors', result.errors)
    print_stat('connection error rate', percent(result.errors, result.connections))
    print_stat('--> correct', result.correct)
    print_stat('--> flipped', result.flipped)
    print_stat('--> ambiguous, anchors disagree', result.ambiguous)
    print_stat('--> anchored on one side only', result.one_sided)
    print_stat('skipped, no anchoring SNP in block', result.no_anchor)
    print_stat('connections whose anchor is also a target', result.anchor_is_target)
    for name, count in sorted(result.by_type.items()):
        print_stat(f'by type: {name}', count)


def report(target, gt, overlapping_sites, within, junctions, chromosome_wide,
           new_connections=None, sv=None, genes=None, names=('truth', 'query'),
           windowed=None, rare=None, rare_sv=None):
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

    for prefix, _, header in WINDOW_SECTIONS:
        result = (windowed or {}).get(prefix)
        if result is not None:
            report_windows(result, header)

    if sv is not None:
        report_anchored(sv, 'SV PHASING (SNP-anchored)', 'non-SNPs')
    if rare is not None:
        report_anchored(rare, 'RARE VARIANT PHASING (SNP-anchored)', 'rare variants')
    if rare_sv is not None:
        report_anchored(rare_sv, 'RARE SV PHASING (SNP-anchored)', 'rare non-SNPs')

    print()
    print('BETWEEN PHASE BLOCKS:'.rjust(LABEL_WIDTH), '-' * COUNT_WIDTH)
    print_stat('block junctions assessed', junctions.junctions)
    print_stat('junction errors', junctions.errors)
    print_stat('junction error rate', percent(junctions.errors, junctions.junctions))
    print_stat('--> at a singleton target block', junctions.at_singleton)
    print_stat('--> of those, errors', junctions.errors_at_singleton)
    print_stat('excluding those, junctions', junctions.junctions - junctions.at_singleton)
    print_stat('--> errors', junctions.errors - junctions.errors_at_singleton)
    print_stat('--> error rate',
               percent(junctions.errors - junctions.errors_at_singleton,
                       junctions.junctions - junctions.at_singleton))

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
        print_stat('target blocks holding one scorable variant', n.singleton_target_blocks)
        print_stat('target blocks with interleaved baseline blocks', n.interleaved_blocks)
        print_stat('joins the target made over all variants', n.structural_joins)
        print_stat('--> not assessed, no truth variant',
                   n.structural_joins - n.junctions - n.no_truth_frame)
        print_stat('--> not assessed, no truth frame', n.no_truth_frame)
        print_stat('variants skipped, unphased in baseline', n.skipped_no_baseline_phase)
        print_stat('naive: connections assessed', n.naive_junctions)
        print_stat('naive: connection errors', n.naive_errors)
        print_stat('naive: connection error rate', percent(n.naive_errors, n.naive_junctions))
        print_stat('--> not assessed, no truth frame', n.naive_no_truth_frame)
        print_stat('--> joins where the naive pair differs', n.naive_differs)

    if genes is not None:
        print()
        print('COMPOUND HETEROZYGOSITY (per gene):'.rjust(LABEL_WIDTH), '-' * COUNT_WIDTH)
        print_stat('genes with >=2 heterozygous sites', genes.genes)
        print_stat('--> heterozygous sites in them', genes.sites)
        print_stat('--> of those, scorable', genes.sites_scorable)
        print_stat('connections evaluated', genes.connections)
        print_stat('connection errors', genes.errors)
        print_stat('connection error rate', percent(genes.errors, genes.connections))
        print_stat('connections unresolved', genes.unresolved)
        print_stat('--> genes entirely correct', genes.genes_correct)
        print_stat('--> genes with a wrong connection', genes.genes_with_error)
        print_stat('--> genes only partly resolved', genes.genes_unresolved)
        print_stat('genes with a single heterozygous site', genes.single_site_genes)

    report_switches(chromosome_wide, 'CHROMOSOME-WIDE (blocks concatenated)')


def switch_columns(prefix, result):
    """The ten switch statistics, named under one section prefix."""
    if result is None:
        keys = ('assessed_pairs', 'switches', 'switch_rate', 'switchflip_switches',
                'switchflip_flips', 'switchflip_rate', 'hamming', 'hamming_rate',
                'diff_genotypes', 'diff_genotypes_rate')
        return {f'{prefix}_{k}': '' for k in keys}
    switches, flips = result.switch_flips
    return {
        f'{prefix}_assessed_pairs': result.assessed_pairs,
        f'{prefix}_switches': result.switches,
        f'{prefix}_switch_rate': fraction(result.switches, result.assessed_pairs),
        f'{prefix}_switchflip_switches': switches,
        f'{prefix}_switchflip_flips': flips,
        f'{prefix}_switchflip_rate': fraction(switches + flips, result.assessed_pairs),
        f'{prefix}_hamming': result.hamming,
        f'{prefix}_hamming_rate': fraction(result.hamming, result.covered_variants),
        f'{prefix}_diff_genotypes': result.diff_genotypes,
        f'{prefix}_diff_genotypes_rate': fraction(result.diff_genotypes, result.covered_variants),
    }


def window_columns(prefix, result):
    """The five windowed-Hamming statistics, named under one section prefix."""
    if result is None:
        keys = ('window', 'windows', 'covered_variants', 'hamming', 'hamming_rate')
        return {f'{prefix}_{k}': '' for k in keys}
    return {
        f'{prefix}_window': result.window,
        f'{prefix}_windows': result.windows,
        f'{prefix}_covered_variants': result.covered_variants,
        f'{prefix}_hamming': result.hamming,
        f'{prefix}_hamming_rate': fraction(result.hamming, result.covered_variants),
    }


SV_TYPES = ('DEL', 'INS', 'Multi', 'Other')


def anchored_columns(prefix, result):
    """The statistics of one anchored section, named under one section prefix.

    ``sv`` keeps the column names it has always had, so rows written before and
    after the rare sections existed still concatenate.
    """
    keys = ('total', 'evaluated', 'connections', 'errors', 'error_rate', 'correct',
            'flipped', 'ambiguous', 'one_sided', 'no_anchor', 'anchor_is_target')
    if result is None:
        return {f'{prefix}_{k}': '' for k in keys + SV_TYPES}
    out = {
        f'{prefix}_total': result.total,
        f'{prefix}_evaluated': result.svs,
        f'{prefix}_connections': result.connections,
        f'{prefix}_errors': result.errors,
        f'{prefix}_error_rate': fraction(result.errors, result.connections),
        f'{prefix}_correct': result.correct,
        f'{prefix}_flipped': result.flipped,
        f'{prefix}_ambiguous': result.ambiguous,
        f'{prefix}_one_sided': result.one_sided,
        f'{prefix}_no_anchor': result.no_anchor,
        f'{prefix}_anchor_is_target': result.anchor_is_target,
    }
    out.update({f'{prefix}_{name}': result.by_type.get(name, 0) for name in SV_TYPES})
    return out


def collect_summary(args, target, gt, overlapping_sites, within, junctions, chromosome_wide,
                    new_connections=None, sv=None, genes=None, windowed=None,
                    rare=None, rare_sv=None):
    """Every statistic the report prints, as one flat row.

    Built from the same result objects report() reads, so the two cannot disagree.
    Tests that did not run still contribute their columns, with empty values, so the
    header stays identical between runs and the rows concatenate.
    """
    het_gt, het_target = set(gt.keys()), set(target.keys())
    row = {
        'label': args.label or '',
        'sample': args.sample or '',
        'chromosomes': ','.join(dict.fromkeys(target.chromosome.tolist())),
        'vcf': args.vcf,
        'gt_vcf': args.gt_vcf,
        'baseline_vcf': args.baseline_vcf or '',
        'annotations': args.annotations or '',
        'match': args.match,
        'only_snvs': int(args.only_snvs),
        'keep_duplicate_positions': int(args.keep_duplicate_positions),
        'min_sv_length': args.min_sv_length,

        'truth_het': len(gt),
        'truth_all': gt.n_total,
        'query_het': len(target),
        'query_all': target.n_total,
        'union_het': len(het_gt | het_target),
        'union_all': len(gt.all_keys | target.all_keys),
        'intersection_het': len(het_gt & het_target),
        'intersection_all': len(gt.all_keys & target.all_keys),
        'common_het': len(overlapping_sites),
        'truth_blocks': within.blocks_gt,
        'truth_covered': within.covered_gt,
        'query_blocks': within.blocks_target,
        'query_covered': within.covered_target,
        'intersection_blocks': within.intersection_blocks,
        'intersection_covered': within.covered_variants,
    }
    row.update(switch_columns('within', within))
    row.update(switch_columns('chrwide', chromosome_wide))
    for prefix, _, _ in WINDOW_SECTIONS:
        row.update(window_columns(prefix, (windowed or {}).get(prefix)))

    clean_junctions = junctions.junctions - junctions.at_singleton
    clean_errors = junctions.errors - junctions.errors_at_singleton
    row.update({
        'junction_assessed': junctions.junctions,
        'junction_errors': junctions.errors,
        'junction_rate': fraction(junctions.errors, junctions.junctions),
        'junction_at_singleton': junctions.at_singleton,
        'junction_errors_at_singleton': junctions.errors_at_singleton,
        'junction_assessed_intact': clean_junctions,
        'junction_errors_intact': clean_errors,
        'junction_rate_intact': fraction(clean_errors, clean_junctions),
    })

    if new_connections is None:
        row.update({f'newconn_{k}': '' for k in (
            'assessed', 'errors', 'rate', 'baseline_blocks', 'target_blocks', 'joins_between',
            'singleton_baseline_blocks', 'singleton_target_blocks', 'interleaved_blocks',
            'joins_total', 'not_assessed_no_truth_variant',
            'not_assessed_no_truth_frame', 'skipped_unphased_baseline',
            'naive_assessed', 'naive_errors', 'naive_rate',
            'naive_no_truth_frame', 'naive_differs')})
    else:
        n = new_connections
        row.update({
            'newconn_assessed': n.junctions,
            'newconn_errors': n.errors,
            'newconn_rate': fraction(n.errors, n.junctions),
            'newconn_baseline_blocks': n.baseline_blocks,
            'newconn_target_blocks': n.target_blocks,
            'newconn_joins_between': n.baseline_blocks - n.target_blocks,
            'newconn_singleton_baseline_blocks': n.singleton_baseline_blocks,
            'newconn_singleton_target_blocks': n.singleton_target_blocks,
            'newconn_interleaved_blocks': n.interleaved_blocks,
            'newconn_joins_total': n.structural_joins,
            'newconn_not_assessed_no_truth_variant':
                n.structural_joins - n.junctions - n.no_truth_frame,
            'newconn_not_assessed_no_truth_frame': n.no_truth_frame,
            'newconn_skipped_unphased_baseline': n.skipped_no_baseline_phase,
            'newconn_naive_assessed': n.naive_junctions,
            'newconn_naive_errors': n.naive_errors,
            'newconn_naive_rate': fraction(n.naive_errors, n.naive_junctions),
            'newconn_naive_no_truth_frame': n.naive_no_truth_frame,
            'newconn_naive_differs': n.naive_differs,
        })

    row.update(anchored_columns('sv', sv))
    row['max_af'] = getattr(args, 'max_af', '') if rare is not None else ''
    row.update(anchored_columns('rare', rare))
    row.update(anchored_columns('raresv', rare_sv))

    if genes is None:
        row.update({f'gene_{k}': '' for k in (
            'genes', 'sites', 'sites_scorable', 'single_site_genes', 'connections', 'errors',
            'error_rate', 'unresolved', 'correct', 'with_error', 'partly_resolved')})
    else:
        row.update({
            'gene_genes': genes.genes,
            'gene_sites': genes.sites,
            'gene_sites_scorable': genes.sites_scorable,
            'gene_single_site_genes': genes.single_site_genes,
            'gene_connections': genes.connections,
            'gene_errors': genes.errors,
            'gene_error_rate': fraction(genes.errors, genes.connections),
            'gene_unresolved': genes.unresolved,
            'gene_correct': genes.genes_correct,
            'gene_with_error': genes.genes_with_error,
            'gene_partly_resolved': genes.genes_unresolved,
        })
    return row


def write_summary_tsv(path, row):
    """Write the header and the single summary row for this run."""
    with open(path, 'w') as out:
        print(*row.keys(), sep='\t', file=out)
        print(*row.values(), sep='\t', file=out)


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
    parser.add_argument('--hamming_window', type=int, default=100000, metavar='BP',
                        help='Additionally score the Hamming distance in windows of this many '
                             'bp, which unlike the block-wise Hamming distance does not grow '
                             'with phase block length and so is comparable between call sets '
                             'whose blocks differ in length. Pick it well below the smallest '
                             'block N50 being compared. 0 turns the test off. (default: '
                             '%(default)s)')
    parser.add_argument('--hamming_window_variants', type=int, default=100, metavar='N',
                        help='The same test with the window fixed at this many consecutive '
                             'comparable variants instead of a span. 0 turns it off. '
                             '(default: %(default)s)')
    parser.add_argument('--switch_error_bed',
                        help='Write within-block switch error intervals to this BED file')
    parser.add_argument('--junction_error_bed',
                        help='Write misoriented block junctions to this BED file')
    parser.add_argument('--new_connection_bed',
                        help='Write misoriented new connections to this BED file')
    parser.add_argument('--sv_error_bed',
                        help='Write SV-to-anchor connections that disagree to this BED file')
    parser.add_argument('--evaluate_sv', action='store_true', default=False,
                        help='Additionally score how each non-SNP is placed relative to the '
                             'nearest SNP on either side. Works for any phasing tool that '
                             'writes FORMAT/PS (HP tags are not read).')
    parser.add_argument('--rare_variants', action='store_true', default=False,
                        help='Read INFO/AF and additionally score how each rare variant, and '
                             'each rare non-SNP, is placed relative to the nearest SNP on '
                             'either side -- the same test --evaluate_sv applies to non-SNPs. '
                             'Needs a VCF carrying every variant, not one already filtered to '
                             'the rare ones, or there are no common variants left to anchor on.')
    parser.add_argument('--max_af', type=float, default=0.01, metavar='AF',
                        help='With --rare_variants, the allele frequency below which a variant '
                             'counts as rare. A variant with no INFO/AF is never a target, but '
                             'stays available as an anchor (default: %(default)s)')
    parser.add_argument('--rare_error_bed',
                        help='Write rare-variant connections that disagree to this BED file')
    parser.add_argument('--rare_sv_error_bed',
                        help='Write rare-SV connections that disagree to this BED file')
    parser.add_argument('--min_sv_length', type=int, default=0,
                        help='With --evaluate_sv, only score non-SNPs whose longest allele is at '
                             'least this many bp (default: all non-SNPs)')
    parser.add_argument('--baseline_vcf', required=False,
                        help='Phasing the target VCF was built on, e.g. the pre-methylation '
                             'call set. Given this, the error rate of the connections the '
                             'target newly established is reported.')
    parser.add_argument('--annotations', required=False,
                        help='4-column BED (chrom, start, end, gene) of exons; one gene may span '
                             'many rows. Scores whether the heterozygous sites sharing a gene are '
                             'put in the right cis/trans arrangement, which is what a compound '
                             'heterozygous call rests on. May be genome-wide.')
    parser.add_argument('--summary_tsv',
                        help='Write every summary statistic as one row to this TSV. The header is '
                             'the same whichever tests were run, so rows from different runs '
                             'concatenate. Rates are fractions, empty when undefined.')
    parser.add_argument('--label', default='',
                        help='Free-form run identifier carried in the --summary_tsv row')
    parser.add_argument('--gene_tsv',
                        help='With --annotations, write per-gene results to this TSV')
    parser.add_argument('--gene_error_bed',
                        help='With --annotations, write wrongly connected within-gene pairs here')
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
    # the block-length-independent counterpart of the block-wise Hamming distance
    sizes = {'winbp': args.hamming_window, 'winvar': args.hamming_window_variants}
    windowed = {prefix: evaluate_windowed_hamming(overlapping_sites, target, gt,
                                                  sizes[prefix], unit)
                for prefix, unit, _ in WINDOW_SECTIONS if sizes[prefix] > 0}

    new_connections = None
    if args.baseline_vcf:
        # only the joins the target added on top of the baseline
        baseline = load_phasing(args.baseline_vcf, **load_kwargs)
        baseline_of = baseline_block_map(target, baseline, match=args.match)
        new_connections = evaluate_new_junctions(overlapping_sites, target, gt, baseline_of)

    sv = None
    if args.evaluate_sv:
        sv = evaluate_svs(overlapping_sites, target, gt, min_sv_length=args.min_sv_length)

    rare = rare_sv = None
    if args.rare_variants:
        rare = evaluate_rare(overlapping_sites, target, gt, args.max_af)
        rare_sv = evaluate_rare(overlapping_sites, target, gt, args.max_af,
                                min_sv_length=args.min_sv_length, svs_only=True)

    genes = None
    if args.annotations:
        exons = load_annotations(args.annotations, set(target.chromosome.tolist()))
        genes = evaluate_genes(overlapping_sites, target, gt, exons)

    report(target, gt, overlapping_sites, within, junctions, chromosome_wide,
           new_connections, sv, genes, windowed=windowed, rare=rare, rare_sv=rare_sv)

    if args.summary_tsv:
        write_summary_tsv(args.summary_tsv,
                          collect_summary(args, target, gt, overlapping_sites, within, junctions,
                                          chromosome_wide, new_connections, sv, genes, windowed,
                                          rare, rare_sv))

    if args.switch_error_bed:
        write_bed(args.switch_error_bed, within.switch_positions, 'switch')
    if args.junction_error_bed:
        write_bed(args.junction_error_bed, junctions.positions, 'junction')
    if args.new_connection_bed and new_connections is not None:
        write_bed(args.new_connection_bed, new_connections.positions, 'new_connection')
    if args.sv_error_bed and sv is not None:
        write_bed(args.sv_error_bed, sv.positions, 'sv_connection')
    if args.rare_error_bed and rare is not None:
        write_bed(args.rare_error_bed, rare.positions, 'rare_connection')
    if args.rare_sv_error_bed and rare_sv is not None:
        write_bed(args.rare_sv_error_bed, rare_sv.positions, 'rare_sv_connection')
    if args.gene_error_bed and genes is not None:
        write_bed(args.gene_error_bed, genes.positions, 'gene_connection')
    if args.gene_tsv and genes is not None:
        with open(args.gene_tsv, 'w') as out:
            print('chrom', 'gene', 'het_sites', 'connections', 'errors', 'unresolved',
                  sep='\t', file=out)
            for row in genes.per_gene:
                print(*row, sep='\t', file=out)


if __name__ == '__main__':
    main(sys.argv[1:])
