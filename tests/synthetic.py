"""Synthetic data builders for the longhap tests.

Kept out of ``conftest.py`` on purpose: the bare module name ``conftest``
collides with ``benchmarking/tests/conftest.py`` when both suites are collected
in one run, so test modules import their helpers from here instead.

The central idea is that a test declares a *truth* -- a list of heterozygous
variants and, for each, which allele sits on haplotype 0 and which on
haplotype 1.  ``simulate_read`` renders a read from one of those two haplotypes,
producing the query sequence and the CIGAR a real aligner would emit, so a test
can assert that longhap recovers the haplotype it was handed.
"""
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pysam

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

CHROM = 'chr1'


# --------------------------------------------------------------------------- #
# reference
# --------------------------------------------------------------------------- #
def make_sequence(length, seed=0):
    """Deterministic pseudo-random sequence with no homopolymer longer than 3.

    Long homopolymers change which realignment penalties longhap picks, so the
    default backdrop deliberately avoids them; a test that wants one inserts it.
    """
    rng = random.Random(seed)
    seq = []
    for _ in range(length):
        base = rng.choice('ACGT')
        while len(seq) >= 3 and seq[-1] == seq[-2] == seq[-3] == base:
            base = rng.choice('ACGT')
        seq.append(base)
    return ''.join(seq)


# --------------------------------------------------------------------------- #
# variants
# --------------------------------------------------------------------------- #
@dataclass
class Variant:
    """One heterozygous site.

    ``pos`` is 1-based (VCF convention).  ``hap`` gives the allele index carried
    by haplotype 0 and haplotype 1, indexing into ``[ref] + alt``, so ``(0, 1)``
    means haplotype 0 is REF and haplotype 1 is ALT.
    """
    pos: int
    ref: str
    alt: list
    hap: tuple = (0, 1)
    vid: str = '.'

    @property
    def alleles(self):
        return [self.ref] + list(self.alt)

    @property
    def end(self):
        """0-based exclusive end of the REF allele."""
        return self.pos - 1 + len(self.ref)


def snv(pos, ref, alt, hap=(0, 1), vid='.'):
    return Variant(pos, ref, [alt], hap, vid)


def ref_snv(ref_seq, pos, hap=(0, 1), vid='.'):
    """A SNV at 1-based ``pos`` whose REF matches ``ref_seq`` and ALT differs."""
    base = ref_seq[pos - 1]
    return Variant(pos, base, ['ACGT'['ACGT'.index(base) - 1]], hap, vid)


def ref_insertion(ref_seq, pos, inserted='TTAG', hap=(0, 1), vid='.'):
    """A left-anchored insertion at 1-based ``pos``."""
    anchor = ref_seq[pos - 1]
    return Variant(pos, anchor, [anchor + inserted], hap, vid)


def ref_deletion(ref_seq, pos, length=4, hap=(0, 1), vid='.'):
    """A left-anchored deletion of ``length`` bp after the anchor at ``pos``."""
    return Variant(pos, ref_seq[pos - 1:pos - 1 + length + 1], [ref_seq[pos - 1]],
                   hap, vid)


def alternating_snvs(ref_seq, positions, pattern=None):
    """SNVs at ``positions`` with an explicit per-site haplotype assignment.

    ``pattern`` is a string of '0'/'1': '0' means haplotype 0 carries REF, '1'
    means haplotype 0 carries ALT.  Defaults to all-'0', i.e. haplotype 0 is the
    reference haplotype throughout.
    """
    pattern = pattern or '0' * len(positions)
    out = []
    for pos, bit in zip(positions, pattern):
        hap = (0, 1) if bit == '0' else (1, 0)
        out.append(ref_snv(ref_seq, pos, hap=hap))
    return out


VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID={chrom},length={length}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


# --------------------------------------------------------------------------- #
# reads
# --------------------------------------------------------------------------- #
def simulate_read(ref_seq, start, end, variants, hap):
    """Render one read of haplotype ``hap`` spanning ``[start, end)`` (0-based).

    Returns ``(query_sequence, cigartuples)``.  Variants are applied
    left-anchored, exactly as VCF spells them, so an insertion becomes M then I
    and a deletion becomes M then D.  Only variants fully inside the span are
    applied; a variant straddling an edge is left as reference.
    """
    query = []
    cigar = []

    def add(op, n):
        if n <= 0:
            return
        if cigar and cigar[-1][0] == op:
            cigar[-1] = (op, cigar[-1][1] + n)
        else:
            cigar.append((op, n))

    pos = start
    for var in sorted(variants, key=lambda v: v.pos):
        vstart = var.pos - 1
        if vstart < pos or var.end > end:
            continue
        add(0, vstart - pos)
        query.append(ref_seq[pos:vstart])

        allele = var.alleles[var.hap[hap]]
        ref = var.ref
        shared = min(len(ref), len(allele))
        add(0, shared)
        query.append(allele[:shared])
        if len(allele) > len(ref):            # insertion
            add(1, len(allele) - len(ref))
            query.append(allele[shared:])
        elif len(ref) > len(allele):          # deletion
            add(2, len(ref) - len(allele))
        pos = var.end

    add(0, end - pos)
    query.append(ref_seq[pos:end])
    return ''.join(query), [tuple(c) for c in cigar]


@dataclass
class ReadSpec:
    """One simulated alignment."""
    start: int
    end: int
    hap: int
    name: str = None
    mapq: int = 60
    base_quality: int = 40
    is_reverse: bool = False
    is_supplementary: bool = False
    is_secondary: bool = False
    flags_extra: int = 0
    methylation: dict = field(default=None)  # {ref_pos_0based: probability}


def _add_methylation(aln, seq, spec, ref_seq):
    """Attach MM/ML tags for the C bases named in ``spec.methylation``.

    Keys are 0-based reference positions; the value is P(methylated).  Positions
    whose read base is not a C are skipped, which keeps a test honest about what
    a real basecaller would emit.
    """
    pairs = []
    for ref_pos, prob in sorted(spec.methylation.items()):
        q = ref_pos - spec.start   # simulate_read keeps M-only reads collinear
        if not 0 <= q < len(seq) or seq[q] != 'C':
            continue
        pairs.append((q, prob))
    if not pairs:
        return
    c_positions = [i for i, b in enumerate(seq) if b == 'C']
    deltas, probs, prev_rank = [], [], -1
    for q, prob in pairs:
        rank = c_positions.index(q)
        deltas.append(rank - prev_rank - 1)
        probs.append(min(255, int(round(prob * 255))))
        prev_rank = rank
    aln.set_tag('MM', 'C+m?,' + ','.join(str(d) for d in deltas) + ';')
    aln.set_tag('ML', probs)


def cpg_positions(ref_seq, lo, hi, count):
    """Reference positions of C bases in ``[lo, hi)``, evenly sampled."""
    cs = [i for i in range(lo, min(hi, len(ref_seq))) if ref_seq[i] == 'C']
    if len(cs) <= count:
        return cs
    step = len(cs) // count
    return cs[::step][:count]


# --------------------------------------------------------------------------- #
# a whole locus, wired together
# --------------------------------------------------------------------------- #
@dataclass
class Locus:
    ref_path: str
    ref_seq: str
    vcf_path: str
    bam_path: str
    variants: list
    chrom: str = CHROM

    @property
    def truth(self):
        """Per-variant (hap0 allele index, hap1 allele index)."""
        return [v.hap for v in self.variants]


class _NoBar:
    """Stand-in for tqdm: iterable when wrapping one, inert when used as a bar."""

    def __init__(self, iterable=None, **kwargs):
        self._iterable = iterable

    def __iter__(self):
        return iter(self._iterable if self._iterable is not None else ())

    def update(self, *args, **kwargs):
        pass

    def close(self):
        pass
