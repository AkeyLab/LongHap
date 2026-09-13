"""Fixtures for the longhap test suite.

Everything is synthetic and hermetic: a small reference contig, a bgzipped +
tabixed VCF, and a coordinate-sorted BAM, all built in ``tmp_path`` with pysam.
No external binaries beyond what pysam already bundles.

The builders these fixtures wrap live in ``synthetic.py`` -- see the note at the
top of that file for why they are not in this one.
"""
import pysam
import pytest

import longhap
from longhap import LongHap
from synthetic import (
    CHROM,
    Locus,
    ReadSpec,
    VCF_HEADER,
    _NoBar,
    _add_methylation,
    make_sequence,
    simulate_read,
)


# --------------------------------------------------------------------------- #
# reference
# --------------------------------------------------------------------------- #
@pytest.fixture
def write_reference(tmp_path):
    """Write and faidx-index a single-contig FASTA; returns (path, sequence)."""
    def _write(sequence, chrom=CHROM, name='ref.fa'):
        path = tmp_path / name
        with open(path, 'w') as fh:
            fh.write(f'>{chrom}\n')
            for i in range(0, len(sequence), 60):
                fh.write(sequence[i:i + 60] + '\n')
        pysam.faidx(str(path))
        return str(path), sequence
    return _write


# --------------------------------------------------------------------------- #
# variants
# --------------------------------------------------------------------------- #
@pytest.fixture
def write_vcf(tmp_path):
    """Write a bgzipped, tabixed VCF.

    ``records`` are ``Variant`` objects, or ``(variant, gt_string)`` pairs when a
    test needs a genotype longhap is meant to skip (``1/1``, ``./.``, ...).
    Extra header lines go in ``extra_header``; per-record INFO in ``infos``.
    """
    def _write(records, sequence_length=10000, chrom=CHROM, sample='S1',
               name='variants.vcf', extra_header='', infos=None):
        raw = tmp_path / name
        header = VCF_HEADER.format(chrom=chrom, length=sequence_length, sample=sample)
        if extra_header:
            header = header.replace('#CHROM', extra_header + '#CHROM')
        lines = [header]
        for n, rec in enumerate(records):
            var, gt = rec if isinstance(rec, tuple) else (rec, '0/1')
            info = (infos or {}).get(n, '.')
            lines.append('\t'.join([chrom, str(var.pos), var.vid, var.ref,
                                    ','.join(var.alt), '.', 'PASS', info,
                                    'GT', gt]) + '\n')
        raw.write_text(''.join(lines))
        gz = str(raw) + '.gz'
        pysam.tabix_compress(str(raw), gz, force=True)
        pysam.tabix_index(gz, preset='vcf', force=True)
        return gz
    return _write


# --------------------------------------------------------------------------- #
# reads
# --------------------------------------------------------------------------- #
@pytest.fixture
def write_bam(tmp_path):
    """Build a coordinate-sorted, indexed BAM from ``ReadSpec``s."""
    def _write(ref_seq, specs, variants=(), chrom=CHROM, name='reads.bam'):
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': chrom, 'LN': len(ref_seq)}]}
        unsorted = tmp_path / ('unsorted.' + name)
        with pysam.AlignmentFile(str(unsorted), 'wb', header=header) as out:
            for n, spec in enumerate(specs):
                seq, cigar = simulate_read(ref_seq, spec.start, spec.end,
                                           variants, spec.hap)
                a = pysam.AlignedSegment()
                a.query_name = spec.name or f'read{n}'
                a.query_sequence = seq
                a.reference_id = 0
                a.reference_start = spec.start
                a.mapping_quality = spec.mapq
                a.cigartuples = cigar
                a.query_qualities = pysam.qualitystring_to_array(
                    chr(33 + spec.base_quality) * len(seq))
                flag = spec.flags_extra
                if spec.is_reverse:
                    flag |= 16
                if spec.is_secondary:
                    flag |= 256
                if spec.is_supplementary:
                    flag |= 2048
                a.flag = flag
                if spec.methylation:
                    _add_methylation(a, seq, spec, ref_seq)
                out.write(a)
        path = tmp_path / name
        pysam.sort('-o', str(path), str(unsorted))
        pysam.index(str(path))
        return str(path)
    return _write


# --------------------------------------------------------------------------- #
# methylation
# --------------------------------------------------------------------------- #
@pytest.fixture
def write_methylation_bed(tmp_path):
    """Write a pb-CpG-tools style pileup BED.

    longhap reads this with ``skiprows=7``, so seven lines are burnt at the top
    exactly as the real tool emits them.  ``sites`` is ``[(start, coverage,
    ratio)]``; ``start`` is the 0-based reference position of the C.
    """
    def _write(sites, chrom=CHROM, name='methylation.bed', header_lines=7):
        path = tmp_path / name
        rows = ['#pb-CpG-tools placeholder header line\n'] * header_lines
        for start, coverage, ratio in sites:
            mod = int(round(coverage * ratio / 100))
            rows.append('\t'.join(str(f) for f in (
                chrom, start, start + 1, 1.0, 'Total', coverage, mod,
                coverage - mod, ratio)) + '\n')
        path.write_text(''.join(rows))
        return str(path)
    return _write


# --------------------------------------------------------------------------- #
# a whole locus, wired together
# --------------------------------------------------------------------------- #
@pytest.fixture
def make_locus(write_reference, write_vcf, write_bam):
    """Reference + VCF + BAM in one call.

    ``read_specs`` defaults to a tiling of 400 bp reads at 100 bp spacing,
    alternating haplotype, which is enough coverage to phase a locus whose
    variants are spaced under 200 bp apart.
    """
    def _make(variants, length=4000, read_specs=None, seed=0, genotypes=None,
              extra_header='', infos=None, read_length=400, step=100):
        ref_path, ref_seq = write_reference(make_sequence(length, seed))
        # ``variants`` may be a callable so a test can build REF alleles that
        # actually match the reference it was just handed.
        if callable(variants):
            variants = variants(ref_seq)
        records = variants if genotypes is None else list(zip(variants, genotypes))
        vcf_path = write_vcf(records, sequence_length=length,
                             extra_header=extra_header, infos=infos)
        if read_specs is None:
            read_specs = []
            for i, start in enumerate(range(0, length - read_length, step)):
                read_specs.append(ReadSpec(start, start + read_length, i % 2))
        bam_path = write_bam(ref_seq, read_specs, variants)
        return Locus(ref_path, ref_seq, vcf_path, bam_path, list(variants))
    return _make


@pytest.fixture
def run_longhap(tmp_path):
    """Construct a LongHap over a ``Locus`` and run the full pipeline.

    Returns the instance so a test can inspect internals as well as outputs.
    Pass ``phase=False`` to stop after transition inference, ``write=False`` to
    stop before the output files are produced.
    """
    def _run(locus, phase=True, write=True, **kwargs):
        opts = dict(vcf_f=locus.vcf_path, bam=locus.bam_path, chrom=locus.chrom,
                    reference_path=locus.ref_path,
                    output_vcf=str(tmp_path / 'phased.vcf.gz'),
                    output_blocks=str(tmp_path / 'blocks.bed'))
        opts.update(kwargs)
        lh = LongHap(**opts)
        lh.infer_variant_transitions()
        lh.infer_methylation_transitions()
        if phase:
            lh.phase()
            if write:
                lh.write_results()
        return lh
    return _run


@pytest.fixture
def read_phased_vcf():
    """Read a phased VCF back as ``[(pos, gt_string, ps)]``."""
    def _read(path):
        from cyvcf2 import VCF
        out = []
        for v in VCF(path):
            gt = v.genotypes[0]
            sep = '|' if gt[2] else '/'
            try:
                ps = v.format('PS')[0][0]
            except (KeyError, TypeError):
                ps = None
            out.append((v.POS, f'{gt[0]}{sep}{gt[1]}', ps))
        return out
    return _read


@pytest.fixture(autouse=True)
def _quiet_progress(monkeypatch):
    """Keep longhap's progress bars out of the test output."""
    monkeypatch.setattr(longhap, 'tqdm', _NoBar)
