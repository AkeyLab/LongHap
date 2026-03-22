#!/usr/bin/env python3
__version__ = "0.1.3"
import argparse
import os
import sys
from cyvcf2 import VCF, Writer
import pysam
from collections import deque, defaultdict
from tqdm import tqdm
import numpy as np
from pyfaidx import Fasta
import logging
import pandas as pd
from scipy.sparse import csc_array
import json
import re
from scipy.special import logsumexp
import parasail
logging.getLogger(__name__)


class NestedDict:
    def __init__(self):
        pass

    @classmethod
    def nested_defaultdict(cls):
        return defaultdict(list)


class LongHap:
    def __init__(self, vcf_f, bam, chrom, reference_path, output_vcf, output_blocks, output_bam=None,
                 output_read_assignments=None, methylation_calls_f=None,
                 output_transition_matrix=None, output_read_states=None,
                 output_variant_read_mapping=None, output_allele_coverage=None, output_transition_matrix_meth=None,
                 output_differentially_methylated_sites=None,
                 output_unphaseable_variants=None, snvs_only=False, multiallelics=False,
                 use_all_methylated_sites=False, max_meth_distance=5000, error_rate=1e-3, llr_thresh=4,
                 sample=None, force=False, max_allele_length=50000, min_allele_count=2, min_base_quality=0, min_mapq=20,
                 flank_snv=33, flank_indel=100, seqtech='pacbio'):
        self.chrom = chrom
        self.snvs_only = snvs_only
        self.multiallelics = multiallelics
        self.vcf_f = vcf_f
        self.bam = bam
        self.reference_path = reference_path
        self.methylation_calls_f = methylation_calls_f
        self.prev_methylations = None
        self.use_all_methylated_sites = use_all_methylated_sites
        self.max_meth_distance = max_meth_distance
        self.error_rate = error_rate
        self.llr_thresh = llr_thresh
        self.output_vcf = output_vcf
        self.output_blocks = output_blocks
        self.output_bam = output_bam
        self.output_read_assignments = output_read_assignments
        self.output_transition_matrix = output_transition_matrix
        self.output_allele_coverage = output_allele_coverage
        self.output_read_states = output_read_states
        self.output_variant_read_mapping = output_variant_read_mapping
        self.output_transition_matrix_meth = output_transition_matrix_meth
        self.output_differentially_methylated_sites = output_differentially_methylated_sites
        self.output_unphaseable_variants = output_unphaseable_variants
        self.sample = sample
        self.force = force
        self.max_allele_length = max_allele_length
        self.min_allele_count = min_allele_count
        self.min_base_quality = min_base_quality
        self.min_mapq = min_mapq
        self.flank_snv = flank_snv
        self.flank_indel = flank_indel
        self.seqtech = seqtech

        self.methylation_read_assignments = defaultdict(list)

        if self.chrom is None:
            variant_calls = VCF(self.vcf_f)
            self.chrom = variant_calls.seqnames[0]
            variant_calls.close()
        if self.sample is None:
            variant_calls = VCF(self.vcf_f)
            self.sample = variant_calls.samples[0]
            variant_calls.close()
        # TODO HOW DO ITERATE OVER MULTIPLE CHROM
        logging.info(f'Loading reference genome from {self.reference_path}')
        self.reference = self.read_reference_fasta(self.reference_path)
        self.reference = self.reference[self.chrom]
        logging.info(f'Loading variant calls from {self.vcf_f}')
        if not os.path.isfile(self.vcf_f):
            logging.error(f"VCF file {self.vcf_f} does not exist.")
            sys.exit(1)
        self.idx_variant_mapping, self.variant_idx_mapping, self.variant_type = self.get_heterozygous_variants()
        self.num_variants = len(self.idx_variant_mapping)
        # intialize transition_matrix
        self.transition_matrix = np.zeros((2, 2, self.num_variants - 1)) + 1e-20
        self.allele_coverage = np.zeros((2, self.num_variants))
        logging.info(f'Loading read alignments from {self.bam}')
        if not os.path.isfile(self.bam):
            logging.error(f"BAM file {self.bam} does not exist.")
            sys.exit(1)
        self.alignments = pysam.AlignmentFile(self.bam, 'rb')
        self.prev_state = dict()
        self.read_states = defaultdict(dict)
        self.variant_read_mapping = defaultdict(list)

        self.methylation_calls = pd.DataFrame(columns=['chrom', 'start', 'end', 'score', 'hap', 'coverage', 'mod_count',
                                                       'unmod_count', 'ratio'])
        self.differentially_methylated_sites = []

        self.haplotypes = np.zeros((2, self.transition_matrix.shape[2] + 1))
        self.delta = np.zeros((2, self.transition_matrix.shape[2] + 1))
        self.block_ends = []
        self.unphaseable = np.array([])
        self.phaseable = np.arange(self.num_variants)

    def infer_variant_transitions(self):
        """
        Infer transition matrix from read alignments
        """
        # check if transition matrix and read states already exist
        if (self.output_transition_matrix is not None and
                os.path.isfile(self.output_transition_matrix) and
                self.output_read_states is not None and
                os.path.isfile(self.output_read_states) and
                self.output_variant_read_mapping is not None and
                os.path.isfile(self.output_variant_read_mapping) and
                self.output_unphaseable_variants is not None and
                os.path.isfile(self.output_unphaseable_variants) and
                self.output_allele_coverage is not None and
                os.path.isfile(self.output_allele_coverage) and
                not self.force):
            logging.info(f'Loading allele coverage from {self.output_allele_coverage}')
            self.allele_coverage = np.load(self.output_allele_coverage)['arr_0']
            logging.info(f'Loading transition matrix from {self.output_transition_matrix}')
            self.transition_matrix = np.load(self.output_transition_matrix)['arr_0']
            logging.info(f'Loading read_states from {self.output_read_states}')
            self.read_states = defaultdict(dict, json.load(open(self.output_read_states, 'r')))
            self.variant_read_mapping = defaultdict(list, json.load(open(self.output_variant_read_mapping, 'r')))
            self.unphaseable = np.load(self.output_unphaseable_variants)['arr_0']
            self.phaseable = self.phaseable[~np.isin(self.phaseable, self.unphaseable)]
        else:
            logging.info('Inferring transition matrix from variant data')
            self.create_directed_graph_of_heterozygous_variants_from_reads()

            logging.info('Rephasing complex variants and variants with low MAC')
            self.rephase_difficult_variants()

            for i in range(self.transition_matrix.shape[2]):
                self.transition_matrix[:, :, i] = self.mirror_transition(self.transition_matrix[:, :, i],
                                                                         normalized=False)
            self.transition_matrix /= self.transition_matrix.sum(axis=1, keepdims=True)
            self.connect_phase_blocks()

            if self.output_allele_coverage is not None:
                np.savez(self.output_allele_coverage, self.allele_coverage)
            if self.output_transition_matrix is not None:
                np.savez(self.output_transition_matrix, self.transition_matrix)
            if self.output_read_states is not None:
                # np.save(output_read_states, read_states)
                json.dump(self.read_states, open(self.output_read_states, 'w'))
            if self.output_variant_read_mapping is not None:
                json.dump(self.variant_read_mapping, open(self.output_variant_read_mapping, 'w'))
            if self.output_unphaseable_variants is not None:
                np.savez(self.output_unphaseable_variants, self.unphaseable)

    def infer_methylation_transitions(self):
        """
        Infer transition matrix from methylation calls
        """
        # check if necessary information for methylation phasing are given and if intermediate files can be re-used
        if (self.methylation_calls_f is not None and
            (self.output_transition_matrix_meth is None or not os.path.isfile(self.output_transition_matrix_meth) or
                self.force)):
            logging.info(f'Loading methylation calls from {self.methylation_calls_f}')
            if not os.path.isfile(self.methylation_calls_f):
                logging.error(f"Methylation calls file {self.methylation_calls_f} does not exist.")
                sys.exit(1)
            self.methylation_calls = pd.read_csv(self.methylation_calls_f, sep='\t',
                                                 names=['chrom', 'start', 'end', 'score', 'hap',
                                                        'coverage', 'mod_count', 'unmod_count',
                                                        'ratio'], engine='pyarrow', skiprows=7)
            # get putative differentially methylated sites
            self.methylation_calls = self.methylation_calls[(self.methylation_calls.chrom == self.chrom) &
                                                            (self.methylation_calls.coverage >= 10) &
                                                            (self.methylation_calls.ratio > 20) &
                                                            (self.methylation_calls.ratio < 80)]
            logging.info('Complementing variant transition matrix with methylation data')
            # fill in transition matrix
            self.get_methylation_transitions_helper()
            if len(self.differentially_methylated_sites) > 0:
                self.differentially_methylated_sites = pd.concat(
                    self.differentially_methylated_sites).sort_values(['chrom', 'start', 'hap']).drop_duplicates()
            else:
                self.differentially_methylated_sites = pd.DataFrame()
        elif self.methylation_calls_f is not None and os.path.isfile(self.output_transition_matrix_meth):
            logging.info(
                f'Loading methylation complemented transition matrix from {self.output_transition_matrix_meth}')
            self.transition_matrix = np.load(self.output_transition_matrix_meth)['arr_0']

    def phase(self):
        """
        Phase variants using Viterbi algorithm
        """
        logging.info("Performing backtracing")
        # do backtracing
        self.calculate_forward_path_probabilities()

    def write_results(self):
        """
        Write phasing results to output files
        """
        # write phased vcf
        if self.output_vcf:
            logging.info(f"Writing phased VCF to {self.output_vcf}")
            self.write_phased_vcf()
        if self.output_blocks:
            # write haplotype blocks
            logging.info(f"Writing phase block coordinates to {self.output_blocks}")
            self.write_haplotype_blocks()
        # write haplotagged bam
        if self.output_bam:
            logging.info(f"Writing haplotagged bam to {self.output_bam}")
            self.haplotag_reads()

        if self.output_differentially_methylated_sites and len(self.differentially_methylated_sites) > 0:
            logging.info(f"Writing used differentially methylated sites to "
                         f"{self.output_differentially_methylated_sites}")
            self.differentially_methylated_sites.to_csv(self.output_differentially_methylated_sites, sep='\t',
                                                        index=False)
        elif self.output_differentially_methylated_sites and len(self.differentially_methylated_sites) == 0:
            logging.info(f"Did not identify any differentially methylated sites. Just creating "
                         f"{self.output_differentially_methylated_sites}")
            f = open(self.output_differentially_methylated_sites, 'w')
            f.write('chrom\tstart\tend\tscore\thap\tcoverage\tmod_count\tunmod_count\tratio\n')
            f.close()

    @staticmethod
    def read_reference_fasta(filepath):
        """
        Parse reference fasta file
        :param filepath: str, file to path
        :return: dict, reference sequence
        """
        if not os.path.isfile(filepath):
            logging.error(f"Reference fasta file {filepath} does not exist.")
            sys.exit(1)
        if not os.path.isfile(f'{filepath}.fai'):
            logging.error(f"Index for reference fasta file {filepath} does not exist. "
                          f"Index with samtools faidx {filepath}")
            sys.exit(1)
        fasta = Fasta(filepath)
        return fasta

    def get_heterozygous_variants(self):
        """
        Parse heterozygous variants from VCF
        :return: dict, for each heterozygous variant index it records "POS", "REF", and "ALT", "gt" fields
        """
        variant_type = []
        variant_calls = VCF(self.vcf_f)

        # indicates VCF was generated by vg/pggb and contains node traversal info for each allele
        has_allele_traversal_info = variant_calls.contains('AT')
        has_allele_level_info = variant_calls.contains('LV')  # has bubble info level

        idx_variant_mapping = dict()
        variant_idx_mapping = dict()
        i = 0
        # only consider heterozygous variants
        for variant in variant_calls(self.chrom):
            if (variant.gt_types == 1 and
                    -1 not in variant.genotypes[0][:2] and
                    (variant.is_snp or (variant.is_indel and not self.snvs_only)) and
                    (len(variant.ALT) == 1 or self.multiallelics) and  # only include multiallelics if specified
                    np.max([len(variant.REF), max([len(a) for a in variant.ALT])]) <= self.max_allele_length):
                # skip complex variant that was split into 2
                if has_allele_level_info and "_" in variant.ID:
                    continue
                idx_variant_mapping[i] = {"POS": variant.POS, 'REF': variant.REF.upper(),
                                          'ALT': [a.upper() for a in variant.ALT],
                                          'alleles': [variant.REF.upper()] + [a.upper() for a in variant.ALT],
                                          'gt': sorted(variant.genotypes[0][:2]),
                                          'ID': variant.ID}
                variant_idx_mapping[f'{variant.POS}_{variant.ID}'] = {"idx": i, 'REF': variant.REF.upper(),
                                                                      'ALT': [a.upper() for a in variant.ALT],
                                                                      'alleles': [variant.REF.upper()] +
                                                                                 [a.upper() for a in variant.ALT],
                                                                      'gt': sorted(variant.genotypes[0][:2]),
                                                                      'POS': variant.POS}
                if has_allele_traversal_info:
                    idx_variant_mapping[i]['AT'] = [re.split(">|<", at)[1:]
                                                    for at in variant.INFO['AT'].split(',')]
                    variant_idx_mapping[f'{variant.POS}_{variant.ID}']['AT'] = [re.split(">|<", at)[1:]
                                                                                for at in variant.INFO['AT'].split(',')]
                if len(idx_variant_mapping[i]["ALT"]) > 1:
                    variant_type.append('MULTI')
                elif len(idx_variant_mapping[i]["REF"]) == len(idx_variant_mapping[i]["ALT"][0]) == 1:
                    variant_type.append('SNP')
                elif len(idx_variant_mapping[i]["REF"]) > len(idx_variant_mapping[i]["ALT"][0]):
                    variant_type.append('DEL')
                elif len(idx_variant_mapping[i]["REF"]) < len(idx_variant_mapping[i]["ALT"][0]):
                    variant_type.append('INS')
                else:
                    variant_type.append("Other")

                i += 1

        variant_type = np.array(variant_type)
        return idx_variant_mapping, variant_idx_mapping, variant_type

    @staticmethod
    def get_query_and_reference_indices_for_variant_cigar(cigar, position, read_start, ref_offset, query_offset,
                                                          operation, length):
        """
        This function returns the query and reference index of matching pairs in an alignment at a specified position
        leveraging the CIGAR string.
        :param cigar: list, cigar tuples returned by pysam.AlignmentSegment.cigartuples
        :param position: int, target position in reference coordinates
        :param read_start: int, read start in reference coordinates
        :param ref_offset: int, reference sequence index offset (i.e., max index that has already been parsed)
        :param query_offset: int, query sequence index offset (i.e., max index that has already been parsed)
        :param operation: int, last CIGAR operation
        :param length: int, remaining length of last cigar operation
        :return: tuple, query and reference index for target position and last CIGAR operation and remaining length
        """
        offset = position - read_start
        while cigar:
            # do not count soft clips and insertions
            if length + ref_offset > offset and operation != 4 and operation != 1:
                break
            # match
            if operation == 0 or operation == 7:
                ref_offset += length
                query_offset += length
            # mismatch
            elif operation == 8:
                ref_offset += length
                query_offset += length
            # insertion, gap in reference
            elif operation == 1:
                query_offset += length
            # deletion, gap in query
            elif operation == 2:
                ref_offset += length
            # soft clip, query does not appear in alignment
            elif operation == 4:
                query_offset += length
            operation, length = cigar.popleft()
        # last operation is soft clip --> don't increment reference offset
        # TODO make sure that this works as intended
        if operation != 2:
            query_offset += offset - ref_offset
        # last operation is a match, mismatch, or deletion
        if operation != 1 and operation != 4:
            length -= offset - ref_offset
            ref_offset = offset
        # last operation is an insertion or soft clip --> do not increment reference offset
        else:
            query_offset += length
            length = 0
        return query_offset, ref_offset, operation, length

    def get_state_at_variant(self, query_sequence, base_qualities, cigar, read_start, ref_offset, query_offset,
                             operation, length, variant):
        """
        Get allelic state of read at variant. If state cannot be determined return None so that read-variant pair
        is considered for realignment
        :param query_sequence: str, read query sequence
        :param base_qualities: list-like, read base qualities
        :param cigar: list-like, subsequent cigar operation
        :param read_start: int, 0-index start coordinate of read
        :param ref_offset: int, number of bases already parsed in reference from read_start
        :param query_offset: int, 0-index number of bases already parsed in read
        :param operation: int, current CIGAR operation
        :param length: int, length of current CIGAR operation
        :param variant: dict, focal variant
        :return: tuple, (state, query_offset, ref_offset, operation, length)
        """
        allele_ref = variant['REF']
        allele_a = variant['alleles'][variant['gt'][0]]
        allele_b = variant['alleles'][variant['gt'][1]]
        position = variant['POS'] - 1

        qpos, r_idx, operation, length = self.get_query_and_reference_indices_for_variant_cigar(cigar, position,
                                                                                                read_start, ref_offset,
                                                                                                query_offset, operation,
                                                                                                length)
        if qpos == -1:
            return -1

        state = None
        # SNV
        if len(allele_ref) == len(allele_a) == len(allele_b) == 1:
            if self.seqtech == "pacbio":
                allele_q = query_sequence[qpos]
                base_quality = base_qualities[qpos]
                if base_quality <= self.min_base_quality:
                    state = -1
                elif allele_q == allele_a:
                    state = 0
                elif allele_q == allele_b:
                    state = 1
                else:
                    # do realignment
                    pass
            # There is a small benefit to this but it's slow because I do more realignments. I need to do something else
            elif self.seqtech == "ont":
                query_length = len(query_sequence)
                upstream_base = self.reference[r_idx + read_start - 1].seq.upper()
                downstream_base = self.reference[r_idx + read_start + 1].seq.upper()

                allele_a = upstream_base + allele_a + downstream_base
                allele_b = upstream_base + allele_b + downstream_base
                allele_q = query_sequence[np.max([0, qpos - 1]): np.min([query_length, qpos + 2])]
                base_quality = np.min(base_qualities[np.max([0, qpos - 1]): np.min([query_length, qpos + 2])])
                if base_quality <= self.min_base_quality:
                    # do realignment
                    state = -1
                elif allele_q == allele_a:
                    state = 0
                elif allele_q == allele_b:
                    state = 1
                else:
                    # do realignment
                    pass

        # insertion
        elif allele_ref == allele_a and len(allele_ref) < len(allele_b):
            cigar_next = cigar.copy()
            q_after, r_next, operation_next, length_next = (
                self.get_query_and_reference_indices_for_variant_cigar(cigar_next, position + 1, read_start,
                                                                       r_idx, qpos, operation, length))

            if qpos == -1 or q_after == -1 or q_after <= qpos:
                state = -1
            insertion_seq = query_sequence[qpos + 1:q_after]
            if insertion_seq == "":
                state = 0
            elif insertion_seq == allele_b[len(allele_ref):]:
                state = 1
            else:
                # do realignment
                pass

        # deletion
        elif allele_ref == allele_a and len(allele_ref) > len(allele_b):
            cigar_next = cigar.copy()
            q_after, r_next, operation_next, length_next = (
                self.get_query_and_reference_indices_for_variant_cigar(cigar_next, position + len(allele_ref), read_start,
                                                                       r_idx, qpos, operation, length))
            deleted = np.max([0, (r_next - r_idx) - (q_after - qpos)])
            # nothing is deleted --> matches ref
            if deleted == 0:
                state = 0
            # deletion matches length different of ref and allele b --> allele b
            elif deleted == len(allele_ref) - len(allele_b):
                state = 1
            else:
                # do realignment
                pass
        # complex variant --> do a realignment later
        else:
            pass

        return state, qpos, r_idx, operation, length

    @staticmethod
    def get_adaptive_gap_penalties(cigartuples):
        matches = np.sum([length for op, length in cigartuples if op in (0, 7)])
        mismatches = np.sum([length for op, length in cigartuples if op == 8])
        ins = np.sum([length for op, length in cigartuples if op == 1])
        dels = np.sum([length for op, length in cigartuples if op == 2])
        total = matches + mismatches + ins + dels
        if total == 0:
            return 5, 1
        indel_rate = (ins + dels) / total
        mismatch_rate = mismatches / total
        if indel_rate < 0.005:
            gap_open, gap_extend = 7, 2
        elif indel_rate < 0.02:
            gap_open, gap_extend = 5, 1
        elif indel_rate < 0.05:
            gap_open, gap_extend = 4, 1
        else:
            gap_open, gap_extend = 3, 1
        return gap_open, gap_extend, indel_rate, mismatch_rate

    def realign_around_variant(self, query_sequence, variant, qpos,
                               gap_open, gap_extend, homopolymer=False):
        """
        Realign sub-reads around complex variants
        :param query_sequence:  str, read query sequence
        :param variant: dict, complex variant
        :param qpos: int, index of variant in read
        :param gap_open: int, gap open penalty for alignment
        :param gap_extend: int, gap extend penalty for alignment
        :param homopolymer: boolean, whether the variant is in a homopolymer region
        :return: int, allelic state
        """
        # extract variant info
        allele_ref = variant['REF']
        allele_a = variant['alleles'][variant['gt'][0]]
        allele_b = variant['alleles'][variant['gt'][1]]
        position = variant['POS'] - 1
        max_allele_len = np.max([len(allele_ref), len(allele_a), len(allele_b)])
        if max_allele_len == 1:
            flank = self.flank_snv
        else:
            flank = self.flank_indel
        # for SNPs choose smaller window and penalize gaps more
        if max_allele_len == 1:
            gap_open = int(gap_open * 2)
        if homopolymer:
            gap_open -= 2
        # get read window
        start = np.max([0, qpos - flank])
        end = np.min([len(query_sequence),
                      qpos + max_allele_len + flank])

        read_window = query_sequence[start:end]
        # get corresponding reference window with allele A and B inserted respectively
        flank_up = qpos - start
        flank_down = end - (qpos + np.max([len(allele_ref), len(allele_a), len(allele_b)]))
        upstream = self.reference[position - flank_up: position].seq.upper()

        downstream_a = self.reference[position + len(allele_ref): position + len(allele_ref) + flank_down +
                                                                  (max_allele_len - len(allele_a))].seq.upper()

        ref_allele_a = str(upstream + allele_a + downstream_a)

        downstream_b = self.reference[position + len(allele_ref): position + len(allele_ref) + flank_down +
                                                                  (max_allele_len - len(allele_b))].seq.upper()

        ref_allele_b = str(upstream + allele_b + downstream_b)

        aln_a = parasail.sg_stats_striped_sat(read_window, ref_allele_a, gap_open, gap_extend, parasail.dnafull)
        aln_b = parasail.sg_stats_striped_sat(read_window, ref_allele_b, gap_open, gap_extend, parasail.dnafull)

        # Compare alignment scores
        score_a = aln_a.score
        score_b = aln_b.score

        if score_a > score_b:
            state = 0
        elif score_b > score_a:
            state = 1
        else:
            state = -1

        return state

    @staticmethod
    def mirror_transition(t, normalized=True):
        """
        Mirror certain transition to resolve uncertain transition and create symmetric transition matrix
        :param t: 2x2 np.array, transition matrix
        :param normalized: boolean, whether the transition_matrix has been normalized to values between 0 and 1
        :return: 2x2 np.array, mirrored, symmetric transition matrix
        """
        # # step from alternative allele is more certain than from reference allele
        if (np.abs(np.log(t[0, 0] / t[0, 1])) < np.abs(np.log(t[1, 0] / t[1, 1])) and
                (t[1, :].max() > 1 or t[0, :].max() <= 1 or normalized) and np.unique(t.argmax(axis=1)).shape[0] > 1):
            # mirror transitions of alternative allele for reference allele
            t[0, :] = t[1, :][::-1]
        # step from reference allele is more certain than from alternative allele
        elif (np.abs(np.log(t[0, 0] / t[0, 1])) > np.abs(np.log(t[1, 0] / t[1, 1])) and
              (t[0, :].max() > 1 or t[1, :].max() <= 1 or normalized) and np.unique(t.argmax(axis=1)).shape[0] > 1):
            # mirror transitions of reference allele for alternative allele
            t[1, :] = t[0, :][::-1]
        # mirror transitions of reference allele for alternative allele
        # as it has the higher confidence despite lower certainty
        elif (np.abs(np.log(t[0, 0] / t[0, 1])) < np.abs(np.log(t[1, 0] / t[1, 1])) and
                (t[1, :].max() == 1 and t[0, :].max() > 1 and not normalized)):
            t[1, :] = t[0, :][::-1]
        # step from alternative allele has higher confidence than from reference allele
        elif (np.abs(np.log(t[0, 0] / t[0, 1])) > np.abs(np.log(t[1, 0] / t[1, 1])) and
              (t[0, :].max() == 1 and t[1, :].max() > 1 and not normalized)):
            t[0, :] = t[1, :][::-1]
        elif np.all(t[0, :] == t[1, :]) and not normalized:
            t[:, :] = 1e-20
        elif np.all(t[0, :] == t[1, :]) and normalized:
            t[:, :] = 0.5
        elif np.unique(t.argmax(axis=1)).shape[0] == 1:
            diff = np.abs(t[:, 0] - t[:, 1])
            # cov = t.sum(axis=1)
            if diff[0] > diff[1]:
                t[1, :] = t[0, :][::-1]
            elif diff[0] < diff[1]:
                t[0, :] = t[1, :][::-1]
            elif np.abs(np.log(t[0, 0] / t[0, 1])) < np.abs(np.log(t[1, 0] / t[1, 1])):
                t[0, :] = t[1, :][::-1]
            elif np.abs(np.log(t[0, 0] / t[0, 1])) > np.abs(np.log(t[1, 0] / t[1, 1])):
                t[1, :] = t[0, :][::-1]
            else:
                if not normalized:
                    t[:, :] = 1e-20
                else:
                    t[:, :] = 0.5

        return t

    def get_allele_transitions_from_known_read_states(self, idx_var_a, idx_var_b):
        """
        Count occurrences of combinations of two adjacent heterozygous variants in overlapping reads from
        known read states
        :param idx_var_a: int, index of first heterozygous variant
        :param idx_var_b: int, index of second heterozygous variant
        """
        transition_matrix = np.zeros((2, 2)) + 1e-20
        reads_a = self.variant_read_mapping[str(idx_var_a)]
        reads_b = self.variant_read_mapping[str(idx_var_b)]
        reads = list(set(reads_a) & set(reads_b))
        for read in reads:
            state_a = self.read_states[read][str(idx_var_a)]
            state_b = self.read_states[read][str(idx_var_b)]
            if state_a != -1 and state_b != -1:
                transition_matrix[state_a, state_b] += 1
        return transition_matrix

    @staticmethod
    def get_methylation_status_at_site_cigar(pos, cigar, read_start, ref_offset, query_offset, operation, length,
                                             mm_tags, ml_tags):
        """
        Get methylation probability for a read at a given site
        :param pos: int, target position in reference coordinates
        :param cigar: list, CIGAR tuples
        :param read_start: int, reference coordinate of alignment start
        :param ref_offset: int, reference sequence index offset (i.e., max index that has already been parsed)
        :param query_offset: int, query sequence index offset (i.e., max index that has already been parsed)
        :param operation: int, last CIGAR operation
        :param length: int, remaining length of last cigar operation
        :param mm_tags: np.array, sites with modification information
        :param ml_tags: np.array, modification probabilities corresponding to sites in mm_tags
        :return: tuple, methylation probability, reference index, query index,
                        last CIGAR operation, remaining length of last CIGAR operation
        """
        q_idx, r_idx, operation, length = LongHap.get_query_and_reference_indices_for_variant_cigar(cigar, pos,
                                                                                                    read_start,
                                                                                                    ref_offset,
                                                                                                    query_offset,
                                                                                                    operation, length)
        if operation == 2 or operation == 1 or operation == 4 or r_idx + read_start != pos:
            meth = -1
        else:
            # get methylation probabilities for target sites
            indices = np.where(mm_tags == q_idx)[0]
            if indices.size > 0:
                meth = ml_tags[indices[0]]
            else:
                meth = -1
        return meth, r_idx, q_idx, operation, length

    @staticmethod
    def calculate_transition_probability_from_methylation_helper(hap, i):
        """
        Calculate probabilities of phase transitions between variants based on methylation data
        :param hap: np.array, haplotype inferred based on methylation data
        :param i: int, number of subsequent transitions that were inferred
        :return: np.array, modified phase transition matrix
        """
        transition_matrix = np.zeros((2, 2)) + 1e-20
        if np.round(hap).astype(int)[i] == 1 and np.round(hap).astype(int)[i + 1] == 1:
            transition_matrix[1, 1] = np.max([hap[i], 1e-20]) * np.max([hap[i + 1], 1e-20])
            transition_matrix[1, 0] = np.max([1 - np.max([hap[i], 1e-20]) * np.max([hap[i + 1], 1e-20]), 1e-20])
        elif np.round(hap).astype(int)[i] == 0 and np.round(hap).astype(int)[i + 1] == 0:
            transition_matrix[0, 0] = np.max([1 - np.max([hap[i], 1e-20]) * np.max([hap[i + 1], 1e-20]), 1e-20])
            transition_matrix[0, 1] = np.max([hap[i], 1e-20]) * np.max([hap[i + 1], 1e-20])
        elif np.round(hap).astype(int)[i] == 0 and np.round(hap).astype(int)[i + 1] == 1:
            transition_matrix[0, 0] = np.max([1 - np.max([1 - hap[i], 1e-20]) * np.max([hap[i + 1], 1e-20]), 1e-20])
            transition_matrix[0, 1] = np.max([1 - hap[i], 1e-20]) * np.max([hap[i + 1], 1e-20])
        elif np.round(hap).astype(int)[i] == 1 and np.round(hap).astype(int)[i + 1] == 0:
            transition_matrix[1, 1] = np.max([1 - np.max([hap[i], 1e-20]) * np.max([1 - hap[i + 1], 1e-20]), 1e-20])
            transition_matrix[1, 0] = np.max([hap[i], 1e-20]) * np.max([1 - hap[i + 1], 1e-20])

        transition_matrix = LongHap.mirror_transition(transition_matrix, normalized=True)
        return transition_matrix

    def calculate_transition_probability_from_methylation(self, hap1, hap2, hap1_cov, hap2_cov, i, idx):
        """
        Calculate probabilities of phase transitions between variants based on methylation data
        as average of both phases.
        :param hap1: np.array, haplotype A inferred based on methylation data
        :param hap2: np.array, haplotype B inferred based on methylation data
        :param hap1_cov: np.array, coverage of haplotype A inferred based on methylation data
        :param hap2_cov: np.array, coverage of haplotype B inferred based on methylation data
        :param i: int, index of target transition in current hap that was inferred from methylation dats
        :param idx: int, index of target transition that was inferred from methylation dats
        :return: np.array, modified phase transition matrix
        """
        if (hap1_cov[i] < self.min_allele_count or hap1_cov[i + 1] < self.min_allele_count or
                hap2_cov[i] < self.min_allele_count or hap2_cov[i + 1] < self.min_allele_count):
            pass
        else:
            t1 = self.calculate_transition_probability_from_methylation_helper(hap1, i)
            t2 = self.calculate_transition_probability_from_methylation_helper(hap2, i)
            cov_hap1 = (hap1_cov[i] + hap1_cov[i + 1]) / 2
            cov_hap2 = (hap2_cov[i] + hap2_cov[i + 1]) / 2
            self.transition_matrix[:, :, idx] = (t1 * cov_hap1 + t2 * cov_hap2) / (cov_hap1 + cov_hap2 + 1e-100)

    @staticmethod
    def get_read_methylation(methylations, positions, cigar, read_start, read_end, ref_offset, query_offset, operation,
                             length, mm_tags, ml_tags):
        """
        Get methylation probabilities for a given read
        :param methylations: dict, previously inferred methylation probabilities
        :param positions: list-like, positions for which to extract methylation probabilities
        :param cigar: list, CIGAR tuples
        :param read_start: int, reference coordinate of alignment start
        :param read_end: int, reference coordinate of alignment end
        :param ref_offset: int, reference sequence index offset (i.e., max index that has already been parsed)
        :param query_offset: int, query sequence index offset (i.e., max index that has already been parsed)
        :param operation: int, last CIGAR operation
        :param length: int, remaining length of last cigar operation
        :param mm_tags: np.array, sites with modification information
        :param ml_tags: np.array, modification probabilities corresponding to sites in mm_tags
        :return:
        """
        read_methylation = []
        for m_pos in positions:
            if m_pos in methylations:
                meth = methylations[m_pos]
            elif read_start <= m_pos <= read_end:
                (meth, ref_offset, query_offset, operation,
                 length) = LongHap.get_methylation_status_at_site_cigar(m_pos, cigar, read_start, ref_offset,
                                                                        query_offset, operation, length, mm_tags,
                                                                        ml_tags)
            else:
                meth = -1
            read_methylation.append(meth)
        read_methylation = np.array(read_methylation)
        return read_methylation, cigar, ref_offset, query_offset, operation, length

    @staticmethod
    def get_read_info(read):
        """
        Get read meta information
        :param read: pysam.AlignedSegment, current alignement
        :return: tupe, (CIGAR tuples, read start, read end, reference offset, query offset, last cigar operation,
                        last cigar operation length, methylation tags, methylation likelihood)
        """
        # pos_mapping_meth = deque(read.get_aligned_pairs())
        cigar = deque(read.cigartuples)
        read_start = read.reference_start
        read_end = read.reference_end
        ref_offset = 0
        query_offset = 0
        operation = None
        length = 0
        # get methylated Cytosines
        if read.is_forward:
            modified_bases = read.modified_bases_forward
            key = ('C', 0, 'm')
        else:
            modified_bases = read.modified_bases
            key = ('C', 1, 'm')

        # TODO: why do I use forwared for forward and normal for reverse? Why do I need to subtract -1 for reverse?
        # I confirmed that this extraction scheme only yields C bases
        if len(modified_bases) > 0 and key in modified_bases:
            modified_bases = modified_bases[key]
            mm_tags, ml_tags = zip(*modified_bases)
            mm_tags = np.array(mm_tags) - (0 if read.is_forward else 1)
            ml_tags = np.array(ml_tags) / 256
        else:
            mm_tags = np.array([])
            ml_tags = np.array([])
        return cigar, read_start, read_end, ref_offset, query_offset, operation, length, mm_tags, ml_tags

    @staticmethod
    def calculate_probability_of_reads_belonging_to_haplotype_based_on_methylation(meth_states_hap, meth_hap,
                                                                                   unmeth_hap, methylation_probs,
                                                                                   unmethylated_probs, unassigned,
                                                                                   diff_meth):
        """
        Calculate the probability that a reads comes from a given haplotype based on methylation patterns
        :param meth_states_hap: np.array, methylation pattern of haplotypes
        :param meth_hap: np.array, methylated sites in haplotype
        :param unmeth_hap: np.array, unmethylated sites in haplotype
        :param methylation_probs: scipy.sparse.csc_array, methylation probabilities for each read (rows)
                                    at each site (cols)
        :param unmethylated_probs: scipy.sparse.csc_array, non-methylation probabilities for each read (rows)
                                    at each site (cols)
        :param unassigned: np.array, unassigned reads
        :param diff_meth: np.array, differentially methylated sites
        :return: np.array, probabilities that unassigned reads come from given haplotype
        """
        # get nonzero entries
        nonzero_meth_hap = methylation_probs[unassigned][:,
                           diff_meth[meth_states_hap[diff_meth] == 1]].nonzero()
        nonzero_unmeth_hap = unmethylated_probs[unassigned][:,
                             diff_meth[meth_states_hap[diff_meth] == 0]].nonzero()
        # calculate probability of read coming from haplotype based on haplotype methylation and read methylation
        # using methylated sites
        methylation_probs[unassigned[nonzero_meth_hap[0]],
                          diff_meth[meth_states_hap[diff_meth] == 1][nonzero_meth_hap[1]]] += \
            meth_hap[diff_meth][meth_states_hap[diff_meth] == 1][nonzero_meth_hap[1]]
        # calculate probability of read coming from haplotype based on haplotype methylation and read methylation
        # using unmethylated sites
        unmethylated_probs[unassigned[nonzero_unmeth_hap[0]],
                           diff_meth[meth_states_hap[diff_meth] == 0][nonzero_unmeth_hap[1]]] += \
            unmeth_hap[diff_meth][meth_states_hap[diff_meth] == 0][nonzero_unmeth_hap[1]]
        p_hap = (methylation_probs[unassigned][:, diff_meth[meth_states_hap[diff_meth] == 1]].sum(axis=1) +
                 unmethylated_probs[unassigned][:, diff_meth[meth_states_hap[diff_meth] == 0]].sum(axis=1))
        return p_hap

    @staticmethod
    def get_diff_methylation_sites_per_hap(meth_calls, meth_probs, reads, hap):
        """
        Get differentially methylated sites per haplotype
        :param meth_calls: pd.DataFrame, methylation calls
        :param meth_probs: scipy.sparse.csc_array, methylation probabilities for each read (rows)
                                    at each site (cols)
        :param reads: list, read IDs
        :param hap: int, haplotype
        :return: pd.DataFrame, methylation calls with differentially methylated sites per haplotype
        """
        df_hap = meth_calls.copy()
        c_meth_probs = meth_probs[reads]
        cov = np.zeros((c_meth_probs.shape[1]))
        mod = np.zeros((c_meth_probs.shape[1]))
        var_id, var_cov = np.unique(c_meth_probs.nonzero()[1], return_counts=True)
        cov[var_id] = var_cov
        mod_condition_mask = 10 ** c_meth_probs.data > 0.5
        for i in range(c_meth_probs.shape[1]):
            # Get the start and end indices in the data array for the current column
            start_idx = c_meth_probs.indptr[i]
            end_idx = c_meth_probs.indptr[i + 1]
            mod[i] = np.sum(mod_condition_mask[start_idx:end_idx])

        df_hap['hap'] = hap
        df_hap['coverage'] = cov
        df_hap['mod_count'] = mod
        df_hap['unmod_count'] = cov - mod
        df_hap['ratio'] = (mod / cov) * 100
        return df_hap

    def get_methylation_transitions(self, idx_var_a, idx_var_b):
        """
        Infer transition probabilities between variant A and B from methylation data
        :param idx_var_a: int, index of variant A
        :param idx_var_b: int, index of variant B
        :return: np.array and dict if prev_methylation is not None, filled in transition matrix
                    and inferred methylation states of reads and extracted meta information
        """
        # 1-indexed start and end positions that are within the reach of reads defined by max_read_length
        pos_a = self.idx_variant_mapping[self.phaseable[idx_var_a]]["POS"]
        pos_b = self.idx_variant_mapping[self.phaseable[idx_var_b]]["POS"]

        start = np.max([0, pos_a - self.max_meth_distance])
        end = np.min([pos_b + self.max_meth_distance, len(self.reference)])

        # get putatively differentially methylated sites in target region
        c_methylation_calls = self.methylation_calls[(self.methylation_calls.start >= start - 1) &
                                                     (self.methylation_calls.end <= end - 1)]
        if c_methylation_calls.shape[0] > 25000 and not self.use_all_methylated_sites:
            c_methylation_calls = c_methylation_calls[(c_methylation_calls.ratio > 45) &
                                                      (c_methylation_calls.ratio < 55) &
                                                      (c_methylation_calls.coverage > 15) &
                                                      (c_methylation_calls.coverage < 45)]
        if c_methylation_calls.shape[0] > 25000 and not self.use_all_methylated_sites:
            c_methylation_calls = c_methylation_calls.sample(n=25000)
        positions = c_methylation_calls.start.values
        methylation_probs = defaultdict(list)
        read_variant_states = []
        read_ids = []
        variant_homopolymer_mapping = {}
        for read in self.alignments.fetch(self.chrom, start, end + 1):
            # only consider primary alignments
            if (read.is_secondary or read.is_duplicate or read.is_unmapped or
                    read.is_qcfail or read.mapping_quality < self.min_mapq or not read.has_tag("MM")):
                continue
            read_name = read.query_name
            if self.prev_methylations is not None:
                if read_name in self.prev_methylations:
                    (methylations, cigar, read_start, read_end, ref_offset,
                     query_offset, operation, length, mm_tags, ml_tags) = self.prev_methylations[read_name]
                    (read_methylation, cigar, ref_offset,
                     query_offset, operation, length) = self.get_read_methylation(methylations, positions, cigar,
                                                                                  read_start, read_end, ref_offset,
                                                                                  query_offset, operation,
                                                                                  length, mm_tags, ml_tags)
                else:
                    (cigar, read_start, read_end, ref_offset, query_offset,
                     operation, length, mm_tags, ml_tags) = self.get_read_info(read)
                    (read_methylation, cigar, ref_offset,
                     query_offset, operation, length) = self.get_read_methylation({}, positions, cigar,
                                                                                  read_start, read_end, ref_offset,
                                                                                  query_offset, operation,
                                                                                  length, mm_tags, ml_tags)
                self.prev_methylations[read_name] = ({m_pos: m for m_pos, m in zip(positions, read_methylation)},
                                                           cigar, read_start, read_end, ref_offset, query_offset,
                                                           operation, length, mm_tags, ml_tags)
            else:
                (cigar, read_start, read_end, ref_offset, query_offset,
                 operation, length, mm_tags, ml_tags) = self.get_read_info(read)
                (read_methylation, cigar, ref_offset,
                 query_offset, operation, length) = self.get_read_methylation({}, positions, cigar,
                                                                              read_start, read_end, ref_offset,
                                                                              query_offset, operation,
                                                                              length, mm_tags, ml_tags)

            # get allele states
            read_idx_states = np.zeros(idx_var_b - idx_var_a + 1).astype(int) - 1

            for i, n in enumerate(self.phaseable[idx_var_a: idx_var_b + 1]):
                if self.idx_variant_mapping[n]['POS'] - 1 < read_start:
                    continue
                if self.idx_variant_mapping[n]['POS'] - 1 >= read_end:
                    break
                if str(n) not in self.read_states[read_name]:
                    # this can only really happen if a read only overlaps one variant and thus wasn't looked at before
                    # I think it might still be informative for methylation patterns though
                    cigar = deque(read.cigartuples)
                    gap_open, gap_extend, indel_rate, mismatch_rate = self.get_adaptive_gap_penalties(cigar)
                    pos = self.idx_variant_mapping[n]['POS'] - 1
                    read_start = read.reference_start
                    read_sequence = read.query_sequence
                    read_base_qualities = read.query_qualities
                    r_idx = 0
                    q_idx = 0
                    operation = None
                    length = 0
                    if i not in variant_homopolymer_mapping:
                        variant_homopolymer_mapping[i] = self.is_homopolymer(self.reference[pos - 5: pos + 5].seq.upper(),
                                                                             k=4)
                    state, qpos, _, _, _ = self.get_state_at_variant(read_sequence, read_base_qualities, cigar,
                                                                     read_start, r_idx, q_idx, operation, length,
                                                                     self.idx_variant_mapping[n])
                    if state is None:
                        state = self.realign_around_variant(read_sequence, self.idx_variant_mapping[n],
                                                            qpos, gap_open, gap_extend,
                                                            variant_homopolymer_mapping[n])
                        self.read_states[read_name][str(n)] = state
                    read_idx_states[i] = state
                else:
                    read_idx_states[i] = self.read_states[read_name][str(n)]
            # methylation_probs.append(read_methylation)
            methylation_probs['ml'].append(read_methylation[read_methylation != -1])
            methylation_probs['row'].append([len(methylation_probs['ml']) - 1] *
                                            read_methylation[read_methylation != -1].shape[0])
            methylation_probs['col'].append(np.where(read_methylation != -1)[0])
            read_variant_states.append(read_idx_states)
            read_ids.append(read_name)
        # subtract 1e-100 to get non-zero entries
        if len(methylation_probs) == 0:
            return None
        unmethylated_probs = csc_array((np.log10(1 - np.concatenate(methylation_probs['ml']) + 1e-20),
                                        (np.concatenate(methylation_probs['row']).astype(int),
                                         np.concatenate(methylation_probs['col']).astype(int))),
                                       shape=(len(methylation_probs['ml']), positions.shape[0]))
        methylation_probs = csc_array((np.log10(np.concatenate(methylation_probs['ml']) + 1e-20),
                                       (np.concatenate(methylation_probs['row']).astype(int),
                                        np.concatenate(methylation_probs['col']).astype(int))),
                                      shape=(len(methylation_probs['ml']), positions.shape[0]))

        read_variant_states = np.vstack(read_variant_states)

        # get last two allele states in last block
        unique_states_last_known, counts = np.unique(read_variant_states[(read_variant_states[:, 0] != -1) |
                                                                         (read_variant_states[:, 1] != -1), :2],
                                                     axis=0, return_counts=True)
        # determine the haplotypes
        if len(unique_states_last_known) == 1:
            hap1 = unique_states_last_known[np.argsort(counts)[-1]]
            # complement
            hap2 = np.where(hap1 == 0, 1, hap1)
            hap2 = np.where(hap1 == 1, 0, hap2)
        elif len(unique_states_last_known) > 1:
            hap1 = unique_states_last_known[np.argsort(counts)[-1]]
            idx_hap2 = -2
            hap2 = unique_states_last_known[np.argsort(counts)[idx_hap2]]
            while (((hap2[0] == hap1[0] and hap2[0] != -1) or (hap2[1] == hap1[1] and hap2[1] != -1))
                   and idx_hap2 - 1 >= -len(counts)):
                idx_hap2 -= 1
                hap2 = unique_states_last_known[np.argsort(counts)[idx_hap2]]
            # complement missing states
            hap1 = np.where((hap1 == -1) & (hap2 == 1), 0, hap1)
            hap1 = np.where((hap1 == -1) & (hap2 == 0), 1, hap1)
            hap2 = np.where((hap2 == -1) & (hap1 == 1), 0, hap2)
            hap2 = np.where((hap2 == -1) & (hap1 == 0), 1, hap2)
        # for some transitions reads cannot be anchored because no read with MM tag is overlapping the previous two
        # positions
        elif len(unique_states_last_known) == 0:
            return None

        # get the reads matching haplotype
        reads_hap1 = np.where(((read_variant_states[:, 0] == hap1[0]) & (read_variant_states[:, 1] == hap1[1])) |
                              ((read_variant_states[:, 0] == hap1[0]) & (hap1[0] != -1) &
                               (read_variant_states[:, 1] == -1)) |
                              ((read_variant_states[:, 0] == -1) & (read_variant_states[:, 1] == hap1[1]) &
                               (hap1[1] != -1)))[0]
        reads_hap2 = np.where(((read_variant_states[:, 0] == hap2[0]) & (read_variant_states[:, 1] == hap2[1])) |
                              ((read_variant_states[:, 0] == hap2[0]) & (hap2[0] != -1) &
                               (read_variant_states[:, 1] == -1)) |
                              ((read_variant_states[:, 0] == -1) & (read_variant_states[:, 1] == hap2[1]) &
                               (hap2[1] != -1)))[0]

        unassigned = np.where((read_variant_states[:, 0] == -1) & (read_variant_states[:, 1] == -1))[0]
        prev_unassigned = []

        while len(unassigned) > 0 and not np.array_equal(unassigned, prev_unassigned):
            meth_hap1 = methylation_probs[reads_hap1].sum(axis=0)
            unmeth_hap1 = unmethylated_probs[reads_hap1].sum(axis=0)
            # methylated sites must have a probability greater than 0.5 and an LLR > 3
            meth_states_hap1 = np.where((meth_hap1 > np.log10(0.5)) & (meth_hap1 - unmeth_hap1 > 3), 1, 0)
            # unmethylated sites must have a probability greater than 0.5 and an LLR < -3
            # meth_states_hap1 = np.where((unmeth_hap1 > np.log10(0.5)) & (unmeth_hap1 - meth_hap1 > 3), 0, meth_states_hap1)

            meth_hap2 = methylation_probs[reads_hap2].sum(axis=0)
            unmeth_hap2 = unmethylated_probs[reads_hap2].sum(axis=0)

            # methylated sites must have a probability greater than 0.5 and an LLR > 3
            meth_states_hap2 = np.where((meth_hap2 > np.log10(0.5)) & (meth_hap2 - unmeth_hap2 > 3), 1, 0)
            # unmethylated sites must have a probability greater than 0.5 and an LLR < -3
            # meth_states_hap2 = np.where((unmeth_hap2 > np.log10(0.5)) & (unmeth_hap2 - meth_hap2 > 3), 0, meth_states_hap2)

            # find differentially methylated sites
            diff_meth = np.where((meth_states_hap1 != meth_states_hap2) & # hap1 and hap2 must have different state
                                 (meth_states_hap1 >= 0) &  # hap1 must be defined
                                 (meth_states_hap2 >= 0))[0]  # hap1 must be defined

            methylation_per_read_hap1 = self.get_diff_methylation_sites_per_hap(c_methylation_calls, methylation_probs,
                                                                                reads_hap1, "hap1")
            methylation_per_read_hap2 = self.get_diff_methylation_sites_per_hap(c_methylation_calls, methylation_probs,
                                                                                reads_hap2, "hap2")

            self.differentially_methylated_sites.append(pd.concat([methylation_per_read_hap1.iloc[diff_meth],
                                                                   methylation_per_read_hap2.iloc[diff_meth],
                                                                   c_methylation_calls.iloc[diff_meth]]))
            # calculate probabilities that methylation pattern of read match either haplotype
            p_hap1 = (
                self.calculate_probability_of_reads_belonging_to_haplotype_based_on_methylation(meth_states_hap1,
                                                                                                meth_hap1, unmeth_hap1,
                                                                                                methylation_probs.copy(),
                                                                                                unmethylated_probs.copy(),
                                                                                                unassigned, diff_meth))
            p_hap2 = (
                self.calculate_probability_of_reads_belonging_to_haplotype_based_on_methylation(meth_states_hap2,
                                                                                                meth_hap2, unmeth_hap2,
                                                                                                methylation_probs.copy(),
                                                                                                unmethylated_probs.copy(),
                                                                                                unassigned, diff_meth))

            prev_unassigned = unassigned.copy()
            # assign reads based A > B and A > 0.5
            reads_hap1 = np.concatenate([reads_hap1, np.array(unassigned)[(p_hap1 > p_hap2) &
                                                                          (p_hap1 > np.log10(0.5))]])
            reads_hap2 = np.concatenate([reads_hap2, np.array(unassigned)[(p_hap2 > p_hap1) &
                                                                          (p_hap2 > np.log10(0.5))]])
            unassigned = np.array(unassigned)[~((p_hap1 > p_hap2) & (p_hap1 > np.log10(0.5))) &
                                              ~((p_hap2 > p_hap1) & (p_hap2 > np.log10(0.5)))]
        self.methylation_read_assignments['hap1'].extend([(read_ids[i], (idx_var_a, hap1[0]), (idx_var_b, hap1[-1]))
                                                          for i in reads_hap1])
        self.methylation_read_assignments['hap2'].extend([(read_ids[i], (idx_var_a, hap2[0]), (idx_var_b, hap2[-1]))
                                                          for i in reads_hap2])
        # calculate most common variant states at sites flanking uncertain transition
        read_variant_states = np.where(read_variant_states == -1, np.nan, read_variant_states)
        read_variant_cov = np.where(np.isnan(read_variant_states), np.nan, 1)
        hap1 = np.nanmean(read_variant_states[reads_hap1], axis=0)
        hap2 = np.nanmean(read_variant_states[reads_hap2], axis=0)
        hap1_cov = np.nansum(read_variant_cov[reads_hap1], axis=0)
        hap2_cov = np.nansum(read_variant_cov[reads_hap2], axis=0)

        # complement missing states
        hap1 = np.where(np.isnan(hap1) & (hap2 == 1), 0, hap1)
        hap1 = np.where(np.isnan(hap1) & (hap2 == 0), 1, hap1)
        hap2 = np.where(np.isnan(hap2) & (hap1 == 1), 0, hap2)
        hap2 = np.where(np.isnan(hap2) & (hap1 == 0), 1, hap2)

        # update transition matrix at uncertain transition
        # --> fraction of sites in state X at site A * fraction of sites in state Y at site B
        # get all transitions
        c_transitions = self.phaseable[idx_var_a: idx_var_b + 1]
        # number of uncertain transitions
        n_uncertain_transitions = idx_var_b - idx_var_a - 3
        offset = n_uncertain_transitions + 1

        for i in range(2, offset + 1):
            self.calculate_transition_probability_from_methylation(hap1, hap2, hap1_cov, hap2_cov, i, c_transitions[i])

    def get_methylation_transitions_helper(self):
        """
        Identify transitions that are ambiguous based on variant data and try to fill in transitions based on
        methylation data.
        """
        uncertain_transitions = self.get_uncertain_transitions(self.transition_matrix[:, :,
                                                               self.phaseable[self.phaseable < self.num_variants - 2]])
        n_uncertain_transitions = len(uncertain_transitions)
        pbar = tqdm(total=len(uncertain_transitions))
        while uncertain_transitions:
            # include last certain transition
            i = uncertain_transitions.popleft()
            idx_var_a = i - 2
            # search for next certain transition
            while uncertain_transitions and uncertain_transitions[0] == i + 1:
                i = uncertain_transitions.popleft()
                pbar.update(1)
            idx_var_b = i + 2
            if idx_var_a < 0 or idx_var_b > self.phaseable.shape[0] - 1:
                continue
            # Not saving previously observed read methylations saves a ton of memory and may actually run faster
            self.get_methylation_transitions(idx_var_a, idx_var_b)
            pbar.update(1)
        self.transition_matrix /= self.transition_matrix.sum(axis=1, keepdims=True)

        uncertain_transitions = self.get_uncertain_transitions(self.transition_matrix[:, :,
                                                               self.phaseable[self.phaseable < self.num_variants - 2]])
        logging.info(f'Filled in {n_uncertain_transitions - len(uncertain_transitions)}/{n_uncertain_transitions} ' +
                     'uncertain transitions with methylation data.')
        if self.output_transition_matrix_meth is not None:
            np.savez(self.output_transition_matrix_meth, self.transition_matrix)


    @staticmethod
    def get_uncertain_transitions(transition_matrix):
        """
        Get indices of transitions that are ambiguous
        :return: list-like, indices of uncertain transitions
        """
        uncertain_transitions = []
        for i in range(transition_matrix.shape[2]):
            if ((transition_matrix[0, :, i].min() == 0.5 and transition_matrix[0, :, i].max() == 0.5) and
                    (transition_matrix[1, :, i].min() == 0.5 and transition_matrix[1, :, i].max() == 0.5)):
                uncertain_transitions.append(i)
        uncertain_transitions = deque(uncertain_transitions)
        return uncertain_transitions

    def update_transition_matrix_considering_adjacent_variants(self, idx_a, idx_b, idx_c, min_cov=1):
        """
        Re-calculate transition probabilities between variant A and B and B and C using "long-range"
        information between variant A and C
        :param idx_a: int, index of variant A
        :param idx_b: int, index of variant B
        :param idx_c: int, index of variant C
        :param min_cov: int, minimum coverage to consider transition between variant A and C informative
        :return: (np.array, np.array, np.array, np.array, np.array),
                    (original transition between A and B, original transition between B and C,
                    new transition between A and B, new transition between B and C, transition between A and C)
        """
        t3 = np.zeros((2, 2)) + 1e-20
        new_t1 = np.zeros((2, 2)) + 1e-20
        new_t2 = np.zeros((2, 2)) + 1e-20
        # transitions cannot be updated
        if (idx_a < 0 or idx_c > self.num_variants - 2 or
                # self.variant_type[idx_a - 1] != 'SNP' or self.variant_type[idx_a + 1] != 'SNP' or
                (np.unique(self.transition_matrix[:, :, idx_a].argmax(axis=1)).shape[0] == 2 and
                 np.unique(self.transition_matrix[:, :, idx_b].argmax(axis=1)).shape[0] == 2)):
            t1 = np.zeros((2, 2)) + 1e-20
            t2 = np.zeros((2, 2)) + 1e-20
            t1 = t1 / t1.sum(axis=1, keepdims=True)
            t2 = t2 / t2.sum(axis=1, keepdims=True)
            t3 = t3 / t3.sum(axis=1, keepdims=True)
            return t1, t2, new_t1, new_t2, t3
        # get tranisiton between A and C
        t3 = self.get_allele_transitions_from_known_read_states(idx_a, idx_c)

        t1 = self.transition_matrix[:, :, idx_a]
        t2 = self.transition_matrix[:, :, idx_b]
        t1 = t1 / t1.sum(axis=1, keepdims=True)
        t2 = t2 / t2.sum(axis=1, keepdims=True)
        # check for minimal coverage
        if t3.sum(axis=1).max() >= min_cov:

            t3 = self.mirror_transition(t3, normalized=False)

            t3 = t3 / t3.sum(axis=1, keepdims=True)

            # paths directly involving A1 to A2 divided by all possible paths exciting A1
            new_t1[0, 0] = (t1[0, 0] * t2[0, 0] * t3[0, 0] + t1[0, 0] * t2[0, 1] * t3[0, 1]) / \
                           (t1[0, 0] * t2[0, 0] * t3[0, 0] + t1[0, 0] * t2[0, 1] * t3[0, 1] +
                            t1[0, 1] * t2[1, 0] * t3[0, 0] + t1[0, 1] * t2[1, 1] * t3[0, 1])

            # A1 to B2 must be the complement of A1 to A2
            new_t1[0, 1] = (t1[0, 1] * t2[1, 0] * t3[0, 0] + t1[0, 1] * t2[1, 1] * t3[0, 1]) / \
                           (t1[0, 0] * t2[0, 0] * t3[0, 0] + t1[0, 0] * t2[0, 1] * t3[0, 1] +
                            t1[0, 1] * t2[1, 0] * t3[0, 0] + t1[0, 1] * t2[1, 1] * t3[0, 1])
            # paths directly involving the edge from B1 to B2 divided by all possible paths exciting B1
            new_t1[1, 1] = (t1[1, 1] * t2[1, 0] * t3[1, 0] + t1[1, 1] * t2[1, 1] * t3[1, 1]) / \
                           (t1[1, 1] * t2[1, 0] * t3[1, 0] + t1[1, 1] * t2[1, 1] * t3[1, 1] +
                            t1[1, 0] * t2[0, 0] * t3[1, 0] * t1[1, 0] * t2[0, 1] * t3[1, 1])
            # B1 to A2 must be the complement of B1 to B2
            new_t1[1, 0] = (t1[1, 0] * t2[0, 0] * t3[1, 0] + t1[1, 0] * t2[0, 1] * t3[1, 1]) / \
                           (t1[1, 1] * t2[1, 0] * t3[1, 0] + t1[1, 1] * t2[1, 1] * t3[1, 1] +
                            t1[1, 0] * t2[0, 0] * t3[1, 0] * t1[1, 0] * t2[0, 1] * t3[1, 1])

            # paths directly involving A2 to A3 divided by all paths traversing A2
            new_t2[0, 0] = (t1[0, 0] * t2[0, 0] * t3[0, 0] + t1[1, 0] * t2[0, 0] * t3[1, 0]) / \
                           (t1[0, 0] * t2[0, 0] * t3[0, 0] + t1[1, 0] * t2[0, 0] * t3[1, 0] +
                            t1[0, 0] * t2[0, 1] * t3[0, 1] + t1[1, 0] * t2[0, 1] * t3[1, 1])
            new_t2[0, 1] = (t1[0, 0] * t2[0, 1] * t3[0, 1] + t1[1, 0] * t2[0, 1] * t3[1, 1]) / \
                           (t1[0, 0] * t2[0, 0] * t3[0, 0] + t1[1, 0] * t2[0, 0] * t3[1, 0] +
                            t1[0, 0] * t2[0, 1] * t3[0, 1] + t1[1, 0] * t2[0, 1] * t3[1, 1])
            # paths directly involving the edge from B2 to B3 divided by all paths traversing B2
            new_t2[1, 1] = (t1[1, 1] * t2[1, 1] * t3[1, 1] + t1[0, 1] * t2[1, 1] * t3[0, 1]) / \
                           (t1[1, 1] * t2[1, 1] * t3[1, 1] + t1[0, 1] * t2[1, 1] * t3[0, 1] +
                            t1[1, 1] * t2[1, 0] * t3[1, 0] + t1[0, 1] * t2[1, 0] * t3[0, 0])
            new_t2[1, 0] = (t1[1, 1] * t2[1, 0] * t3[1, 0] + t1[0, 1] * t2[1, 0] * t3[0, 0]) / \
                           (t1[1, 1] * t2[1, 1] * t3[1, 1] + t1[0, 1] * t2[1, 1] * t3[0, 1] +
                            t1[1, 1] * t2[1, 0] * t3[1, 0] + t1[0, 1] * t2[1, 0] * t3[0, 0])
            new_t1 += 1e-20
            new_t2 += 1e-20
        else:
            t3 = np.zeros((2, 2)) + 1e-20
            t3 = t3 / t3.sum(axis=1, keepdims=True)

        return t1, t2, new_t1, new_t2, t3

    def get_new_transition_conditional_on_beliefs(self, beliefs, layer_idx, layer_idx_mapping, n_preceding, n_succeeding):
        """
        Update transition matrix based on inferred beliefs at each variant position
        :param beliefs: np.array, inferred beliefs at each variant position
        :param layer_idx: np.array, variant indices associated with layers
        :param layer_idx_mapping: dict, mapping from variant index to layer index
        :param n_preceding: int, number of preceding variants to consider
        :param n_succeeding: int, number of succeeding variants to consider
        """
        n_preceding = max(n_preceding, 1)
        n_succeeding = max(n_succeeding, 1)
        parents = layer_idx[n_preceding - 1: layer_idx.shape[0] - n_succeeding]
        children = layer_idx[n_preceding: layer_idx.shape[0] - n_succeeding + 1]

        parents_idx = np.array([layer_idx_mapping[p] for p in parents])
        children_idx = np.array([layer_idx_mapping[c] for c in children])

        beliefs_parents = beliefs[parents_idx]
        beliefs_children = beliefs[children_idx]

        new_transitions = np.ones((2, 2, len(parents)))

        mask = np.abs(beliefs_parents[:, 0] - 0.5) != 0.0

        t_rr = np.maximum((beliefs_parents[mask, 0] + beliefs_children[mask, 0] - 1) /
                          (2 * beliefs_parents[mask, 0] - 1), 0)
        t_ra = np.maximum((beliefs_parents[mask, 0] + beliefs_children[mask, 1] - 1) /
                          (2 * beliefs_parents[mask, 0] - 1), 0)
        t_ar = np.maximum((beliefs_parents[mask, 1] + beliefs_children[mask, 0] - 1) /
                          (2 * beliefs_parents[mask, 1] - 1), 0)
        t_aa = np.maximum((beliefs_parents[mask, 1] + beliefs_children[mask, 1] - 1) /
                          (2 * beliefs_parents[mask, 1] - 1), 0)

        new_transitions[0, 0, mask] = t_rr
        new_transitions[0, 1, mask] = t_ra
        new_transitions[1, 0, mask] = t_ar
        new_transitions[1, 1, mask] = t_aa

        new_transitions = new_transitions / new_transitions.sum(axis=1, keepdims=True) + 1e-300
        new_transitions = np.nan_to_num(new_transitions, nan=0.5)
        new_transitions *= np.mean(self.allele_coverage)
        for i in range(new_transitions.shape[2]):
            new_transitions[:, :, i] = self.mirror_transition(new_transitions[:, :, i], normalized=False)
        self.transition_matrix[:, :, parents] = new_transitions

    def loopy_belief_propagation(self, layer_idx, n_preceding, n_succeeding, normalized=False, damping=0.2, tol=1e-10,
                                 max_iters=500):
        """
        Perform loopy belief propagation to infer marginal probabilities at each variant position in a subgraph of
        difficult to phase variants
        :param layer_idx: np.array, variant indices associated with layers
        :param n_preceding: int, number of preceding variants to consider
        :param n_succeeding: int, number of succeeding variants to consider
        :param normalized: boolean, transitions have been normalized
        :param damping: float, damping factor for message updates
        :param tol: float, tolerance for convergence
        :param max_iters: int, maximum number of iterations
        """
        n_layers = layer_idx.shape[0]
        # --- Build parents/children and log-factors ------------------------------------------------
        parents = defaultdict(list)
        children = defaultdict(list)
        neighbors = defaultdict(set)  # undirected graph

        # log_factors = np.zeros((n_layers, n_layers, 2, 2))
        log_factors = {}

        edges = []
        coverage = {}

        layer_idx_mapping = {layer_idx[i]: i for i in range(len(layer_idx))}
        # direct adjacent transitions
        for i, idx_a in enumerate(layer_idx[:-1]):
            t = self.transition_matrix[:, :, idx_a]
            t = self.mirror_transition(t, normalized=normalized)
            # normalize rows in linear space, then take log
            cov = t.sum(axis=1).max()
            t = t / t.sum(axis=1, keepdims=True) + 1e-300

            log_factors[(i, i + 1)] = np.log(t)
            log_factors[(i + 1, i)] = np.log(t).T

            parents[i + 1].append(i)
            neighbors[i].add(i + 1)
            neighbors[i + 1].add(i)
            edges.append((i, i + 1))
            edges.append((i + 1, i))
            coverage[(i, i + 1)] = cov
            coverage[(i + 1, i)] = cov

        # add long-range upstream transitions
        for i, idx_a in enumerate(layer_idx[:n_preceding]):
            for idx_b in layer_idx[i + 2:]:
                t = self.get_allele_transitions_from_known_read_states(idx_a, idx_b)
                t = self.mirror_transition(t, normalized=False)
                cov = t.sum(axis=1).max()
                t = t / t.sum(axis=1, keepdims=True) + 1e-300
                # skip trivial 0.5 matrix if desired
                if np.allclose(t, 0.5):
                    continue
                ix_b = layer_idx_mapping[idx_b]

                log_factors[(i, ix_b)] = np.log(t)
                log_factors[(ix_b, i)] = np.log(t).T

                parents[ix_b].append(i)
                neighbors[i].add(ix_b)
                neighbors[ix_b].add(i)
                edges.append((i, ix_b))
                edges.append((ix_b, i))
                edges.append((ix_b, i))
                coverage[(i, ix_b)] = cov
                coverage[(ix_b, i)] = cov

        # add long-range downstream transitions
        for idx_a in layer_idx[n_preceding:-n_succeeding] if n_succeeding > 0 else layer_idx[n_preceding:]:
            for idx_b in layer_idx[-n_succeeding:] if n_succeeding > 0 else []:
                # ensure at least distance 2 (your previous guard)
                if (idx_b - idx_a) < 2:
                    continue
                t = self.get_allele_transitions_from_known_read_states(idx_a, idx_b)
                t = self.mirror_transition(t, normalized=False)
                cov = t.sum(axis=1).max()
                t = t / t.sum(axis=1, keepdims=True) + 1e-300
                if np.allclose(t, 0.5):
                    continue
                ix_a = layer_idx_mapping[idx_a]
                ix_b = layer_idx_mapping[idx_b]

                log_factors[(ix_a, ix_b)] = np.log(t)
                log_factors[(ix_b, ix_a)] = np.log(t).T

                parents[ix_b].append(ix_a)
                neighbors[ix_a].add(ix_b)
                neighbors[ix_b].add(ix_a)
                edges.append((ix_a, ix_b))
                edges.append((ix_b, ix_a))
                coverage[(ix_a, ix_b)] = cov
                coverage[(ix_b, ix_a)] = cov

        # build children mapping
        for child, ps in parents.items():
            for p in ps:
                children[p].append(child)

        # unary log-potentials: root fixed to REF (if you want)
        phi = np.log(np.ones((n_layers, 2)))
        phi[0, :] = np.array([0, -np.inf])

        # initialize messages (i->j) as log(1)=0
        messages = np.zeros((n_layers, n_layers, 2))
        new_messages = np.zeros((n_layers, n_layers, 2))
        edges = np.asarray(edges, dtype=np.int32)
        edges_i = edges[:, 0]
        edges_j = edges[:, 1]
        median_coverage =  np.median([v for v in coverage.values()])
        factors = np.stack([(coverage[i, j] / (coverage[i, j] + median_coverage)) * log_factors[i, j]
                            for i, j in zip(edges_i, edges_j)], axis=0)

        for it in range(max_iters):
            new_messages[:] = 0.0
            total_incoming = messages.sum(axis=0)
            total = phi + total_incoming

            incoming = total[edges_i] - messages[edges_j, edges_i, :]
            m_vals = logsumexp(incoming[:, :, None] + factors, axis=1)
            m_vals -= logsumexp(m_vals, axis=1)[:, None]

            if damping:
                m_vals = damping * messages[edges_i, edges_j] + (1 - damping) * m_vals

            new_messages[edges_i, edges_j] = m_vals
            delta = np.max(np.abs(new_messages[edges_i, edges_j] - messages[edges_i, edges_j]))

            messages, new_messages = new_messages, messages
            if delta < tol:
                break

        # compute node beliefs and pairwise marginals
        total_incoming = messages.sum(axis=0)
        beliefs = phi + total_incoming
        beliefs -= logsumexp(beliefs, axis=1, keepdims=True)
        beliefs = np.exp(beliefs)

        # compute posterior transition probabilities and DAG with long-range edges projected into pairwise transitions
        self.get_new_transition_conditional_on_beliefs(beliefs, layer_idx, layer_idx_mapping, n_preceding,
                                                       n_succeeding)

    def rephase_difficult_variants(self, n_preceding=2, n_succeeding=2, normalized=False, damping=0.0):
        """
        Rephase difficult variants considering long-range phase information of adjacent variants upstream and downstream.
        :param n_preceding: int, number of upstream variants to consider
        :param n_succeeding: int, number of downstream variants to consider
        """
        # indices of all difficult variants
        vars_to_rephase = np.where((self.variant_type[self.phaseable] != 'SNP') |
                                   (self.allele_coverage[:, self.phaseable].min(axis=0) < self.min_allele_count))[0]
        p_idx_a = -1
        # find first difficult variant
        for idx_a in tqdm(vars_to_rephase):
            if idx_a <= p_idx_a:
                continue
            first_var = idx_a
            last_var = idx_a
            if first_var < n_preceding:
                continue
            # find last difficult variant in consecutive series
            while (last_var + 1 < self.phaseable.shape[0] - n_succeeding - 2 and
                   ((self.variant_type[self.phaseable[last_var + 1]] != 'SNP') or
                    (self.allele_coverage[:, self.phaseable[last_var + 1]].min() < self.min_allele_count))):
                last_var += 1
            if last_var > self.phaseable.shape[0] - n_succeeding - 2:
                continue
            layer_idx = np.array([self.phaseable[i] for i in range(first_var - n_preceding,
                                                                   last_var + n_succeeding + 1)])
            self.loopy_belief_propagation(layer_idx, n_preceding, n_succeeding, normalized=normalized,
                                          damping=damping)
            # self.loopy_belief_propagation_sparse(layer_idx, n_preceding, n_succeeding, normalized=False, damping=0.0)
            p_idx_a = last_var

    def connect_phase_blocks(self):
        """
        Connect phase block by checking if a transition between the last variant of block A and the first variant of
        block B can be established. All variants inbetween will be marked as unphaseable.
        """
        norm_trans_mat = self.transition_matrix / self.transition_matrix.sum(axis=1, keepdims=True)
        uncertain_transitions = (
            self.get_uncertain_transitions(norm_trans_mat[:, :,
                                           self.phaseable[self.phaseable < self.num_variants - 2]]))
        new_unphaseable = []
        n_blocks = 1
        connected = 0
        while uncertain_transitions:
            # last variant of previous phase block
            idx_a = uncertain_transitions.popleft()

            idx_b = idx_a + 1
            c_unphaseable = []
            n_blocks += 1
            # find first variant of next phaseblock
            while len(uncertain_transitions) > 0 and idx_b == uncertain_transitions[0]:
                if idx_b >= self.phaseable.shape[0]:
                    break
                c_unphaseable.append(self.phaseable[idx_b])
                idx_b += 1
                uncertain_transitions.popleft()
            # check if last and first variant are within read length
            if idx_b - idx_a == 1 and idx_b < self.phaseable.shape[0]:
                c_unphaseable.append(self.phaseable[idx_b])
                idx_b += 1
            idx_a_p = idx_a
            idx_b_p = idx_b
            if idx_b >= self.phaseable.shape[0]:
                continue
            idx_a = self.phaseable[idx_a]
            idx_b = self.phaseable[idx_b]
            if idx_b < self.num_variants and idx_b_p - idx_a_p > 2:
                t = self.get_allele_transitions_from_known_read_states(idx_a, idx_b)
                t = self.mirror_transition(t, normalized=False)

                # can connect the blocks
                # TODO only having the first conditions gives a lower error
                if ((t / t.sum(axis=1, keepdims=True)).max() > 0.5 and
                        self.transition_matrix[:, :, self.phaseable[idx_b_p - 1]].max() == 0.5):
                    self.transition_matrix[:, :, idx_a] = t
                    new_unphaseable.extend(c_unphaseable)
                    connected += 1
            elif idx_b_p - idx_a_p == 2 and 0 < idx_a < self.num_variants - 2 and idx_a_p > 0:
                t1, t2, new_t1, new_t2, t3 = self.update_transition_matrix_considering_adjacent_variants(
                    self.phaseable[idx_a_p - 1], idx_a, self.phaseable[idx_a_p + 1])
                new_t1 = new_t1 / new_t1.sum(axis=1, keepdims=True)
                new_t2 = new_t2 / new_t2.sum(axis=1, keepdims=True)
                new_t1 = self.mirror_transition(new_t1, normalized=True)
                new_t2 = self.mirror_transition(new_t2, normalized=True)
                # t1 --> idx_a - 1 to idx_a
                # t2 --> idx_a to idx_a + 1
                # t3 --> idx_a - 1 to idx_a + 1
                # new_t1 and new_t2 are uninformative but t3 is informative --> use t3 to connect idx_a - 1 and idx_a + 1
                # and mark idx_a as uninphaseable
                if ((np.all(new_t1) == 0.5 or np.all(new_t2) == 0.5 or
                     (np.unique(new_t1.argmax(axis=1)).shape[0] == 1 and new_t1.max(axis=1).min() > 0) or
                     (np.unique(new_t2.argmax(axis=1)).shape[0] == 1 and new_t2.max(axis=1).min() > 0) or
                     new_t2[new_t1[0, :].argmax(), :].argmax() != t3.argmax()) and
                        t3.max() > 0.5):
                    self.transition_matrix[:, :, self.phaseable[idx_a_p - 1]] = t3
                    new_unphaseable.append(idx_a)
                    connected += 1
                # new t1 and t2 are informative and in agreement with t3 --> use new t1 and new t2 to connect
                elif (new_t1.max() > 0.5 and new_t2.max() > 0.5 and
                      new_t2[new_t1[0, :].argmax(), :].argmax() == t3.argmax()):
                    self.transition_matrix[:, :, self.phaseable[idx_a_p - 1]] = new_t1
                    self.transition_matrix[:, :, idx_a] = new_t2
                    connected += 1
                # else cannot connect the blocks
                # consider idx_a --> idx_a +1 --> idx_a + 2
                elif 0 < idx_a < self.num_variants - 3:
                    t1, t2, new_t1, new_t2, t3 = self.update_transition_matrix_considering_adjacent_variants(
                        idx_a, self.phaseable[idx_a_p + 1], self.phaseable[idx_a_p + 2])
                    new_t1 = new_t1 / new_t1.sum(axis=1, keepdims=True)
                    new_t2 = new_t2 / new_t2.sum(axis=1, keepdims=True)
                    new_t1 = self.mirror_transition(new_t1, normalized=True)
                    new_t2 = self.mirror_transition(new_t2, normalized=True)

                    if ((np.all(new_t1) == 0.5 or np.all(new_t2) == 0.5 or
                         (np.unique(new_t1.argmax(axis=1)).shape[0] == 1 and new_t1.max(axis=1).min() > 0) or
                         (np.unique(new_t2.argmax(axis=1)).shape[0] == 1 and new_t2.max(axis=1).min() > 0) or
                         new_t2[new_t1[0, :].argmax(), :].argmax() != t3.argmax()) and
                            t3.max() > 0.5):
                        self.transition_matrix[:, :, idx_a] = t3
                        new_unphaseable.append(self.phaseable[idx_a_p + 1])
                        connected += 1
                    elif (new_t1.max() > 0.5 and new_t2.max() > 0.5 and
                          new_t2[new_t1[0, :].argmax(), :].argmax() == t3.argmax()):

                        self.transition_matrix[:, :, idx_a] = new_t1
                        self.transition_matrix[:, :, self.phaseable[idx_a_p + 1]] = new_t2
                        connected += 1

        self.phaseable = self.phaseable[~np.isin(self.phaseable, new_unphaseable)]
        self.unphaseable = np.sort(np.concatenate([new_unphaseable, self.unphaseable])).astype(int)
        self.transition_matrix /= self.transition_matrix.sum(axis=1, keepdims=True)

        logging.info(f'Connected {connected}/{n_blocks} phase blocks, leaving {len(new_unphaseable)} '
                     f'variants unphased.')

    @staticmethod
    def is_homopolymer(seq, k=4):
        """
        Check if a sequence is a homopolymer of length k or longer
        :param seq: str, sequence
        :param k: int, minimum homopolymer length
        :return: boolean, True if homopolymer of length k or longer
        """
        if len(seq) < k:
            return False
        count = 1
        prev_base = seq[0]
        for base in seq[1:]:
            if base == prev_base:
                count += 1
                if count >= k:
                    return True
            else:
                count = 1
                prev_base = base
        return False

    def create_directed_graph_of_heterozygous_variants_from_reads(self):
        """
        Create phase transition matrix based on overlap in heterozygous variants
        """
        min_var_idx = 0
        warned = False
        # iterate over all reads
        if self.seqtech == 'ont':
            strands = np.zeros((2, 2, self.num_variants))
        variant_homopolymer_mapping = {}
        for read in tqdm(self.alignments.fetch(self.chrom), total=self.alignments.count(self.chrom)):
            if (read.is_secondary or read.is_duplicate or read.is_unmapped or read.is_qcfail or
                    read.mapping_quality < self.min_mapq):
                continue
            prev_state = -1
            cigar = deque(read.cigartuples)
            gap_open, gap_extend, indel_rate, mismatch_rate = self.get_adaptive_gap_penalties(cigar)
            read_start = read.reference_start
            read_end = read.reference_end
            read_name = read.query_name
            read_sequence = read.query_sequence
            read_base_qualities = read.query_qualities
            r_idx = 0
            q_idx = 0
            operation = None
            length = 0
            found_first_var = False
            query_idx_vars = {}
            if min_var_idx > 0 and read_start < self.idx_variant_mapping[min_var_idx - 1]["POS"] - 1 and not warned:
                warned = True
                logging.warning('BAM is not sorted by position. Suppressing future warnings.')
            # iterate over all variants starting with the first variant that was covered by the previous read
            # --> assumes bam is position sorted
            for i in range(min_var_idx, self.num_variants):
                var = self.idx_variant_mapping[i]
                pos = var["POS"]
                # variant is not yet covered by read
                if pos - 1 < read_start:
                    continue
                # variant is beyond read end
                if pos - 1 >= read_end:
                    break
                if not found_first_var:
                    found_first_var = True
                    min_var_idx = i
                if i not in variant_homopolymer_mapping:
                    variant_homopolymer_mapping[i] = self.is_homopolymer(self.reference[pos - 5: pos + 5].seq.upper(),
                                                                         k=4)
                state, q_idx, r_idx, operation, length = self.get_state_at_variant(read_sequence,
                                                                                   read_base_qualities, cigar,
                                                                                   read_start, r_idx,
                                                                                   q_idx, operation, length, var)
                # fill transition matrix
                if state != -1 and state is not None:
                    self.allele_coverage[state, i] += 1
                    if self.seqtech == 'ont' and read.is_forward:
                        strands[state, 0, i] += 1
                    elif self.seqtech == 'ont' and not read.is_forward:
                        strands[state, 1, i] += 1
                    if prev_state != -1 and prev_state is not None:
                        self.transition_matrix[prev_state, state, i - 1] += 1

                self.read_states[read_name][str(i)] = state
                self.variant_read_mapping[str(i)].append(read_name)
                query_idx_vars[str(i)] = q_idx
                prev_state = state

            # realign uncertain variants
            idx_to_realign = [int(n) for n, s in self.read_states[read_name].items() if s is None]
            for n in sorted(idx_to_realign):
                state = self.realign_around_variant(read_sequence, self.idx_variant_mapping[n],
                                                    query_idx_vars[str(n)], gap_open, gap_extend, variant_homopolymer_mapping[n])
                self.read_states[read_name][str(n)] = state
                if state != -1:
                    if self.seqtech == 'ont' and read.is_forward:
                        strands[state, 0, n] += 1
                    elif self.seqtech == 'ont' and not read.is_forward:
                        strands[state, 1, n] += 1

                    self.allele_coverage[state, n] += 1

                    if str(n - 1) in self.read_states[read_name] and self.read_states[read_name][str(n - 1)] != -1:
                        self.transition_matrix[self.read_states[read_name][str(n - 1)], state, n - 1] += 1
        if self.seqtech == 'ont':
            strands /= strands.sum(axis=1, keepdims=True)
            self.unphaseable = np.unique(np.sort(np.where(strands == 0)[2]))
            self.phaseable = self.phaseable[~np.isin(self.phaseable, self.unphaseable)]
            logging.info(f'Marked {self.unphaseable.shape[0]} variants as not phaseable '
                         f'as they are likely false positive calls only found on one strand')

            i = 0
            while i < self.phaseable.shape[0] - 1:
                if self.phaseable[i + 1] - self.phaseable[i] > 1:
                    self.transition_matrix[:, :,
                        self.phaseable[i]] = self.get_allele_transitions_from_known_read_states(self.phaseable[i],
                                                                                                self.phaseable[i + 1])
                i += 1

    @staticmethod
    def backtrace(delta, hap, phi, last_position, first_position):
        """
        Get most likely haplotype (sequence of phases/states) for current block
        :param delta: np.array, delta matrix with likelihoods
        :param hap: np.array, phased haplotype
        :param phi: np.array, point to previous state with the highest likelihood
        :param last_position: int, end position of current block in hap (backtracing starts here)
        :param first_position: int, start position of current block in hap (backtracing ends here)
        :return: np.array, modified version of hap with haplotypes information for current block
        """
        hap[last_position] = np.argmax(delta[:, last_position])
        for m in range(last_position - 1, first_position - 1, -1):
            hap[m] = phi[hap[m + 1], m + 1]
        return hap

    def calculate_forward_path_probabilities(self):
        """
        Find haplotypes with the highest likelihood for each haplotype block. We use a Viterbi-like approach
        :return: np.array, list, np.array: (2, L) phases for the L heterozygous variants,
                                            end indices for each haplotype block,
                                            and state probabilities at each variants
        """
        # pick sites that are can be phased
        transition_matrix = self.transition_matrix[:, :, self.phaseable[:-1]]
        delta = np.zeros((2, transition_matrix.shape[2] + 1)) - 1
        # argmax for each state
        phi = np.zeros((2, transition_matrix.shape[2] + 1)).astype(int) - 1

        # initialize paths
        hap = np.zeros(transition_matrix.shape[2] + 1).astype(int) - 1
        start_prev_block = 0

        # initial probabilities
        delta[:, 0] = np.log([0.5, 0.5])
        phi[:, 0] = [0, 1]
        transition_matrix = np.log(np.where(transition_matrix == 0, 1e-20, transition_matrix))
        for l in range(1, delta.shape[1]):
            if (transition_matrix[0, 0, l - 1] == transition_matrix[0, 1, l - 1] and
                    transition_matrix[1, 0, l - 1] == transition_matrix[1, 1, l - 1]):

                if start_prev_block == l - 1:
                    start_prev_block = l
                else:
                    # backtrace current block
                    hap = self.backtrace(delta, hap, phi, l - 1, start_prev_block)
                    start_prev_block = l
                    self.block_ends.append(self.phaseable[l - 1])

                # initialize new block
                delta[:, l] = np.log([0.5, 0.5])
                phi[:, l] = [0, 1]
                continue

            # calculate path probabilities
            delta[0, l] = np.max(delta[:, l - 1] + transition_matrix[0, :, l - 1])
            delta[1, l] = np.max(delta[:, l - 1] + transition_matrix[1, :, l - 1])

            # record pointer
            phi[0, l] = np.argmax(delta[:, l - 1] + transition_matrix[0, :, l - 1])
            phi[1, l] = np.argmax(delta[:, l - 1] + transition_matrix[1, :, l - 1])

        # backtrace last block
        hap = self.backtrace(delta, hap, phi, l, start_prev_block)
        # hap = np.argmax(delta, axis=0)
        self.block_ends.append(self.phaseable[l])
        hap_0 = hap.copy()
        hap_1 = np.where(hap_0 == 0, 1, 0)
        hap_1 = np.where(hap_0 == -1, -1, hap_1)

        # add delta's for phaseable variants
        self.delta[:, self.phaseable] = delta
        self.haplotypes[:, self.phaseable] = np.vstack([hap_0, hap_1])
        self.haplotypes = np.where(self.haplotypes[0, :] == self.haplotypes[1, :], -1, self.haplotypes).astype(int)

    def write_phased_vcf(self):
        """
        Write phased VCF. Only heterozygous variants that we were able to phase will be written with phase information.
        All other variants will be written without phase information.
        """
        # TODO add command line to header
        vcf = VCF(self.vcf_f)
        has_allele_level_info = vcf.contains('LV')  # has bubble info level
        vcf.add_format_to_header({'ID': "PS", "Number": 1, 'Type': "Integer", "Description": "Phase block ID"})
        # initialize phased VCF
        vcf_phased = Writer(self.output_vcf, vcf, 'wz')
        v_idx = 0
        block_id = 0
        block_ends = deque(self.block_ends)
        block_end = block_ends.popleft()
        for v in vcf(self.chrom):

            if (v.gt_types != 1 or not (v.is_snp or (v.is_indel and not self.snvs_only)) or
                    (len(v.ALT) > 1 and not self.multiallelics) or  # include multiallelics only if specified
                    v.genotypes[0][0] == -1 or
                    v.genotypes[0][1] == -1 or
                    np.max([len(v.REF), max([len(a) for a in v.ALT])]) > self.max_allele_length):
                gt = v.genotypes[0][:2]
                gt.append(False)
                v.genotypes = [gt]
                vcf_phased.write_record(v)
            else:
                if self.haplotypes[0, v_idx] == -1:
                    gt = self.idx_variant_mapping[v_idx]['gt'].copy()
                    gt.append(False)
                    v.genotypes = [gt]
                    vcf_phased.write_record(v)
                else:
                    # increment phase set ID
                    if v_idx > block_end:
                        block_id += 1
                        block_end = block_ends.popleft()
                    # assign phase set
                    v.set_format("PS", np.array([block_id]))
                    gt = [self.idx_variant_mapping[v_idx]['gt'][i] for i in self.haplotypes[:, v_idx]]
                    gt.append(True)
                    v.genotypes = [gt]
                    vcf_phased.write_record(v)
                v_idx += 1

    def write_haplotype_blocks(self):
        """
        Write haplotype block boundaries
        """
        start = 0
        haplotype_blocks = open(self.output_blocks, 'w')
        for end in self.block_ends:
            haplotype_blocks.write(f'{self.chrom}\t{self.idx_variant_mapping[start]["POS"]}\t'
                                   f'{self.idx_variant_mapping[end]["POS"]}\n')
            start = end + 1
        haplotype_blocks.close()

    def get_methylation_based_haplotag(self, read_name, meth_hap):
        """
        Get haplotype tag for reads based on methylation information that cannot be tagged otherwiser.
        :param read_name: str, read name
        :param meth_hap: str, indicating which methylation haplotype to use ('hap1' or 'hap2')
        :return: int, list, HP tag for read, variant idx
        """
        # get indices of read in list
        read_idx = [i for i, r in enumerate(self.methylation_read_assignments[meth_hap]) if r[0] == read_name]
        # get variant indices
        idx = np.array([self.methylation_read_assignments[meth_hap][r_idx][i][0] for i in range(1, 3) for r_idx in read_idx[:1]])
        # get haplotype states
        states = np.array([self.methylation_read_assignments[meth_hap][r_idx][i][1] for i in range(1, 3) for r_idx in read_idx[:1]])
        # compare to inferred haplotypes
        if np.all(self.haplotypes[0][idx] == states):
            return 1, idx
        elif np.all(self.haplotypes[1][idx] == states):
            return 2, idx
        else:
            return None, None

    def haplotag_reads(self):
        """
        Haplotag reads in bam file
        """
        # TODO add PG entry
        bam = pysam.AlignmentFile(self.output_bam, "wb", template=self.alignments)
        read_assignments = open(self.output_read_assignments, 'w')
        reads_0 = []
        reads_1 = []
        hp1 = 0
        hp2 = 0
        methylation_reads_hap1 = [r[0] for r in self.methylation_read_assignments['hap1']]
        methylation_reads_hap2 = [r[0] for r in self.methylation_read_assignments['hap2']]
        for read in self.alignments.fetch(self.chrom):
            read_name = read.query_name
            if read_name in self.read_states and len(self.read_states[read_name]) > 0:
                states = np.vstack([(k, v) for k, v in self.read_states[read_name].items()]).astype(int)
                idx = states[states[:, 1] != -1][:, 0]
                states = states[states[:, 1] != -1][:, 1]
                if states.shape[0] > 1:
                    # naive bayes calculator
                    # P(R | H) = P(H | R) x P(R)
                    # --> probability that read comes from haplotype
                    # X = probability of haplotype X given all reads times
                    # the probability of read matching inferred haplotype
                    prob_read = lambda hap: np.log(np.where(self.haplotypes[hap, idx] == states, 1 - self.error_rate,
                                                              self.error_rate))
                    prob_haplotype = lambda hap: (
                        np.concatenate([[self.delta[hap, idx[0]]],
                                        np.log(self.transition_matrix[self.haplotypes[hap, idx[:-1]],
                                                 self.haplotypes[hap, idx[1:]],
                                                 idx[:-1]])]))
                    # log(P(R) x P(H|R)) = log(P(R)) + log(P(H|R))
                    prob_0 = (prob_read(0) + prob_haplotype(0)).sum()
                    prob_1 = (prob_read(1) + prob_haplotype(1)).sum()
                    # calculate LLR  for coming from either haplotype
                    # log(P(R|H_A) / P(R|H_B)) = log(P(R|H_A)) - log(P(R|H_B))
                    if prob_0 - prob_1 >= self.llr_thresh:
                        reads_0.append(read_name)
                        read.set_tag("HP", 1)
                        if np.max(idx) > np.max(self.block_ends):
                            block = np.max(self.block_ends) + 1
                        else:
                            block = np.where(self.block_ends >= np.max(idx))[0][0]
                        read_assignments.write(f'{read_name}\tH1\t{block}\n')
                        hp1 += 1
                    elif prob_1 - prob_0 >= self.llr_thresh:
                        reads_1.append(read_name)
                        read.set_tag("HP", 2)
                        if np.max(idx) > np.max(self.block_ends):
                            block = np.max(self.block_ends) + 1
                        else:
                            block = np.where(self.block_ends >= np.max(idx))[0][0]
                        read_assignments.write(f'{read_name}\tH2\t{block}\n')
                        hp2 += 1
                    # check if methylation information allows tagging
                    elif read_name in methylation_reads_hap1:
                        tag, idx = self.get_methylation_based_haplotag(read_name, 'hap1')
                        if tag is None:
                            read_assignments.write(f'{read_name}\tnone\tnone\n')
                        elif tag == 1:
                            reads_0.append(read_name)
                            read.set_tag("HP", 1)
                            if np.max(idx) > np.max(self.block_ends):
                                block = np.max(self.block_ends) + 1
                            else:
                                block = np.where(self.block_ends >= np.max(idx))[0][0]
                            read_assignments.write(f'{read_name}\tH1\t{block}\n')
                            hp1 += 1
                        elif tag == 2:
                            reads_1.append(read_name)
                            read.set_tag("HP", 2)
                            if np.max(idx) > np.max(self.block_ends):
                                block = np.max(self.block_ends) + 1
                            else:
                                block = np.where(self.block_ends >= np.max(idx))[0][0]
                            read_assignments.write(f'{read_name}\tH2\t{block}\n')
                            hp2 += 1
                    elif read_name in methylation_reads_hap2:
                        tag, idx = self.get_methylation_based_haplotag(read_name, 'hap2')
                        if tag is None:
                            read_assignments.write(f'{read_name}\tnone\tnone\n')
                        elif tag == 1:
                            reads_0.append(read_name)
                            read.set_tag("HP", 1)
                            if np.max(idx) > np.max(self.block_ends):
                                block = np.max(self.block_ends) + 1
                            else:
                                block = np.where(self.block_ends >= np.max(idx))[0][0]
                            read_assignments.write(f'{read_name}\tH1\t{block}\n')
                            hp1 += 1
                        elif tag == 2:
                            reads_1.append(read_name)
                            read.set_tag("HP", 2)
                            if np.max(idx) > np.max(self.block_ends):
                                block = np.max(self.block_ends) + 1
                            else:
                                block = np.where(self.block_ends >= np.max(idx))[0][0]
                            read_assignments.write(f'{read_name}\tH2\t{block}\n')
                            hp2 += 1
                    else:
                        read_assignments.write(f'{read_name}\tnone\tnone\n')
                # check if methylation information allows tagging
                elif read_name in methylation_reads_hap1:
                    tag, idx = self.get_methylation_based_haplotag(read_name, 'hap1')
                    if tag is None:
                        read_assignments.write(f'{read_name}\tnone\tnone\n')
                    elif tag == 1:
                        reads_0.append(read_name)
                        read.set_tag("HP", 1)
                        if np.max(idx) > np.max(self.block_ends):
                            block = np.max(self.block_ends) + 1
                        else:
                            block = np.where(self.block_ends >= np.max(idx))[0][0]
                        read_assignments.write(f'{read_name}\tH1\t{block}\n')
                        hp1 += 1
                    elif tag == 2:
                        reads_1.append(read_name)
                        read.set_tag("HP", 2)
                        if np.max(idx) > np.max(self.block_ends):
                            block = np.max(self.block_ends) + 1
                        else:
                            block = np.where(self.block_ends >= np.max(idx))[0][0]
                        read_assignments.write(f'{read_name}\tH2\t{block}\n')
                        hp2 += 1
                elif read_name in methylation_reads_hap2:
                    tag, idx = self.get_methylation_based_haplotag(read_name, 'hap2')
                    if tag is None:
                        read_assignments.write(f'{read_name}\tnone\tnone\n')
                    elif tag == 1:
                        reads_0.append(read_name)
                        read.set_tag("HP", 1)
                        if np.max(idx) > np.max(self.block_ends):
                            block = np.max(self.block_ends) + 1
                        else:
                            block = np.where(self.block_ends >= np.max(idx))[0][0]
                        read_assignments.write(f'{read_name}\tH1\t{block}\n')
                        hp1 += 1
                    elif tag == 2:
                        reads_1.append(read_name)
                        read.set_tag("HP", 2)
                        if np.max(idx) > np.max(self.block_ends):
                            block = np.max(self.block_ends) + 1
                        else:
                            block = np.where(self.block_ends >= np.max(idx))[0][0]
                        read_assignments.write(f'{read_name}\tH2\t{block}\n')
                        hp2 += 1
                else:
                    read_assignments.write(f'{read_name}\tnone\tnone\n')
            # check if methylation information allows tagging
            elif read_name in methylation_reads_hap1:
                tag, idx = self.get_methylation_based_haplotag(read_name, 'hap1')
                if tag is None:
                    read_assignments.write(f'{read_name}\tnone\tnone\n')
                elif tag == 1:
                    reads_0.append(read_name)
                    read.set_tag("HP", 1)
                    if np.max(idx) > np.max(self.block_ends):
                        block = np.max(self.block_ends) + 1
                    else:
                        block = np.where(self.block_ends >= np.max(idx))[0][0]
                    read_assignments.write(f'{read_name}\tH1\t{block}\n')
                    hp1 += 1
                elif tag == 2:
                    reads_1.append(read_name)
                    read.set_tag("HP", 2)
                    if np.max(idx) > np.max(self.block_ends):
                        block = np.max(self.block_ends) + 1
                    else:
                        block = np.where(self.block_ends >= np.max(idx))[0][0]
                    read_assignments.write(f'{read_name}\tH2\t{block}\n')
                    hp2 += 1
            elif read_name in methylation_reads_hap2:
                tag, idx = self.get_methylation_based_haplotag(read_name, 'hap2')
                if tag is None:
                    read_assignments.write(f'{read_name}\tnone\tnone\n')
                elif tag == 1:
                    reads_0.append(read_name)
                    read.set_tag("HP", 1)
                    if np.max(idx) > np.max(self.block_ends):
                        block = np.max(self.block_ends) + 1
                    else:
                        block = np.where(self.block_ends >= np.max(idx))[0][0]
                    read_assignments.write(f'{read_name}\tH1\t{block}\n')
                    hp1 += 1
                elif tag == 2:
                    reads_1.append(read_name)
                    read.set_tag("HP", 2)
                    if np.max(idx) > np.max(self.block_ends):
                        block = np.max(self.block_ends) + 1
                    else:
                        block = np.where(self.block_ends >= np.max(idx))[0][0]
                    read_assignments.write(f'{read_name}\tH2\t{block}\n')
                    hp2 += 1
            else:
                read_assignments.write(f'{read_name}\tnone\tnone\n')
            bam.write(read)
        bam.close()
        read_assignments.close()


def read_phasing(args):
    if args.verbose:
        logging.basicConfig(handlers=[logging.StreamHandler(), logging.FileHandler(args.log, mode='w')],
                            level=logging.INFO,
                            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(handlers=[logging.FileHandler(args.log, mode='w')],
                            level=logging.INFO,
                            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
    if args.ont and not args.pacbio:
        seqtech = 'ont'
        flank_snv = 66
        flank_indel = 200
    elif args.pacbio and not args.ont:
        seqtech = 'pacbio'
        flank_snv = 33
        flank_indel = 100
    else:
        raise ValueError('Please specify sequencing technology either with --ont or --pacbio flag. '
                         'Only supporting PacBio or ONT.')
    if args.flank_snv is not None:
        flank_snv = args.flank_snv
    if args.flank_indel is not None:
        flank_indel = args.flank_indel
    longhap = LongHap(vcf_f=args.vcf, bam=args.bam, chrom=args.chrom, reference_path=args.reference,
                      output_vcf=args.output_vcf, output_read_states=args.output_read_states,
                      output_blocks=args.output_blocks, output_bam=args.output_bam,
                      output_variant_read_mapping=args.output_variant_read_mapping,
                      output_read_assignments=args.output_read_assignments, methylation_calls_f=args.methylation_calls,
                      output_transition_matrix=args.output_transition_matrix,
                      output_allele_coverage=args.output_allele_coverage,
                      output_transition_matrix_meth=args.output_transition_matrix_meth,
                      output_differentially_methylated_sites=args.output_differentially_methylated_sites,
                      output_unphaseable_variants=args.output_unphaseable_variants,
                      snvs_only=args.snvs_only, multiallelics=args.multiallelics,
                      use_all_methylated_sites=args.use_all_methylated_sites, force=args.force,
                      max_allele_length=args.max_allele_length, min_allele_count=args.min_allele_count,
                      min_base_quality=args.min_base_quality, seqtech=seqtech, min_mapq=args.min_mapq,
                      flank_snv=flank_snv, flank_indel=flank_indel,
                      )
    longhap.infer_variant_transitions()
    longhap.infer_methylation_transitions()
    longhap.phase()
    longhap.write_results()


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version=__version__)

    parser.add_argument('--vcf', help='Input VCF with called variants', required=True)
    parser.add_argument('-b', '--bam', help='Sorted alignment bam', required=True)
    parser.add_argument('-r', '--reference', help='Reference fasta. Must be indexed with samtools faidx',
                        required=True)
    parser.add_argument('-c', '--chrom', help='Chromosome', required=True)
    parser.add_argument('-m', '--methylation_calls', help='Methylation calls from pileup model',
                           required=False, default=None)
    parser.add_argument('--snvs_only', help='Whether to phase SNVs only ["False]',
                           default=False, action='store_true')
    parser.add_argument('--multiallelics', help='Also phase multiallelic variants or not [False]',
                           default=False, action='store_true')
    parser.add_argument('--ont', help='Data is Oxford Nanopore data [False]', default=False, action='store_true')
    parser.add_argument('--pacbio', help='Data is PacBio HiFi data [False]', default=False, action='store_true')
    parser.add_argument('--max_allele_length', default=50000, type=int,
                           help='Maximum length of alleles to consider for phasing in bp [50000]',)
    parser.add_argument('--min_allele_count',
                           help='How many examples of the minor allele must be present in the reads to consider the '
                                'variant for phasing [1]', type=int, default=1)
    parser.add_argument('--min_base_quality',
                           help='Minimum base quality to consider a base for phasing. Only affects SNP phasing. '
                                'For HiFi data, all bases should be consider, that is a minimum quality of 0. '
                                'For ONT data, a threshold of 10 is recommended [0]', type=int, default=0)
    parser.add_argument('--min_mapq', help='Minimum mapping quality to consider a read for phasing [20]',
                        type=int, default=20)
    parser.add_argument('--flank_snv', help='Number of flanking bp to use for realignment around '
                                            'uncertain SNVs. Default is 66 for ONT and 33 for PacBio', type=int,
                        default=None)
    parser.add_argument('--flank_indel', help='Number of flanking bp to use for realignment around '
                                              'uncertain indels. Default is 200 for ONT and 100 for PacBio', type=int,
                        default=None)
    parser.add_argument('-o', '--output_vcf', help='Output phased vcf', required=True)
    parser.add_argument('--output_bam', help='Output haplotagged bam')
    parser.add_argument('--output_read_assignments', help='Haplotype assignments for each read')
    parser.add_argument('--output_blocks', help='Haplotype blocks in bed format')
    parser.add_argument('--output_transition_matrix', required=False, default=None,
                           help='If provided transition matrix will be saved to this file as numpy array (.npz). '
                                'Allows faster re-runs.')
    parser.add_argument('--output_transition_matrix_meth', required=False, default=None,
                           help='If provided transition matrix filled in with methylation data will be saved to this '
                                'file as numpy array (.npz). Allows faster re-runs.')
    parser.add_argument('--output_read_states', required=False, default=None,
                           help='If provided read states will be saved to this file as json. '
                                'Allows faster re-runs.')
    parser.add_argument('--output_variant_read_mapping', required=False, default=None,
                           help='If provided read names covering a specific variant will be saved to this file as json.'
                                ' Allows faster re-runs. ')
    parser.add_argument('--output_allele_coverage', required=False, default=None,
                           help='If provided allele coverage will be saved to this file as npy (.npy). '
                                'Sites with one allele absent from reads bill be ignored. Allows faster re-runs.')
    parser.add_argument('--output_unphaseable_variants', required=False, default=None,
                           help='If provided unphaseable variants will be saved to this file as npz. '
                                'Allows faster re-runs.')
    parser.add_argument('--output_differentially_methylated_sites',
                           help='Write differentially methylated files used by longhap to infer transitions to file',
                           required=False, default=None)
    parser.add_argument('--use_all_methylated_sites', default=False, action='store_true',
                           help="Whether to use all methylated sites or not. If False, at most 25,000 methylated sites "
                                "per transition are used. This guarantees fast runtimes and does not seem to "
                                "sacrifice accuracy. [False]")

    parser.add_argument('--force', action='store_true', default=False,
                           help='If transition matrix output is provided and file already exists this file will be '
                                'loaded by default unless --force is set. Then the transition matrix will '
                                'be re-inferred.')
    parser.add_argument('--log', help='Log file', default='longhap.log')
    parser.add_argument('-v', '--verbose', help='Print logging information to stdout', action='store_true',
                           default=False)

    args = parser.parse_args(argv)
    read_phasing(args)


if __name__ == "__main__":
    main(sys.argv[1:])
