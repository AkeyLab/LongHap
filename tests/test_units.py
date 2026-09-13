"""Unit tests for the pure and near-pure helpers.

These need no files: where a method touches ``self``, the test builds a bare
LongHap with ``object.__new__`` and sets only the attributes that method reads.
"""
from collections import deque

import numpy as np
import pysam
import pytest

from longhap import LongHap


@pytest.fixture
def bare():
    """A LongHap with no I/O done; set attributes per test."""
    def _bare(**attrs):
        lh = object.__new__(LongHap)
        for k, v in attrs.items():
            setattr(lh, k, v)
        return lh
    return _bare


# --------------------------------------------------------------------------- #
# CIGAR walking, checked against pysam
# --------------------------------------------------------------------------- #
def aligned_pairs_oracle(cigar, read_start):
    """{reference position: query position} from pysam, for M/=/X columns only."""
    a = pysam.AlignedSegment()
    a.query_name = 'r'
    a.reference_id = 0
    a.reference_start = read_start
    a.cigartuples = cigar
    n_query = sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
    a.query_sequence = 'A' * n_query
    return {r: q for q, r in a.get_aligned_pairs() if q is not None and r is not None}


CIGARS = [
    ([(0, 100)], 1000),
    ([(4, 20), (0, 100)], 1000),
    ([(0, 40), (1, 5), (0, 60)], 1000),
    ([(0, 40), (2, 7), (0, 60)], 1000),
    ([(4, 15), (0, 30), (2, 4), (0, 20), (1, 6), (0, 40), (4, 10)], 1000),
    ([(0, 25), (8, 3), (0, 25), (2, 2), (7, 30)], 1000),
    ([(0, 50), (1, 2), (0, 50), (4, 8)], 1000),
    ([(5, 10), (0, 60), (5, 10)], 1000),
]


@pytest.mark.parametrize('cigar,read_start', CIGARS)
def test_cigar_walker_matches_pysam(cigar, read_start):
    """Every aligned reference column must map to the query index pysam reports.

    The walker is stateful, so positions are fed in ascending order exactly as
    longhap feeds them while scanning a read's variants.
    """
    oracle = aligned_pairs_oracle(cigar, read_start)
    positions = sorted(oracle)[::7]     # sample, keeps the test readable

    walker = deque(cigar)
    ref_offset = query_offset = 0
    operation, length = None, 0
    for pos in positions:
        query_offset, ref_offset, operation, length = (
            LongHap.get_query_and_reference_indices_for_variant_cigar(
                walker, pos, read_start, ref_offset, query_offset, operation, length))
        assert query_offset == oracle[pos], (
            f'cigar={cigar} pos={pos}: walker said q={query_offset}, '
            f'pysam says q={oracle[pos]}')
        assert ref_offset == pos - read_start, (
            f'cigar={cigar} pos={pos}: ref_offset={ref_offset}, '
            f'expected {pos - read_start}')


def walk_all(cigar, read_start):
    """Walk every aligned reference column; return the positions pysam disagrees with."""
    o = aligned_pairs_oracle(cigar, read_start)
    walker = deque(cigar)
    ref_offset = query_offset = 0
    operation, length = None, 0
    wrong = []
    for pos in sorted(o):
        query_offset, ref_offset, operation, length = (
            LongHap.get_query_and_reference_indices_for_variant_cigar(
                walker, pos, read_start, ref_offset, query_offset, operation, length))
        if query_offset != o[pos]:
            wrong.append((pos, query_offset, o[pos]))
    return wrong


def test_cigar_walker_handles_a_hard_clip_before_an_insertion():
    """Hard clips ride on supplementary alignments, which longhap does read."""
    wrong = walk_all([(5, 10), (0, 20), (1, 3), (0, 20)], 100)
    assert not wrong, (
        f'{len(wrong)} positions mis-mapped, first three: {wrong[:3]} '
        '(pos, walker query index, pysam query index)')


def test_cigar_walker_handles_a_hard_clip_without_an_indel():
    """A hard clip alone is harmless: query and reference stay collinear."""
    assert not walk_all([(5, 10), (0, 60), (5, 10)], 100)


def test_cigar_walker_handles_a_reference_skip():
    wrong = walk_all([(0, 20), (3, 10), (0, 20)], 100)
    assert not wrong, f'{len(wrong)} positions mis-mapped, first three: {wrong[:3]}'


def test_cigar_walker_is_consistent_when_called_once_per_position():
    """Restarting the walker for each position must give the same answer."""
    cigar, read_start = CIGARS[4]
    oracle = aligned_pairs_oracle(cigar, read_start)
    for pos in sorted(oracle)[::11]:
        q, r, _, _ = LongHap.get_query_and_reference_indices_for_variant_cigar(
            deque(cigar), pos, read_start, 0, 0, None, 0)
        assert q == oracle[pos], f'pos={pos}: got q={q}, expected {oracle[pos]}'


# --------------------------------------------------------------------------- #
# mirror_transition
# --------------------------------------------------------------------------- #
def test_mirror_transition_copies_the_confident_row():
    t = np.array([[9.0, 1.0], [1.0, 20.0]])
    out = LongHap.mirror_transition(t.copy(), normalized=False)
    assert out[0].tolist() == [20.0, 1.0], out
    assert out[1].tolist() == [1.0, 20.0], out


def test_mirror_transition_collapses_a_tie():
    t = np.array([[5.0, 5.0], [5.0, 5.0]])
    out = LongHap.mirror_transition(t.copy(), normalized=True)
    assert (out == 0.5).all(), out


def test_mirror_transition_leaves_an_agreeing_matrix_symmetric():
    t = np.array([[10.0, 1.0], [1.0, 10.0]])
    out = LongHap.mirror_transition(t.copy(), normalized=False)
    assert out[0].tolist() == out[1][::-1].tolist(), out


def test_mirror_transition_does_not_mutate_its_argument():
    matrix = np.zeros((2, 2, 3)) + 1e-20
    matrix[:, :, 0] = np.array([[9.0, 1.0], [1.0, 20.0]])
    before = matrix[:, :, 0].copy()
    LongHap.mirror_transition(matrix[:, :, 0], normalized=False)
    assert np.array_equal(matrix[:, :, 0], before), (
        f'transition_matrix changed underneath the caller: '
        f'{before.tolist()} -> {matrix[:, :, 0].tolist()}')


# --------------------------------------------------------------------------- #
# transition counting
# --------------------------------------------------------------------------- #
def test_allele_transitions_counted_from_read_states(bare):
    lh = bare(
        variant_read_mapping={'0': ['a', 'b', 'c', 'd'], '1': ['a', 'b', 'c', 'e']},
        read_states={'a': {'0': 0, '1': 0}, 'b': {'0': 0, '1': 0},
                     'c': {'0': 1, '1': 1}, 'd': {'0': 1}, 'e': {'1': 0}},
    )
    t = lh.get_allele_transitions_from_known_read_states(0, 1)
    assert t[0, 0] == pytest.approx(2)
    assert t[1, 1] == pytest.approx(1)
    assert t[0, 1] < 1e-10 and t[1, 0] < 1e-10


def test_allele_transitions_skip_undetermined_states(bare):
    lh = bare(
        variant_read_mapping={'0': ['a', 'b'], '1': ['a', 'b']},
        read_states={'a': {'0': -1, '1': 0}, 'b': {'0': 1, '1': -1}},
    )
    t = lh.get_allele_transitions_from_known_read_states(0, 1)
    assert t.sum() < 1e-10, t


# --------------------------------------------------------------------------- #
# uncertain transitions
# --------------------------------------------------------------------------- #
def test_uncertain_transitions_finds_flat_columns():
    tm = np.zeros((2, 2, 4))
    tm[:, :, 0] = [[0.9, 0.1], [0.1, 0.9]]
    tm[:, :, 1] = 0.5
    tm[:, :, 2] = [[0.2, 0.8], [0.8, 0.2]]
    tm[:, :, 3] = 0.5
    assert list(LongHap.get_uncertain_transitions(tm)) == [1, 3]


def test_uncertain_transitions_tolerates_float_noise():
    """A column a hair off 0.5 is still recognised as uninformative."""
    tm = np.full((2, 2, 1), 0.5)
    tm[0, 0, 0] = 0.5 + 1e-12
    tm[0, 1, 0] = 0.5 - 1e-12
    assert list(LongHap.get_uncertain_transitions(tm)) == [0]


# --------------------------------------------------------------------------- #
# three-variant re-estimation
# --------------------------------------------------------------------------- #
def _three_variant_lh(bare, t_ab, t_bc):
    """A LongHap positioned to re-estimate the A-B and B-C transitions.

    The method bails out when *both* neighbouring transitions already point
    somewhere (their argmax rows disagree), so at least one of them is given a
    flat, uninformative shape here — which is the situation it exists to fix.
    """
    states = {f'r{i}': {'0': i % 2, '1': i % 2, '2': i % 2} for i in range(10)}
    tm = np.zeros((2, 2, 4)) + 1e-20
    tm[:, :, 0] = t_ab
    tm[:, :, 1] = t_bc
    return bare(
        transition_matrix=tm,
        num_variants=5,
        variant_read_mapping={k: list(states) for k in ('0', '1', '2')},
        read_states=states,
    )


def test_three_variant_update_bails_out_when_both_edges_are_informative(bare):
    lh = _three_variant_lh(bare, [[0.8, 0.2], [0.2, 0.8]], [[0.7, 0.3], [0.3, 0.7]])
    _, _, new_t1, new_t2, _ = \
        lh.update_transition_matrix_considering_adjacent_variants(0, 1, 2)
    assert new_t1.sum() < 1e-15 and new_t2.sum() < 1e-15


def test_three_variant_update_row_zero_is_a_distribution(bare):
    """new_t1's REF row is a complementary pair, so it must sum to 1."""
    lh = _three_variant_lh(bare, [[0.8, 0.2], [0.8, 0.2]], [[0.7, 0.3], [0.3, 0.7]])
    _, _, new_t1, _, _ = lh.update_transition_matrix_considering_adjacent_variants(0, 1, 2)
    assert new_t1[0].sum() == pytest.approx(1.0, abs=1e-9), new_t1


def test_three_variant_update_row_one_is_a_distribution(bare):
    lh = _three_variant_lh(bare, [[0.8, 0.2], [0.8, 0.2]], [[0.7, 0.3], [0.3, 0.7]])
    _, _, new_t1, _, _ = lh.update_transition_matrix_considering_adjacent_variants(0, 1, 2)
    assert new_t1[1].sum() == pytest.approx(1.0, abs=1e-9), (
        f'row 1 sums to {new_t1[1].sum()}, row 0 sums to {new_t1[0].sum()}')


def test_three_variant_update_new_t2_rows_are_distributions(bare):
    """new_t2 uses the same construction and is spelled correctly."""
    lh = _three_variant_lh(bare, [[0.8, 0.2], [0.8, 0.2]], [[0.7, 0.3], [0.3, 0.7]])
    _, _, _, new_t2, _ = lh.update_transition_matrix_considering_adjacent_variants(0, 1, 2)
    assert new_t2[0].sum() == pytest.approx(1.0, abs=1e-9), new_t2
    assert new_t2[1].sum() == pytest.approx(1.0, abs=1e-9), new_t2


def test_three_variant_update_ratio_survives_renormalisation(bare):
    """Both call sites renormalise, which happens to hide ISSUE-11.

    Row 1's two entries share the wrong denominator, so their *ratio* is still
    right and a row-wise renormalisation recovers the intended value.  This test
    documents why the typo has not shown up as a phasing error, and will keep
    holding after it is fixed.
    """
    lh = _three_variant_lh(bare, [[0.8, 0.2], [0.8, 0.2]], [[0.7, 0.3], [0.3, 0.7]])
    _, _, new_t1, _, _ = lh.update_transition_matrix_considering_adjacent_variants(0, 1, 2)
    renormalised = new_t1 / new_t1.sum(axis=1, keepdims=True)
    assert renormalised.sum(axis=1) == pytest.approx([1.0, 1.0])
    assert 0.0 <= renormalised.min() and renormalised.max() <= 1.0


# --------------------------------------------------------------------------- #
# Viterbi backtrace
# --------------------------------------------------------------------------- #
def test_backtrace_follows_the_pointers():
    delta = np.array([[-1.0, -2.0, -3.0],
                      [-5.0, -0.5, -0.2]])
    phi = np.array([[0, 1, 1],
                    [1, 0, 0]])
    hap = np.zeros(3, dtype=int) - 1
    out = LongHap.backtrace(delta, hap, phi, 2, 0)
    assert out[2] == 1                 # argmax of delta[:, 2]
    assert out[1] == phi[1, 2]
    assert out[0] == phi[out[1], 1]


def test_backtrace_leaves_positions_outside_the_block_untouched():
    delta = np.zeros((2, 5))
    delta[1, 4] = 1.0
    phi = np.zeros((2, 5), dtype=int)
    hap = np.zeros(5, dtype=int) - 1
    out = LongHap.backtrace(delta, hap, phi, 4, 3)
    assert out[0] == out[1] == out[2] == -1, out


# --------------------------------------------------------------------------- #
# methylation-derived transitions
# --------------------------------------------------------------------------- #
def test_methylation_transition_helper_is_a_distribution():
    hap = np.array([0.95, 0.93, 0.9])
    t = LongHap.calculate_transition_probability_from_methylation_helper(hap, 0)
    assert t.sum(axis=1) == pytest.approx([1.0, 1.0], abs=1e-6)


def test_methylation_transition_helper_prefers_staying_when_states_agree():
    """Both sites read as haplotype 1, so 1->1 should beat 1->0."""
    hap = np.array([0.99, 0.99])
    t = LongHap.calculate_transition_probability_from_methylation_helper(hap, 0)
    t = t / t.sum(axis=1, keepdims=True)
    assert t[1, 1] > t[1, 0], t


def test_methylation_transition_helper_prefers_switching_when_states_differ():
    hap = np.array([0.01, 0.99])
    t = LongHap.calculate_transition_probability_from_methylation_helper(hap, 0)
    t = t / t.sum(axis=1, keepdims=True)
    assert t[0, 1] > t[0, 0], t


# --------------------------------------------------------------------------- #
# misc helpers
# --------------------------------------------------------------------------- #
def test_methylation_haplotag_declines_an_unknown_read(bare):
    lh = bare(methylation_read_assignments={'hap1': []},
              haplotypes=np.array([[0, 1, 0], [1, 0, 1]]))
    tag, idx = lh.get_methylation_based_haplotag('missing-read', 'hap1')
    assert tag is None, f'unknown read was assigned to haplotype {tag}'
