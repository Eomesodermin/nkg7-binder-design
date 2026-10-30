"""Tests for the sequence-complexity metrics."""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipeline"))

import sequence_complexity as sc


def test_ala_fraction():
    assert sc.ala_fraction("AAAA") == 1.0
    assert sc.ala_fraction("ACDE") == 0.25


def test_entropy_bounds():
    assert sc.shannon_entropy("AAAAAAAA") == 0.0
    assert math.isclose(sc.shannon_entropy("ACDEFGHIKLMNPQRSTVWY"), math.log2(20))


def test_longest_homopolymer():
    assert sc.longest_homopolymer("ACDEFG") == 1
    assert sc.longest_homopolymer("ACCCCDE") == 4


def test_low_complexity_flags_poly_alanine():
    assert sc.is_low_complexity("DAAAAAAAARAAALRATLEGHRAQGEANVARADALGPELGELARAALEAGLARLEAKLA")


def test_diverse_sequence_is_not_flagged():
    assert not sc.is_low_complexity("MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATP")


def test_summarise_counts_the_funnel():
    rows = [
        {"arm": "a", "design": "d1", "iptm": "0.9", "complex_plddt": "0.8", "pass": "True",
         "binder_seq": "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATP"},
        {"arm": "a", "design": "d2", "iptm": "0.9", "complex_plddt": "0.8", "pass": "True",
         "binder_seq": "A" * 40},
    ]
    funnel = dict(sc.summarise(sc.annotate(rows)))
    assert funnel["designs scored"] == 2
    assert funnel["ipTM >= 0.5"] == 2
    assert funnel["not low-complexity"] == 1
