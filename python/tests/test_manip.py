#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

from seguid.manip import rc
from seguid.manip import complementary
from seguid.manip import rotate
from seguid.manip import rotate_to_min
from seguid.config import set_min_rotation
from seguid.config import _min_rotation


def test_sort_order():

    from string import printable

    assert "".join(sorted(printable)) == ('\t\n\x0b\x0c\r !"'
                                          "#$%&\'"
                                          "()*+,-./0123456789:;<=>?@"
                                          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                          "[\\]^_`"
                                          "abcdefghijklmnopqrstuvwxyz"
                                          "{|}~")

def test_complementary():
    """docstring."""
    seq = "AACCGGTT"
    seq_complementary = "TTGGCCAA"
    assert complementary(seq) == seq_complementary
    assert complementary(seq_complementary) == seq
    assert complementary(complementary(seq)) == seq

    seq = "AACCGGTTxx"
    try:
        print(complementary(seq))
        print("Should not be reached")
    except ValueError:
        pass


def test_rc2():
    """docstring."""
    watson = "ACGTAACCGGTT"
    crick = "AACCGGTTACGT"
    assert rc(watson) == crick
    assert rc(crick) == watson
    assert rc(rc(watson)) == watson


def test_rotate():
    """docstring."""
    seq = "ACGTAACCGGTT"
    n = len(seq)
    assert rotate(seq, amount=0) ==  (seq, "")
    assert rotate(seq, amount=n) == rotate(seq, amount=0)
    assert rotate(seq, amount=2*n) == rotate(seq, amount=0)
    assert rotate(seq, amount=n-1) == rotate(seq, amount=-1)

    assert rotate("", amount=0) == ("", "")
    assert rotate("", amount=1) == ("", "")

    ## Rotate on the complementary strand
    assert complementary(rotate(complementary(seq), amount=+1)[0]) == rotate(seq, amount=+1)[0]

    from seguid.reprutils import tuple_from_repr, repr_from_tuple

    watson = "GATACCA"
    crick =  "CTATGGT"[::-1]

    n = len(watson)
    assert rotate(watson, crick, amount=0) ==  (watson, crick)
    assert rotate(watson, crick, amount=n) == rotate(watson, crick, amount=0)
    assert rotate(watson, crick, amount=2*n) == rotate(watson, crick, amount=0)
    assert rotate(watson, crick, amount=n-1) == rotate(watson, crick, amount=-1)

    assert rotate("","", amount=0) == ("", "")
    assert rotate("", "", amount=1) == ("", "")

    # https://onlinestringtools.com/rotate-string


def test_rotate2():
    """docstring."""
    watson = "ACGTAACCGGTT"
    crick =  "TGCATTGGCCAA"[::-1]

    ## Rotate on Watson, is the opposite rotation on Crick
    assert rc(watson)                 == crick
    assert rc(rotate(crick, amount=-1)[0])      == rotate(watson, amount=+1)[0]
    assert rc(rotate(rc(watson), amount=-1)[0]) == rotate(watson, amount=+1)[0]


def test_rc():
    """docstring."""
    assert rc("GAT") == "ATC"
    assert rc("GTT") == "AAC"

    with pytest.raises(ValueError):
        rc("GTZ")


def test_min_rotation_pydivsufsort():
    """Tests for the pydivsufsort min_rotation"""
    set_min_rotation("pydivsufsort")
    assert _min_rotation("") == 0
    assert _min_rotation("Aa") == 0
    assert rotate_to_min("taaa") == ("aaat", "")
    assert (
        rotate_to_min("abaabaaabaababaaabaaababaab")
        == ("aaabaaababaababaabaaabaabab", "")
    )
    assert (
        rotate_to_min("abaabaaabaababaaabaaaBabaab")
        == ("Babaababaabaaabaababaaabaaa", "")
    )
    set_min_rotation("built-in")


def test_min_rotation_built_in():
    set_min_rotation("built-in")
    assert _min_rotation("Aa") == 0
    assert rotate_to_min("TAAA") == ("AAAT", "")
    assert (
        rotate_to_min("ACAACAAACAACACAAACAAACACAAC")
        == ("AAACAAACACAACACAACAAACAACAC", "")
    )
    assert _min_rotation("") == 0

    watson = "GATACCA"
    crick =  "CTATGGT"[::-1]

    assert rotate_to_min(watson, "") == ("ACCAGAT", "")
    assert rotate_to_min(watson, crick) == ("ACCAGAT", "ATCTGGT")
