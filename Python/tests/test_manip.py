#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

from seguid.manip import rc
from seguid.manip import min_rotation_py
from seguid.manip import complementary
from seguid.manip import rotate

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
    assert rotate(seq,   0) ==        seq
    assert rotate(seq,   n) == rotate(seq,  0)
    assert rotate(seq, 2*n) == rotate(seq,  0)
    assert rotate(seq, n-1) == rotate(seq, -1)

    ## Rotate on the complementary strand
    assert complementary(rotate(complementary(seq), +1)) == rotate(seq, +1)


def test_rotate2():
    """docstring."""
    watson = "ACGTAACCGGTT"
    crick = "AACCGGTTACGT"

    ## Rotate on Watson, is the opposite rotation on Crick
    assert rc(watson)                 == crick
    assert rc(rotate(crick, -1))      == rotate(watson, +1)
    assert rc(rotate(rc(watson), -1)) == rotate(watson, +1)


def test_rc():
    """docstring."""
    assert rc("GAT") == "ATC"
    assert rc("GTT") == "AAC"

    with pytest.raises(ValueError):
        rc("GTZ")

        
if __name__ == "__main__":
    pytest.main([__file__, "-vvv", "-s"])
