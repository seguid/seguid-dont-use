#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

from seguid.asserts import assert_anneal
from seguid.asserts import assert_in_alphabet
from seguid.asserts import assert_table

from seguid.tables import COMPLEMENT_TABLE_DNA
from seguid.tables import COMPLEMENT_TABLE_RNA
from seguid.tables import COMPLEMENT_TABLE_IUPAC
from seguid.tables import TABLE_IUPAC_PROTEIN

def test_assert_in_alphabet():
    """docstring."""
    seq = "ABCDEFGH"
    alphabet = {'A', 'C', 'G', 'T'}

    assert_in_alphabet("ACGT", alphabet = alphabet)
    assert_in_alphabet("AAAA", alphabet = alphabet)
    assert_in_alphabet("", alphabet = alphabet)

    try:
        assert_in_alphabet("x", alphabet = alphabet)
        print("Should not be reached")
    except ValueError:
        pass

    try:
        assert_in_alphabet("ACGTx", alphabet = alphabet)
        print("Should not be reached")
    except ValueError:
        pass


def test_assert_table():
    """docstring."""
    assert_table(str.maketrans("GATC", "CTAG"))
    assert_table(COMPLEMENT_TABLE_DNA)

    try:
        assert_table(str.maketrans("GATx", "CTAG"))
        print("Should not be reached")
    except ValueError:
        pass


if __name__ == "__main__":
    pytest.main([__file__, "-vvv", "-s"])
