#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest

from pathlib import Path
from hashlib import sha1
from base64 import urlsafe_b64encode as b64us

from seguid import slseguid
from seguid import scseguid
from seguid import dlseguid
from seguid import dcseguid

from seguid.reprutils import tuple_from_repr
from seguid.reprutils import repr_from_tuple

from seguid.tables import COMPLEMENT_TABLE_DNA
from seguid.tables import TABLE_IUPAC_PROTEIN

from seguid.asserts import assert_anneal
from seguid.asserts import assert_in_alphabet
from seguid.manip import rc
from seguid.manip import min_rotation_py
from seguid.manip import complementary
from seguid.manip import rotate

from seguid.chksum import seguid


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


def test_min_rotation():
    """Tests for the pydivsufsort min_rotation"""
    from pydivsufsort import min_rotation

    def smallest_rotation(s):
        i = min_rotation(s)
        return s[i:] + s[:i]

    assert smallest_rotation("taaa") == "aaat"
    assert (
        smallest_rotation("abaabaaabaababaaabaaababaab")
        == "aaabaaababaababaabaaabaabab"
    )
    assert (
        smallest_rotation("abaabaaabaababaaabaaaBabaab")
        == "Babaababaabaaabaababaaabaaa"
    )


def test_min_rotation_py():
    """docstring."""

    def smallest_rotation_py(s):
        i = min_rotation_py(s)
        return s[i:] + s[:i]

    assert smallest_rotation_py("TAAA") == "AAAT"
    assert (
        smallest_rotation_py("ACAACAAACAACACAAACAAACACAAC")
        == "AAACAAACACAACACAACAAACAACAC"
    )


def test_seguid():
    assert seguid("AT") == "seguid:Ax/RG6hzSrMEEWoCO1IWMGska+4"
    NP_313053_1 = (
        "MKALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSG"
        "ASRGIRLLQEEEEGLPLVGRVAAGEPLLAQQHIEGHYQVDPSLFKPNADFLLRVSGMSMKD"
        "IGIMDGDLLAVHKTQDVRNGQVVVARIDDEVTVKRLKKQGNKVELLPENSEFKPIVVDLRQ"
        "QSFTIEGLAVGVIRNGDWL"
    )

    assert seguid(NP_313053_1, table=TABLE_IUPAC_PROTEIN) == "seguid:2c4yjE+JqjvzYF1d0OmUh8pCpz8"




def cs(arg):
    return (
        b64us(sha1(arg.encode("ASCII")).digest()).decode("ASCII").rstrip("=")
    )


def test_slseguid():
    assert slseguid("AT") == "slseguid:Ax_RG6hzSrMEEWoCO1IWMGska-4"
    assert cs("AT") in slseguid("AT")


def test_scseguid():
    m13dna = Path("test_data/M13.txt").read_text().strip()
    sc = "scseguid:aAjgnsF9cPI6cu8IQ81sYnstVzU"
    assert scseguid(m13dna) == sc
    assert cs(Path("test_data/M13msg.txt").read_text().strip()) in sc


def test_dlseguid():
    ct = COMPLEMENT_TABLE_DNA
    table = ct | {"\n":"\n", "-":"-"}
    # AT
    # TA

    dlDNA = "AT"
    dlDNA_dlseguid = "AWD-dt5-TEua8RbOWfnctJIu9nA"
    assert dlseguid(dlDNA, rc(dlDNA), 0) == f"dlseguid:{dlDNA_dlseguid}"
    assert dlDNA_dlseguid in slseguid("AT\nTA", table = table)
    assert cs("AT\nTA") == dlDNA_dlseguid

    #  -AT
    #  AT-

    dlDNA2 = ("AT", "TA", 1)
    dlDNA2_dlseguid = "JwB2eUmZkCNjyWAv471JeUbiSDM"
    assert dlseguid(*dlDNA2) == f"dlseguid:{dlDNA2_dlseguid}"
    assert dlDNA2_dlseguid in slseguid("-AT\nAT-", table = table)
    assert cs("-AT\nAT-") == dlDNA2_dlseguid

    # TA-
    # -TA

    dlDNA3 = ("TA", "AT", -1)
    dlDNA3_dlseguid = "bv0UOR12eWrBeaAx79PNZvveviU"
    assert dlseguid(*dlDNA3) == f"dlseguid:{dlDNA3_dlseguid}"
    assert dlDNA3_dlseguid in slseguid("AT-\n-AT", table = table)
    assert cs("AT-\n-AT") == dlDNA3_dlseguid

    # CTATAG
    #   TA

    dlDNA4 = ("CTATAG", "AT", -2)
    dlDNA4_dlseguid = "np3hncfQvOh8rZ8Co1Ts_02NXg4"
    assert dlseguid(*dlDNA4) == f"dlseguid:{dlDNA4_dlseguid}"
    assert dlDNA4_dlseguid in slseguid("--AT--\nGATATC", table = table)
    assert cs("--AT--\nGATATC") == dlDNA4_dlseguid

    #   AT
    # GATATC

    dlDNA5 = ("AT", "CTATAG", 2)
    dlDNA5_dlseguid = "np3hncfQvOh8rZ8Co1Ts_02NXg4"
    assert dlseguid(*dlDNA5) == f"dlseguid:{dlDNA4_dlseguid}"
    assert dlDNA5_dlseguid in slseguid("--AT--\nGATATC", table = table)
    assert cs("--AT--\nGATATC") == dlDNA5_dlseguid

    repr_from_tuple("AT", "CTATAG", 2)


def test_dcseguid():
    pUC19dna = Path("test_data/puc19.txt").read_text().strip()
    dcsg = "dcseguid:zhw8Yrxfo3FO5DDccx4PamBVPCQ"
    assert dcseguid(pUC19dna, rc(pUC19dna)) == dcsg
    w, c = Path("test_data/pUC19msg.txt").read_text().splitlines()
    assert dlseguid(w, c[::-1], 0) == "dlseguid:zhw8Yrxfo3FO5DDccx4PamBVPCQ"


def _test_speed():
    import timeit

    dna500 = """\
    CAGTGAAATCAGAACCCATGAGGGCGGACGAGTCATATCC
    GGTATTAGAGATTTATACAGTCTGGACACCTAGCGAACCG
    ACTTGAACCACCAGGATTGAAGACGAAACCTTAGAGTATA
    GTAATGCCGTACGTGTCGGGGCCCACGCATCTAGGACAGG
    ATCGCATGATGGTGGTTTTAGTTGCCGTTGTACCGGATTT
    CTTAGTAGTATAAGCATGAGGATAAGTGAAACCGGGTGAA
    GGTGGTTTGTGTGAGTGCCTAATAGTCCGACTCCCCGAGG
    GGAGTAGGCACTGCCTTCAGCGTTCAGTTATTGAGCACGT
    CCGCCCGGCGAAAGATGGCTTTGAGCTCCACTGACAGCCA
    GGGACCGCGTGCATGAGGCTAGAGCAGAGTCGTTGACAGT
    GAGATTAGATTGATCATTTTTATCTGAAACGGCAGCATAC
    CGACAGTTGTTCTCAAGCAAAGTGGTCTTGCCTAGATTCA
    ATATTGCCCACAATCAGCTC""".replace(
        "\n", ""
    )

    print("pure Python : ", end="")
    print(
        timeit.timeit(
            "cseguid(dna500, minrotation=min_rotation)",
            globals=globals(),
            number=1000,
        )
    )
    print("pydivsufsort: ", end="")
    print(timeit.timeit("cseguid(dna500)", globals=globals(), number=1000))
