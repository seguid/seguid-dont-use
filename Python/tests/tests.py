#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from hashlib import sha1
from base64 import urlsafe_b64encode as b64us

from seguid import slseguid
from seguid import scseguid
from seguid import dlseguid
from seguid import dcseguid

from seguid.chksum import seguid

from seguid.manip import rc
from seguid.manip import min_rotation_py

from seguid.reprutils import tuple_from_repr
from seguid.reprutils import repr_from_tuple

from seguid.asserts import assert_anneal
from seguid.asserts import assert_in_alphabet

from seguid.tables import COMPLEMENT_TABLE
from seguid.tables import COMPLEMENT_TABLE_RNA
from seguid.tables import COMPLEMENT_TABLE_IUPAC
from seguid.tables import TABLE_IUPAC_PROTEIN

import pytest

def test_assert_anneal():

    tuples = (("AT", "TA", 1),
              ("CTATAG", "AT", -2),
              ("AT", "CTATAG", 2),
              ("AT", "AT", 0))

    for watson, crick, overhang in tuples:
        assert_anneal(watson, crick, overhang)

    tuples = (("AT", "CG", 1),
              ("CTATAG", "AT", -3),
              ("AT", "CTATAG", 1))

    for watson, crick, overhang in tuples:
        with pytest.raises(ValueError):
            assert_anneal(watson, crick, overhang)

    with pytest.raises(AssertionError):
        assert_anneal("AT", "AT", 4)


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


def test_tuple_from_repr():
    """docstring."""

    rpr = """
      -TATGCC
      CATACG-
    """

    assert tuple_from_repr(rpr) == ("TATGCC", "GCATAC", 1)

    rpr = """
       TATGCC--
       -TACGGGG
    """

    assert tuple_from_repr(rpr) == ("TATGCC", "GGGGCAT", -1)

    rpr = """     # This stuff will give ValueError
        TATGCC
        ATACGG
    """

    with pytest.raises(ValueError):
        tuple_from_repr(rpr)

    rpr = """
       TATGCC--
       --ACGGGG
    """

    assert tuple_from_repr(rpr) == ("TATGCC", "GGGGCA", -2)

    rpr = """
       TATGCC
       ATACGG
    """

    assert tuple_from_repr(rpr) == ('TATGCC', 'GGCATA', 0)

    rpr_should_err = """
       - TGCC
       ATACGG
    """

    with pytest.raises(ValueError):
        tuple_from_repr(rpr_should_err)

    rpr_should_err = """
        -TGCC
       ATACGG
    """

    with pytest.raises(ValueError):
        tuple_from_repr(rpr_should_err)

    rpr_should_err = """
      ---TGCC
       ATACGG
    """

    with pytest.raises(ValueError):
        tuple_from_repr(rpr_should_err)

    rpr_should_err = """
      ---TGCC-
      -ATACGG-
    """

    tuple_from_repr(rpr_should_err)


def test_repr_from_tuple():
    assert repr_from_tuple(*("TATGCC", "GGGGCA", -2)) == "TATGCC--\n--ACGGGG"


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
    m13dna = Path("M13.txt").read_text().strip()
    sc = "scseguid:aAjgnsF9cPI6cu8IQ81sYnstVzU"
    assert scseguid(m13dna) == sc
    assert cs(Path("M13msg.txt").read_text().strip()) in sc


def test_dlseguid():
    ct = COMPLEMENT_TABLE
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
    pUC19dna = Path("puc19.txt").read_text().strip()
    dcsg = "dcseguid:zhw8Yrxfo3FO5DDccx4PamBVPCQ"
    assert dcseguid(pUC19dna, rc(pUC19dna)) == dcsg
    w, c = Path("pUC19msg.txt").read_text().splitlines()
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


if __name__ == "__main__":
    pytest.main([__file__, "-vvv", "-s"])
