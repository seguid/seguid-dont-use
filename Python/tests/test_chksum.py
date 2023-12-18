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

from seguid.chksum import seguid

from seguid.manip import reverse
from seguid.manip import rc
from seguid.manip import min_rotation_py
from seguid.manip import rotate_to_min

from seguid.reprutils import repr_from_tuple

from seguid.tables import COMPLEMENT_TABLE_DNA
from seguid.tables import TABLE_IUPAC_PROTEIN


def test_min_rotation():
    """Tests for the pydivsufsort min_rotation"""
    from pydivsufsort import min_rotation

    def smallest_rotation(s):
        i = min_rotation(s)
        return s[i:] + s[:i]

    assert rotate_to_min("taaa") == "aaat"
    assert (
        rotate_to_min("abaabaaabaababaaabaaababaab")
        == "aaabaaababaababaabaaabaabab"
    )
    assert (
        rotate_to_min("abaabaaabaababaaabaaaBabaab")
        == "Babaababaabaaabaababaaabaaa"
    )


def test_min_rotation_py():
    """docstring."""

    assert rotate_to_min("TAAA", min_rotation = min_rotation_py) == "AAAT"
    assert (
        rotate_to_min("ACAACAAACAACACAAACAAACACAAC",
                      min_rotation = min_rotation_py)
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
    assert dlseguid(w, reverse(c), 0) == "dlseguid:zhw8Yrxfo3FO5DDccx4PamBVPCQ"
