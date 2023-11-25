from seguid import lseguid_blunt
from seguid import lseguid_sticky
from seguid import tuple_from_representation
from seguid import useguid
from seguid import cseguid
from seguid import rc
from seguid import _nseguid

from seguid import smallest_rotation
from seguid import smallest_rotation_py

from pathlib import Path

from hashlib import sha1
from base64 import urlsafe_b64encode as b64us

import pytest


def test_rc():
    """docstring."""
    assert rc("Gat") == 'atC'
    assert rc("GTT", strict=True) == 'AAC'

    with pytest.raises(ValueError):
        rc("GTZ", strict=True)


def test_smallest_rotation():
    """docstring."""
    assert smallest_rotation("taaa") == 'aaat'
    assert smallest_rotation("abaabaaabaababaaabaaababaab") == 'aaabaaababaababaabaaabaabab'
    assert smallest_rotation("abaabaaabaababaaabaaaBabaab") == 'Babaababaabaaabaababaaabaaa'


def test_smallest_rotation_py():
    """docstring."""
    assert smallest_rotation_py("taaa") == 'aaat'
    assert smallest_rotation_py("abaabaaabaababaaabaaababaab") == 'aaabaaababaababaabaaabaabab'
    assert smallest_rotation_py("abaabaaabaababaaabaaaBabaab") == 'Babaababaabaaabaababaaabaaa'


def test_tuple_from_representation():
    """docstring."""

    rpr = """
       TATGCC
      catacg
    """

    assert tuple_from_representation(rpr) == ('TATGCC', 'gcatac', 1)

    rpr = """
       TATGCC
       atacgg
    """

    assert tuple_from_representation(rpr) == ('TATGCC', 'ggcata', 0)

    rpr = """
       TATGCC
        tacgggg
    """

    assert tuple_from_representation(rpr) == ('TATGCC', 'ggggcat', -1)

    rpr = """
       TATGCC
         acgggg
    """

    assert tuple_from_representation(rpr) == ('TATGCC', 'ggggca', -2)


def cs(arg):
    return b64us(sha1(arg.encode("ASCII")).digest()).decode("ASCII").rstrip("=")


def test_dcseguid():

    pUC19dna = Path("puc19.txt").read_text().strip()

    sg = 'TjUtGuJLLCwtboFklW_FlQUAyWo'

    assert _nseguid(pUC19dna, circular=True, ds=True) == sg

    assert useguid(Path("pUC19msg.txt").read_text().strip()) == sg


def test_scseguid():

    m13dna = Path("M13.txt").read_text().strip()

    sg = 'MXhEZHw_nj5AlSnYstCXSluezVU'

    assert _nseguid(m13dna, circular=True, ds=False) == sg

    assert cs(Path("M13msg.txt").read_text().strip()) == sg


def test_dlseguid():

    dlDNA = "AT"

    dlDNA_dlseguid = 'AWD-dt5-TEua8RbOWfnctJIu9nA'

    assert _nseguid(dlDNA, circular=False, ds=True) == dlDNA_dlseguid

    assert useguid("AT\nTA") == dlDNA_dlseguid

    assert cs("AT\nTA") == dlDNA_dlseguid

    #  TA
    # TA

    dlDNA2 = ("TA", "TA", 1)

    dlDNA2_dlseguid = 'ZXq8qWnfgGJ69cxyVURZ-m69_8Y'

    assert _nseguid(*dlDNA2, circular=False, ds=True) == dlDNA2_dlseguid

    assert useguid(" TA\nAT") == dlDNA2_dlseguid

    assert cs(" TA\nAT") == dlDNA2_dlseguid

    # TA
    #  TA

    dlDNA3 = ("TA", "TA", -1)

    dlDNA3_dlseguid = 'YrR7e0wVUeGwO8ftuRJLJfKlbPw'

    assert _nseguid(*dlDNA3, circular=False, ds=True) == dlDNA3_dlseguid

    assert useguid("TA\n AT") == dlDNA3_dlseguid

    assert cs("TA\n AT") == dlDNA3_dlseguid

    # CTATAG
    #   TA

    dlDNA4 = ("CTATAG", "AT", -2)

    dlDNA4_dlseguid = '2njX4EGB3Is9Yvz3wmnAc3EagUI'

    assert _nseguid(*dlDNA4, circular=False, ds=True) == dlDNA4_dlseguid

    assert useguid("  AT\nGATATC") == dlDNA4_dlseguid

    assert cs("  AT\nGATATC") == dlDNA4_dlseguid

    #   AT
    # GATATC

    dlDNA5 = ("AT", "CTATAG", 2)

    dlDNA5_dlseguid = '2njX4EGB3Is9Yvz3wmnAc3EagUI'

    assert _nseguid(*dlDNA5, circular=False, ds=True) == dlDNA5_dlseguid

    assert useguid("  AT\nGATATC") == dlDNA5_dlseguid

    assert cs("  AT\nGATATC") == dlDNA5_dlseguid


def test_slseguid():

    slDNA = "AT"

    slDNA_slseguid = "Ax_RG6hzSrMEEWoCO1IWMGska-4"

    assert _nseguid(slDNA, circular=False, ds=False) == slDNA_slseguid

    assert cs("AT") == slDNA_slseguid


x = "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"


def test_lseguid_blunt():
    assert (
        lseguid_blunt(x) ==
        useguid(rc(x)+chr(10)+x[::-1]) ==
        'R7VPKJXvozX-xPk0wFeNxbZd_dM'
    )

def test_lseguid_sticky():

    case1 = tuple_from_representation("""
       TATGCC
      catacg
    """)
    assert case1 == ('TATGCC', 'gcatac', 1)
    assert lseguid_sticky(*case1) == 'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
    assert lseguid_sticky(*case1) == useguid(" gcatac\nCCGTAT")

    case2 = tuple_from_representation("""
      gTATGC
      catacg
    """)
    assert case2 == ('gTATGC', 'gcatac', 0)
    assert lseguid_sticky(*case2) == 'b0Xa5pLe4LNd5T8fhGWHicCI_f4'
    assert lseguid_sticky(*case2) == lseguid_blunt("gcatac") # Here nseguid is different!

    case3 = tuple_from_representation("""
      gTATGC
       atacgc

    """)
    assert case3 == ('gTATGC', 'cgcata', -1)

    assert lseguid_sticky(*case3) == '3PZNLU2tPDYs78AJP3mg5w1uEw4'

    assert lseguid_sticky(*case3) == useguid("cgcata\n CGTATg")


def test_useguid():

    assert useguid(x) == "cl5ukSUdlvZeBaBLEUhxisdRaL8"
    assert useguid(x.upper()) == "cl5ukSUdlvZeBaBLEUhxisdRaL8"
    assert useguid(x.lower()) == "cl5ukSUdlvZeBaBLEUhxisdRaL8"

NP_313053_1 = "MKALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSGASRGIRLLQEEEEGLPLVGRVAAGEPLLAQQHIEGHYQVDPSLFKPNADFLLRVSGMSMKDIGIMDGDLLAVHKTQDVRNGQVVVARIDDEVTVKRLKKQGNKVELLPENSEFKPIVVDLRQQSFTIEGLAVGVIRNGDWL"

def test_cseguid():

    assert cseguid(x) == 'zH8Vfbq2vCWgr-LyPj2pYoMCAJs'
    assert cseguid(x) == cseguid(x.upper())
    assert cseguid(x) == cseguid(x.lower())


def _____test_speed():

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
    ATATTGCCCACAATCAGCTC""".replace("\n", "")

    print("pure Python : ", end="")
    print(timeit.timeit("cseguid(dna500, minrotation=min_rotation)",
                        globals=globals(),
                        number=1000))
    print("pydivsufsort: ", end="")
    print(timeit.timeit("cseguid(dna500)",
                        globals=globals(),
                        number=1000))


if __name__ == "__main__":
    pytest.main([__file__, "-vvv", "-s"])
