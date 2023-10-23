from seguid import lseguid_blunt
from seguid import lseguid_sticky
from seguid import tuple_from_representation
from seguid import useguid
from seguid import cseguid
from seguid import rc
from seguid import nseguid

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
    rprs = {}

    rprs[0] = """

        5'-TATGCC-3'
           |||||
       3'-catacg-5'

      """

    rprs[1] = """

        5'-TATGCC-3'
           |||||
        3-catacg-5'

      """

    rprs[2] = """
           TATGCC
           |||||
          catacg
      """
    rprs[3] = """
       TATGCC
      catacg

    """

    for i, rpr in rprs.items():

        assert tuple_from_representation(rpr) == ('TATGCC', 'gcatac', 1)

    rprs[0] = """

        5'-TATGCC-3'
           ||||||
        3'-atacgg-5'

      """

    rprs[1] = """

        5'-TATGCC-3'
           ||||||
         3-atacgg-5'

      """

    rprs[2] = """
           TATGCC
           ||||||
           atacgg
      """
    rprs[3] = """
       TATGCC
       atacgg
    """

    for i, rpr in rprs.items():
        assert tuple_from_representation(rpr) == ('TATGCC', 'ggcata', 0)

    rprs[0] = """

        5'-TATGCC-3'
           ||||||
        3'- tacgggg-5'

      """

    rprs[1] = """

        5'-TATGCC-3'
           ||||||
         3- tacgggg-5'

      """

    rprs[2] = """
           TATGCC
           ||||||
            tacgggg
      """
    rprs[3] = """
       TATGCC
        tacgggg
    """

    for i, rpr in rprs.items():
        assert tuple_from_representation(rpr) == ('TATGCC', 'ggggcat', -1)

    rprs[0] = """

        5'-TATGCC-3'
           ||||||
        3'-  acgggg-5'

      """

    rprs[1] = """

        5'-TATGCC-3'
           ||||||
         3-  acgggg-5'

      """

    rprs[2] = """
           TATGCC
           ||||||
             acgggg
      """
    rprs[3] = """
       TATGCC
         acgggg
    """

    for i, rpr in rprs.items():
        assert tuple_from_representation(rpr) == ('TATGCC', 'ggggca', -2)



def test_dcseguid():

    pUC19dna = Path("puc19.txt").read_text().strip()

    pUC19_scseguid = "n-NZfWfjHgA7wKoEBU6zfoXib_0"

    assert nseguid(pUC19dna, circular=True, ds=True) == pUC19_scseguid

    assert useguid(smallest_rotation(rc(pUC19dna))) == pUC19_scseguid


def test_scseguid():

    m13dna = Path("M13.txt").read_text().strip()

    m13_scseguid = 'aAjgnsF9cPI6cu8IQ81sYnstVzU'

    assert nseguid(m13dna, circular=True, ds=False) == m13_scseguid

    assert useguid(smallest_rotation(m13dna)) == m13_scseguid


def test_dlseguid():

    dlDNA = "AT"

    dlDNA_dlseguid = 'AWD-dt5-TEua8RbOWfnctJIu9nA'

    assert nseguid(dlDNA, circular=False, ds=True) == dlDNA_dlseguid

    assert useguid("AT\nTA") == dlDNA_dlseguid

    assert b64us(sha1("AT\nTA".encode("ASCII").upper()).digest()).decode("ASCII").rstrip("=") == dlDNA_dlseguid

    #  TA
    # TA

    dlDNA2 = ("TA", "TA", 1)

    dlDNA2_dlseguid = 'ZXq8qWnfgGJ69cxyVURZ-m69_8Y'

    assert nseguid(*("TA", "TA", 1), circular=False, ds=True) == dlDNA2_dlseguid

    assert useguid(" TA\nAT") == dlDNA2_dlseguid

    assert b64us(sha1(" TA\nAT".encode("ASCII").upper()).digest()).decode("ASCII").rstrip("=") == dlDNA2_dlseguid

    # TA
    #  TA

    dlDNA3 = ("TA", "TA", -1)

    dlDNA3_dlseguid = 'YrR7e0wVUeGwO8ftuRJLJfKlbPw'

    assert nseguid(*dlDNA3, circular=False, ds=True) == dlDNA3_dlseguid

    assert useguid("TA\n AT") == dlDNA3_dlseguid

    assert b64us(sha1("TA\n AT".encode("ASCII").upper()).digest()).decode("ASCII").rstrip("=") == dlDNA3_dlseguid

    # CTATAG
    #   TA

    dlDNA4 = ("CTATAG", "AT", -2)

    dlDNA4_dlseguid = '2njX4EGB3Is9Yvz3wmnAc3EagUI'

    assert nseguid(*dlDNA4, circular=False, ds=True) == dlDNA4_dlseguid

    assert useguid("  AT\nGATATC") == dlDNA4_dlseguid

    assert b64us(sha1("  AT\nGATATC".encode("ASCII").upper()).digest()).decode("ASCII").rstrip("=") == dlDNA4_dlseguid

    #   AT
    # GATATC

    dlDNA5 = ("AT", "CTATAG", 2)

    dlDNA5_dlseguid = '2njX4EGB3Is9Yvz3wmnAc3EagUI'

    assert nseguid(*dlDNA5, circular=False, ds=True) == dlDNA5_dlseguid

    assert useguid("  AT\nGATATC") == dlDNA5_dlseguid

    assert b64us(sha1("  AT\nGATATC".encode("ASCII").upper()).digest()).decode("ASCII").rstrip("=") == dlDNA5_dlseguid


def test_slseguid():

    slDNA = "AT"

    slDNA_dlseguid = b64us(sha1("AT".encode("ASCII").upper()).digest()).decode("ASCII").rstrip("=")

    assert nseguid(slDNA, circular=False, ds=False) == slDNA_dlseguid


x = "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"


def test_lseguid_blunt():
    assert (
        lseguid_blunt(x)
        == useguid(rc(x)+chr(10)+x[::-1])
        == 'R7VPKJXvozX-xPk0wFeNxbZd_dM'
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

    assert cseguid(x) == "naaZmDzyMa58OsNXROe5SvjC7WU"
    assert cseguid(x) == cseguid(x.upper())
    assert cseguid(x) == cseguid(x.lower())
