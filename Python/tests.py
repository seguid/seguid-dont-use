from seguid import lseguid_blunt
from seguid import lseguid_sticky
from seguid import tuple_from_representation
from seguid import useguid
from seguid import cseguid
from seguid import rc
from seguid import dsseguid

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
    assert lseguid_sticky(*case2) == lseguid_blunt("gcatac") # Here dsSEGUID is different!

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


def test_tuple_from_representation():



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
