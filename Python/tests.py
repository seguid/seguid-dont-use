x = "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"


def test_lseguid_blunt():
    from seguid import lseguid_blunt
    assert (
        lseguid_blunt(x)
        == lseguid_blunt(x.upper())
        == lseguid_blunt(x.lower())
        == "bHrqalTJ793oAigMQ5_qCttJRTk"
    )

def test_lseguid_sticky():
    from seguid import lseguid_sticky
    from seguid import tuple_from_representation
    from seguid import useguid

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
    assert lseguid_sticky(*case2) == 'RAgd7GiTGrnLcI2VQ55u-lZiGsw'
    assert lseguid_sticky(*case2) == useguid("gcatac") # Here dsSEGUID is different!

    case3 = tuple_from_representation("""
      gTATGC
       atacgc

    """)
    assert case3 == ('gTATGC', 'cgcata', -1)

    assert lseguid_sticky(*case3) == '3PZNLU2tPDYs78AJP3mg5w1uEw4'

    assert lseguid_sticky(*case3) == useguid("cgcata\n CGTATg")


def test_useguid():
    from seguid import useguid
    assert useguid(x) == "cl5ukSUdlvZeBaBLEUhxisdRaL8"
    assert useguid(x.upper()) == "cl5ukSUdlvZeBaBLEUhxisdRaL8"
    assert useguid(x.lower()) == "cl5ukSUdlvZeBaBLEUhxisdRaL8"


def test_cseguid():
    from seguid import cseguid
    assert cseguid(x) == "naaZmDzyMa58OsNXROe5SvjC7WU"
    assert cseguid(x) == cseguid(x.upper())
    assert cseguid(x) == cseguid(x.lower())


def test_tuple_from_representation():

    from seguid import tuple_from_representation

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
