x = "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"


def test_lseguid_blunt():
    from seguids import lseguid_blunt
    assert (
        lseguid_blunt(x)
        == lseguid_blunt(x.upper())
        == lseguid_blunt(x.lower())
        == "bHrqalTJ793oAigMQ5_qCttJRTk"
    )

def test_lseguid_sticky():
    from seguids import lseguid_sticky


def test_useguid():
    from seguids import useguid
    assert useguid(x) == "cl5ukSUdlvZeBaBLEUhxisdRaL8"
    assert useguid(x.upper()) == "cl5ukSUdlvZeBaBLEUhxisdRaL8"
    assert useguid(x.lower()) == "cl5ukSUdlvZeBaBLEUhxisdRaL8"


def test_cseguid():
    from seguids import cseguid
    assert cseguid(x) == "naaZmDzyMa58OsNXROe5SvjC7WU"
    assert cseguid(x) == cseguid(x.upper())
    assert cseguid(x) == cseguid(x.lower())


def test_tuple_from_representation():

    from seguids import tuple_from_representation

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
