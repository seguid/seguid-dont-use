from seguid.tables import COMPLEMENT_TABLE_DNA
# from textwrap import dedent
from inspect import cleandoc
from seguid.asserts import assert_in_alphabet
from seguid.asserts import assert_anneal
from string import whitespace


def tuple_from_repr(
    rpr: str,
    table: dict = COMPLEMENT_TABLE_DNA,
    space: str = "-",
) -> tuple:
    """Generate a tuple from dsDNA text representation.

    This function can generate a tuple (watson, crick, overhang)
    from a dsDNA figure such as the ones depicted below. The resulting
    tuple can be used as an argument for the lSEGUID_sticky or nseguid
    functions. See these functions for the definition of watson, crick and
    overhang.
    ::

              -TATGCC
              catacg-


    The three figures above represent the same dsDNA molecule.


    Examples
    --------
    >>> rpr = \"""
    ...           -TATGCC
    ...           CATACG- \"""
    >>> tuple_from_repr(rpr)
    ('TATGCC', 'GCATAC', 1)
    >>> rpr2 = \"""
    ...                      -TATGCC
    ...                      CATACG- \"""
    >>> tuple_from_repr(rpr2)
    ('TATGCC', 'GCATAC', 1)
    >>> tuple_from_repr(rpr) == tuple_from_repr(rpr2)
    True
    """
    space = "-"

    assert isinstance(space, str)
    assert len(space) == 1

    assert_in_alphabet(rpr, alphabet=set(table.keys()) | set(space) | set(whitespace))

    watson, crickrv = cleandoc(rpr).split()

    assert len(watson) == len(crickrv)
    assert not (watson.startswith(space) and crickrv.startswith(space))
    assert not (watson.endswith(space) and crickrv.endswith(space))

    overhang = (
        len(watson) - len(watson.lstrip(space)) - (len(crickrv) - len(crickrv.lstrip(space)))
    )

    result = watson.strip(space), crickrv.strip(space)[::-1], overhang
    assert_anneal(*result, table=table | {space: space})

    return result


def repr_from_tuple(
    watson: str, crick: str, overhang: int
) -> str:
    """docstring."""
    assert_anneal(watson, crick, overhang)

    msg = (
        f"{overhang*chr(45)}{watson}{chr(45)*(-overhang+len(crick)-len(watson))}"
        "\n"
        f"{-overhang*chr(45)}{crick[::-1]}{chr(45)*(overhang+len(watson)-len(crick))}"
    ).rstrip()

    return msg


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
