from seguid.tables import COMPLEMENT_TABLE_DNA
from textwrap import dedent
from seguid.asserts import assert_in_alphabet
from seguid.asserts import assert_anneal

def tuple_from_repr(
    rpr: str,
    table: dict = COMPLEMENT_TABLE_DNA,
    space: str = "-",
    sep: str = "\n"
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
    # rpr = """
    #   ---TGCC-
    #   -ATACGG-
    # """
    assert isinstance(space, str)
    assert len(space) == 1
    assert isinstance(sep, str)
    assert len(sep) == 1

    ws = " "

    assert_in_alphabet(rpr, alphabet=set(table.keys()) | set(space) | set(sep) | set(ws))

    rpr_dedent = dedent(sep.join(ln for ln in rpr.split(sep) if ln.strip(ws)))

    if sep not in rpr_dedent:
        raise ValueError(f"Expected two non-empty lines separated by {sep}")

    w, c = [x.strip(ws) for x in rpr_dedent.split(sep)]
    assert not (w.startswith(space) and c.startswith(space))
    assert not (w.endswith(space) and c.endswith(space))

    watson, crick = [x.rstrip(space + ws) for x in rpr_dedent.split(sep)]

    overhang = (
        len(watson) - len(watson.lstrip(space)) - (len(crick) - len(crick.lstrip(space)))
    )

    assert not (watson.startswith(space) and crick.startswith(space))


    result = watson.strip(space), crick.strip(space)[::-1], overhang

    assert_anneal(*result, table=table | {c: c for c in space + sep})

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
