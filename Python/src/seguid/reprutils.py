from seguid.tables import COMPLEMENT_TABLE
from textwrap import dedent
from seguid.asserts import assert_in_alphabet
from seguid.asserts import assert_anneal

def tuple_from_repr(
    rpr: str,
    table: dict = COMPLEMENT_TABLE,
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
    >>> s = \"""
    ...           -TATGCC
    ...           CATACG- \"""
    >>> tuple_from_repr(s)
    ('TATGCC', 'GCATAC', 1)
    >>> t = \"""
    ...                      -TATGCC
    ...                      CATACG- \"""
    >>> tuple_from_repr(s)
    ('TATGCC', 'GCATAC', 1)
    >>> tuple_from_repr(s) == tuple_from_repr(t)
    True
    """

    assert isinstance(space, str)
    assert len(space) == 1
    assert isinstance(sep, str)
    assert len(sep) == 1

    assert_in_alphabet(rpr, alphabet=set(table.keys()) | set(space) | set(sep) | set(" "))

    rpr_dedent = dedent(sep.join(ln for ln in rpr.split(sep) if ln.strip()))

    if sep not in rpr_dedent:
        raise ValueError(f"Expected two non-empty lines separated by {sep}")

    watson, crick = [x.rstrip(space) for x in rpr_dedent.split(sep)]

    overhang = (
        len(watson) - len(watson.lstrip(space)) - (len(crick) - len(crick.lstrip(space)))
    )

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