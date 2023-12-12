#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
seguid.

The seguid module provides four functions for calculations of SEGUID checksums
for biological sequences with varying topologies

slseguid (s)ingle-stranded (l)inear SEGUID
scseguid (s)ingle-stranded (c)ircular SEGUID
dlseguid (d)ouble-stranded (l)inear SEGUID
dcseguid (d)ouble-stranded (c)ircular SEGUID

Some auxillary functions are also provided.

The functions can be made to work without external dependencies, but
scseguid and dcseguid are considerably faster with pydivsufsort installed.

"""

__all__ = [
    "tuple_from_repr",
    "repr_from_tuple",
    "seguid",
    "slseguid",
    "scseguid",
    "dlseguid",
    "dcseguid",
]

import hashlib
import base64
from textwrap import dedent
from typing import Callable
from warnings import warn
from array import array
from seguid.manip import rc
from seguid.tables import COMPLEMENT_TABLE
from seguid.asserts import assert_in_alphabet

try:
    from pydivsufsort import min_rotation
except ModuleNotFoundError:
    warn("pydivsufsort not found.", ImportWarning)  # TODO Is this the right way?


def min_rotation_py(s: str,
                    table: dict = COMPLEMENT_TABLE) -> int:
    """Start position for the smallest rotation of a string s (pure Python).

    Algorithm described in:

    Pierre Duval, Jean. 1983. Factorizing Words
    over an Ordered Alphabet. Journal of Algorithms & Computational Technology
    4 (4) (December 1): 363–381. and Algorithms on strings and sequences based
    on Lyndon words, David Eppstein 2011.
    https://gist.github.com/dvberkel/1950267

    This is a pure python implementation, considerably slower than the
    min_rotation function from pydivsufsort that is imported at the top of
    this script.

    Should only be used if pydivsufsort.min_rotation can not be used.

    Note that both functions are case-sensitive and sorts by ASCII-code order
    or "ASCIIbetical" order so:

    - Uppercase come before lowercase letters; for example, "Z" precedes "a"
    - Digits and several punctuation marks come before letters.

    See the last two examples below for an example of the consequences of this.

    Examples
    --------
    >>> min_rotation_py("taaa")
    1
    >>> "taaa"[1:] + "taaa"[:1]
    'aaat'
    >>> s = "abaabaaabaababaaabaaababaab"
    >>> min_rotation_py(s)
    14
    >>> s[14:] + s[:14]
    'aaabaaababaababaabaaabaabab'
    >>> s = "abaabaaabaababaaabaaaBabaab"
    >>> min_rotation_py(s)
    21
    >>> s[21:] + s[:21]
    'Babaababaabaaabaababaaabaaa'
    >>> s= "abaabaaabaaba-baaabaaaBabaab"
    >>> min_rotation(s)
    13
    >>> s[13:] + s[:13]
    '-baaabaaaBabaababaabaaabaaba'
    """

    assert_table(table)
    assert_in_alphabet(s, alphabet = set(table.keys()))

    prev, rep = None, 0
    ds = array("u", 2 * s)
    lens = len(s)
    lends = lens * 2
    old = 0
    k = 0
    w = ""
    while k < lends:
        i, j = k, k + 1
        while j < lends and ds[i] <= ds[j]:
            i = (ds[i] == ds[j]) and i + 1 or k
            j += 1
        while k < i + 1:
            k += j - i
            prev = w
            w = ds[old:k]
            old = k
            if w == prev:
                rep += 1
            else:
                prev, rep = w, 1
            if len(w) * rep == lens:
                return old - i


def assert_anneal(watson: str,  #
           crick: str,
           overhang: int,
           table: dict = COMPLEMENT_TABLE) -> bool:
    """docstring."""

    assert_table(table)
    assert_in_alphabet(watson, alphabet=set(table.keys()))
    assert_in_alphabet(crick, alphabet=set(table.keys()))

    assert isinstance(overhang, int)
    assert -len(watson) < overhang
    assert overhang < len(crick)

    up = f"{overhang*chr(45)}{watson}{chr(45)*(-overhang+len(crick)-len(watson))}"
    dn = f"{-overhang*chr(45)}{rc(crick)}{chr(45)*(overhang+len(watson)-len(crick))}"

    up = watson[max(-overhang, 0): max(overhang + len(watson), len(crick))]
    dn = rc(crick, table=table)[max(overhang, 0): max(overhang + len(watson), len(crick))]

    if up != dn:
        raise ValueError("Mismatched basepairs.")


def tuple_from_repr(
    rpr: str,
    table: dict = COMPLEMENT_TABLE,
    space: str = "-. ",
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
    ...           catacg- \"""
    >>> tuple_from_repr(s)
    ('TATGCC', 'gcatac', 1)
    >>> t = \"""
    ...                      -TATGCC
    ...                      catacg- \"""
    >>> tuple_from_repr(s)
    ('TATGCC', 'gcatac', 1)
    >>> tuple_from_repr(s) == tuple_from_repr(t)
    True
    """

    assert isinstance(space, (str, set))
    assert isinstance(sep, (str, set))

    if not isinstance(space, set):
        space = set(space)

    if not isinstance(sep, set):
        sep = set(sep)

    assert_in_alphabet(rpr, alphabet=set(table.keys()) | space | sep)

    # cleaned_rpr = "".join(c if c in allowed + sep else " " for c in rpr)

    cleaned_rpr = sep.join(ln for ln in cleaned_rpr.split(sep) if ln.strip())

    cleaned_rpr = dedent(cleaned_rpr)

    if sep not in cleaned_rpr:
        raise ValueError(f"Expected two non-empty lines separated by {sep}")

    watson, crick = [x.rstrip() for x in cleaned_rpr.split(sep)]

    overhang = (
        len(watson) - len(watson.lstrip()) - (len(crick) - len(crick.lstrip()))
    )

    result = watson.strip(), crick.strip()[::-1], overhang

    assert_anneal(*result)

    return result


def repr_from_tuple(
    watson: str, crick: str, overhang: int
) -> str:

    assert_anneal(watson, crick, overhang)

    msg = (
        f"{overhang*chr(45)}{watson}{chr(45)*(-overhang+len(crick)-len(watson))}"
        "\n"
        f"{-overhang*chr(45)}{crick[::-1]}{chr(45)*(overhang+len(watson)-len(crick))}"
    ).rstrip()

    return msg

def _seguid(seq: str,
            table: dict = COMPLEMENT_TABLE,
            encoding=base64.standard_b64encode,
            prefix: str = "seguid:") -> str:
    assert isinstance(prefix, str)
    assert callable(encoding)
    assert_in_alphabet(seq, alphabet=set(table.keys()))
    m = hashlib.sha1()
    m.update(seq.encode("ASCII").upper())
    hs = encoding(m.digest())
    return f"{prefix}{hs.decode('ASCII').rstrip('=')}"

def seguid(seq: str,
           table
           prefix: str = "seguid:") -> str:
    """SEGUID checksum for protein or single stranded linear DNA.

    OBSOLETE, use slseguid instead.

    Given a nucleotide or amino-acid sequence `seq`, the function returns
    a string containing the SEquence Globally Unique IDentifier (SEGUID).

    The optional ´encoding` argument expects a function accepting a
    byte string an returning another byte string. Several such functions are
    available from the standard library:

    https://docs.python.org/3/library/base64.html

    The SEGUID is defined as the Base64 encoded sha1 checksum calculated for
    the sequence in upercase with the trailing "=" character removed. This
    means that upper or lower case symbols in `seq` do not affect the result.

    The resulting string is not url-safe as the Base64 encoding it potentially
    produces / and + characters, carrying special meaning in a Uniform Resource
    Locator (URL). It can also not be used as an identifier or variable name in
    programming languanges such as Python.

    This implementation follows the original SEGUID definition by
    Babnigg et al. 2006. For more information:

    Babnigg, G., & Giometti, C. S. (2006). A database of unique protein
    sequence identifiers for proteome studies. Proteomics, 6(16), 4514–4522.
    https://doi.org/10.1002/pmic.200600032

    The checksum is prefixed with "seguid:"

    Examples
    --------
    >>> seguid("AT")
    'seguid:Ax/RG6hzSrMEEWoCO1IWMGska+4'
    >>> seguid("at")
    'seguid:Ax/RG6hzSrMEEWoCO1IWMGska+4'
    """
    return _seguid(seq, encoding=base64.standard_b64encode, prefix=prefix)


def slseguid(seq: str, table, prefix: str = "slseguid:") -> str:
    """SEGUID checksum for single stranded linear DNA (slSEGUID).

    Identical to the seguid function except for that the '+' and '/' characters
    of standard Base64 encoding are replaced by '-' and '_', respectively
    following the standard in RFC 4648 section 5.

    The base64.urlsafe_b64encode from the Python standard libary is used.

    This checksum is applicable to single stranded linear DNA sequences.
    Can also be used for protein sequences.

    The checksum is prefixed with "slseguid:"

    Examples
    --------
    >>> slseguid("AT")
    'slseguid:Ax_RG6hzSrMEEWoCO1IWMGska-4'
    >>> slseguid("at")
    'slseguid:Ax_RG6hzSrMEEWoCO1IWMGska-4'
    """
    return _seguid(seq, encoding=base64.urlsafe_b64encode, prefix=prefix)

def scseguid(seq: str,
             table
             min_rotation: Callable[[str, dict?], int] = min_rotation, # hur fixa kwargs
             prefix="scseguid:") -> str:
    r"""SEGUID checksum for single stranded circular DNA (scSEGUID).

    The scSEGUID is the slSEGUID checksum calculated for the lexicographically
    smallest string rotation of a ssDNA sequence.

    Only defined for circular sequences.

    The srfun argument has to take a string as an argument and
    return another string.

    The checksum is prefixed with "scseguid:"

    Examples
    --------
    >>> scseguid("attt")
    'scseguid:ot6JPLeAeMmfztW1736Kc6DAqlo'
    >>> slseguid("attt")
    'slseguid:ot6JPLeAeMmfztW1736Kc6DAqlo'
    >>> scseguid("ttta")
    'scseguid:ot6JPLeAeMmfztW1736Kc6DAqlo'
    >>> slseguid("ttta")
    'slseguid:8zCvKwyQAEsbPtC4yTV-pY0H93Q'
    """
    # validera ?
    start = min_rotation(seq, table) # TODO tysta varningen
    return slseguid(seq[start:] + seq[:start], prefix=prefix, table)


def dlseguid(watson: str,
             crick: str,
             overhang: int,
             table
             prefix="dlseguid:") -> str:
    r"""SEGUID checksum for double stranded linear DNA (dlSEGUID)

    Calculates the dlSEGUID checksum for a dsDNA sequence defined by two
    strings representing the upper (Watson) and lower (Crick) strand
    complementary DNA strands and an integer value describing the stagger
    between the two strands in the 5' (left) end of the molecule.

    The overhang is defined as the amount of 3' overhang at the start side
    of the molecule. A molecule with 5' overhang has a negative
    overhang value.

    See examples below:

    ::


        dsDNA       overhang

        --nnn...    2
        nnnnn...

        -nnnn...    1
        nnnnn...

        nnnnn...    0
        nnnnn...

        nnnnn...   -1
        -nnnn...

        nnnnn...   -2
        --nnn...


    The algorithm first selects the lexicographically smallest
    of the top or bottom strands.

    For positive overhang, the top strand is is left padded with the number
    of hyphen characters (ASCII 45) indicated by the overhang value.

    For negative overhang the reverse of the bottom strand is similarly padded.

    The string pair is similarly padded on the other side, so that two eaqual
    length strings are formed.

    The two string are joined, separated by a line break (ASCII 10) and the
    lsSEGUID function is used on the resulting string.

    ::

        dsDNA       overhang  dlSEGUID

        -TATGCC     1        Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU
         |||||
        catacg-

        -gcatac     1        Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU
         |||||
        CCGTAT-

    For the linear dsDNA sequence defined by Watson = "TATGCC", Crick ="gcatac"
    and overhang = 1 (see figures above), The "gcatac" strand is selected as
    "gcatac" < "TATGCC".

    Overhang is positive, so the first strand is padded to "-gcatac".

    A string is constructed like so:
    ::

        "-gcatac" + chr(10) + "CCGTAT-"

    The checksum is prefixed with "dlseguid:"

    Examples
    --------
    >>> dlseguid("TATGCC", "gcatac", 1)
    'dlseguid:E7YtPGWjj3qCaPzWurlYBaJy_X4'
    >>> dlseguid("gcatac", "TATGCC", 1)
    'dlseguid:E7YtPGWjj3qCaPzWurlYBaJy_X4'
    >>> slseguid("-gcatac\nCCGTAT-")
    'slseguid:E7YtPGWjj3qCaPzWurlYBaJy_X4'
    >>> dlseguid("gcatac", "GTATGC", 0)
    'dlseguid:b0Xa5pLe4LNd5T8fhGWHicCI_f4'
    >>> slseguid("gcatac\nCGTATG")
    'slseguid:b0Xa5pLe4LNd5T8fhGWHicCI_f4'
    """
    watson, crick = watson.upper(), crick.upper() # bort!
    # validera
    w, c, o = min(
        (
            (watson, crick, overhang),
            (crick, watson, len(watson) - len(crick) + overhang),
        )
    )
    spacer = "-"
    msg = (
        f"{o*spacer}{w}{spacer*(-o+len(c)-len(w))}"
        "\n"
        f"{-o*spacer}{c[::-1]}{spacer*(o+len(w)-len(c))}"
    ).rstrip()
    # lägg till newline & spacer till tabell gör ny tabell

    return slseguid(msg, table=extendedtable, prefix=prefix)


def dcseguid(watson: str,
             crick: str,
             min_rotation: Callable[[str], int] = min_rotation,
             table,
             prefix="dcseguid:") -> str:
    """SEGUID checksum for double stranded circular DNA (dcSEGUID).

    The dcSEGUID is the slSEGUID checksum calculated for the lexicographically
    smallest string rotation of a dsDNA sequence. Only defined for circular
    sequences. The min_rotation argument is a callable that takes a string as
    argument and returns another string.

    The checksum is prefixed with "dcseguid:"
    """
    # validera att watson lika lang som crick
    watson, crick = watson.upper(), crick.upper()
    lw = len(watson)
    x = min_rotation(watson, table)
    y = min_rotation(crick, table )
    #watson2 = rotate(watson, lw - y)
    #watson1 = rotate(watson, x)
    #watson2 = rotate(watson1, lw - y -x)

    w, c = min(
        (watson[x:] + watson[:x], crick[lw - x :] + crick[: lw - x]),
        (crick[y:] + crick[:y], watson[lw - y :] + watson[: lw - y]),
    )

    return dlseguid(w, c, overhang=0, table, prefix=prefix)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
