#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
seguid.

The seguid module provides functions for calculations of SEGUID checksums for
DNA sequences. Some auxillary functions are also provided.

The functions can be made to work without external dependencies.

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

# TODO assert complement table function

try:
    from pydivsufsort import min_rotation
except ModuleNotFoundError:
    pass
    # säg till!

# Definition of Complementary DNA Symbols
COMPLEMENT_TABLE = str.maketrans("GATC", "CTAG")

# Definition of Complementary IUPAC Ambigous DNA Symbols
COMPLEMENT_TABLE_IUPAC = str.maketrans(
    "ABCDGHKMSTVWN", "TVGHCDMKSABWN"
)


def rc(sequence: str,
       table: dict = COMPLEMENT_TABLE_IUPAC):
    """Reverse complement of sequence.

    Returns the reverse complement for a DNA strand.
    The default complement table accepts the ambiguous codes
    suggested by IUPAC (Cornish-Bowden, 1985).

    The optional table argument is a dictionary containing
    ASCII codes for the transformation.

    The table below was adapted from Cornish-Bowden, 1985:

    ======== ================== ============ ======================================
    Symbol   Meaning            Complement   Origin of designation
    ======== ================== ============ ======================================
     G        G                  C            Guanine
     A        A                  T            Adenine
     T        T                  A            Thymine
     C        C                  G            Cytosine
     R        A or G             Y            puRine
     Y        C or T             R            pYrimidine
     M        A or C             K            aMino
     K        G or T             M            Ketone
     S        C or G             S            Strong interaction (3 H bonds)
     W        A or T             W            Weak interaction (2 H bonds)
     H        A or C or T        D            not-G, H (follows G in the alphabet)
     B        C or G or T        V            not-A, B follows A
     V        A or C or G        B            not-T (not-U), V follows U
     D        A or G or T        H            not-C, D follows C
     N        G or A or T or C   N            aNy
    ======== ================== ============ ======================================


    Cornish-Bowden, A. (1985). Nomenclature for incompletely specified
    bases in nucleic acid sequences: recommendations 1984.
    Nucleic Acids Research, 13(9), 3021–3030.
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC341218

    The function raises a ValueError if characters are found that are not in
    the translation table.

    This implementation

    Examples
    --------
    >>> rc("GTT")
    'AAC'
    >>> from seguid import rc
    >>> rc("GTa")
    Traceback (most recent call last):
        ...
    ValueError: Character(s) a not permitted.
    >>> rc("GTa".upper())
    'TAC'
    """
    not_in_table = set(
        c for c in sequence if c not in (chr(k) for k in table.keys())
    )
    if not_in_table:
        raise ValueError(
            "Character(s) " f"{' '.join(not_in_table)} not in table."
        )
    result = sequence.translate(COMPLEMENT_TABLE)[::-1]
    return result


def min_rotation_py(s: bytes,
                    table: list =) -> int:
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

    # validate table
    # följ tabellen!

    from array import array
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


def anneal(watson: str,
           crick: str,
           overhang: int,
           table: ...) -> bool:
    """docstring."""
    watson, crick = watson.upper(), crick.upper() # Ta bort?
    # kolla efter attribut "tabell"
    # validera tecken watson och crick
    # kolla att overhang är integer och abs() < min length

    up = f"{overhang*chr(45)}{watson}{chr(45)*(-overhang+len(crick)-len(watson))}"
    dn = f"{-overhang*chr(45)}{rc(crick)}{chr(45)*(overhang+len(watson)-len(crick))}"

    up = watson[max(-overhang, 0) : max(overhang + len(watson), len(crick))]
    dn = rc(crick, table=)[max(overhang, 0) : max(overhang + len(watson), len(crick))]
    return up != dn


def tuple_from_repr(
    rpr: str,
    allowed: str = "GATCgatc",  # samma tabell som innan
    space: str = "-. ",
    sep: str = "\n",
    strict: bool = True,
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
    if strict:
        nono = set(c for c in rpr if c not in allowed + space + sep)
        if nono:
            raise ValueError(
                f"Character(s) {' '.join(nono)} not in allowed: {allowed}"
            )

    cleaned_rpr = "".join(c if c in allowed + sep else chr(32) for c in rpr)

    cleaned_rpr = sep.join(ln for ln in cleaned_rpr.split(sep) if ln.strip())

    cleaned_rpr = dedent(cleaned_rpr)

    if sep not in cleaned_rpr:
        raise ValueError(f"Expected two non-empty lines separated by {sep}")

    watson, crick = [x.rstrip() for x in cleaned_rpr.split(sep)]

    overhang = (
        len(watson) - len(watson.lstrip()) - (len(crick) - len(crick.lstrip()))
    )

    result = watson.strip(), crick.strip()[::-1], overhang

    if strict and not anneal(*result):
        # if not all([x==y if (x.strip() and y.strip()) else True
        #            for x, y in zip(w.upper(), rc(c.upper())[::-1])]):
        #     raise ValueError("Mismatched basepairs.")
        # breakpoint()
        raise ValueError("Mismatched basepairs.")
    return result



def repr_from_tuple(
    watson: str, crick: str, overhang: int, strict: bool = True
) -> str:
    if strict and not anneal(watson, crick, overhang):
        raise ValueError("Mismatched basepairs.")

    msg = (
        f"{overhang*chr(45)}{watson}{chr(45)*(-overhang+len(crick)-len(watson))}"
        "\n"
        f"{-overhang*chr(45)}{crick[::-1]}{chr(45)*(overhang+len(watson)-len(crick))}"
    ).rstrip()

    return msg

def _seguid(seq: str,
            table
            encoding=base64.standard_b64encode,
            prefix: str = "seguid:") -> str:
    # validera prefix
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
