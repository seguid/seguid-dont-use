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

import hashlib
import base64
from typing import Callable
import warnings

from seguid.manip import rotate
from seguid.tables import COMPLEMENT_TABLE_DNA
from seguid.asserts import assert_in_alphabet
from seguid.asserts import assert_anneal

try:
    from pydivsufsort import min_rotation as mr
except ModuleNotFoundError:
    warnings.warn("pydivsufsort not found.", ImportWarning)  # TODO Is this the right way?
else:
    def min_rotation(s):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = mr(s)
        return result


def _seguid(seq: str,
            table: dict = COMPLEMENT_TABLE_DNA,
            encoding: callable = base64.standard_b64encode,
            prefix: str = "seguid:") -> str:

    assert isinstance(prefix, str)
    assert callable(encoding)

    assert_in_alphabet(seq, alphabet=set(table.keys()))
    m = hashlib.sha1()
    m.update(seq.encode("ASCII").upper())
    hs = encoding(m.digest())

    return f"{prefix}{hs.decode('ASCII').rstrip('=')}"


def seguid(seq: str,
           table: dict = COMPLEMENT_TABLE_DNA,
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
    """
    return _seguid(seq,
                   table=table,
                   encoding=base64.standard_b64encode,
                   prefix=prefix)


def slseguid(seq: str,
             table: dict = COMPLEMENT_TABLE_DNA,
             prefix: str = "slseguid:") -> str:
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
    """
    return _seguid(seq,
                   table=table,
                   encoding=base64.urlsafe_b64encode,
                   prefix=prefix)


def scseguid(seq: str,
             table: dict = COMPLEMENT_TABLE_DNA,
             min_rotation: Callable[[str], int] = min_rotation,
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
    >>> scseguid("ATTT")
    'scseguid:ot6JPLeAeMmfztW1736Kc6DAqlo'
    >>> slseguid("ATTT")
    'slseguid:ot6JPLeAeMmfztW1736Kc6DAqlo'
    >>> scseguid("TTTA")
    'scseguid:ot6JPLeAeMmfztW1736Kc6DAqlo'
    >>> slseguid("TTTA")
    'slseguid:8zCvKwyQAEsbPtC4yTV-pY0H93Q'
    """
    start = min_rotation(seq)

    return slseguid(rotate(seq, start),
                    table=table,
                    prefix=prefix)


def dlseguid(watson: str,
             crick: str,
             overhang: int,
             table: dict = COMPLEMENT_TABLE_DNA,
             prefix="dlseguid:"
             ) -> str:
    r"""SEGUID checksum for double stranded linear DNA (dlSEGUID).

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
    >>> dlseguid("TATGCC", "GCATAC", 1)
    'dlseguid:E7YtPGWjj3qCaPzWurlYBaJy_X4'
    >>> dlseguid("GCATAC", "TATGCC", 1)
    'dlseguid:E7YtPGWjj3qCaPzWurlYBaJy_X4'
    """
    assert_anneal(watson, crick, overhang, table=table)

    w, c, o = min(
        (
            (watson, crick, overhang),
            (crick, watson, len(watson) - len(crick) + overhang),
        )
    )

    space = "-"
    sep = "\n"

    msg = (
        f"{o*space}{w}{space*(-o+len(c)-len(w))}"
        f"{sep}"
        f"{-o*space}{c[::-1]}{space*(o+len(w)-len(c))}"
    ).rstrip()

    extable = table | {space: space, sep: sep}

    return slseguid(msg, table=extable, prefix=prefix)


def dcseguid(watson: str,
             crick: str,
             table: dict = COMPLEMENT_TABLE_DNA,
             min_rotation: Callable[[str], int] = min_rotation,
             prefix="dcseguid:") -> str:
    """SEGUID checksum for double stranded circular DNA (dcSEGUID).

    The dcSEGUID is the slSEGUID checksum calculated for the lexicographically
    smallest string rotation of a dsDNA sequence. Only defined for circular
    sequences. The min_rotation argument is a callable that takes a string as
    argument and returns another string.

    The checksum is prefixed with "dcseguid:"
    """
    assert len(watson) == len(crick)
    ln = len(watson)

    assert_anneal(watson, crick, 0, table=table)

    x = min_rotation(watson)
    y = min_rotation(crick)
    minwatson = rotate(watson, x)
    mincrick = rotate(crick, y)

    w, c = min(
        (minwatson, rotate(crick, ln - x)),
        (mincrick, rotate(watson, ln - y)),
    )

    return dlseguid(w, c, overhang=0, table=table, prefix=prefix)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
