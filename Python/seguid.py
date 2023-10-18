#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MIT License.

Copyright (c) 2023 Björn Johansson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import hashlib
import base64
from textwrap import dedent

# The DNA complement can handle RNA (U>A)

ambiguous_dna_complement = {
    "A": "T", "C": "G", "G": "C",
    "T": "A", "M": "K", "R": "Y",
    "W": "W", "S": "S", "Y": "R",
    "K": "M", "V": "B", "H": "D",
    "D": "H", "B": "V", "X": "X",
    "N": "N", "U": "A",
}

_keys = "".join(ambiguous_dna_complement.keys())
_values = "".join(ambiguous_dna_complement.values())
_complement_table = str.maketrans(_keys + _keys.lower(),
                                  _values + _values.lower())

ambiguous_dna = "".join(ambiguous_dna_complement)

def rc(sequence: str):
    """Reverse complement.

    accepts mixed DNA/RNA
    """
    return sequence.translate(_complement_table)[::-1]


def seguid(seq: str) -> str:
    """Return the SEGUID (string) for a sequence (string).

    Given a nucleotide or amino-acid sequence (or any string),
    returns the SEGUID string (SEquence Globally Unique IDentifier).
    seq type = str.

    Note that the case is not important:

    >>> seguid("ACGTACGTACGT")
    'If6HIvcnRSQDVNiAoefAzySc6i4'
    >>> seguid("aaa")
    'YG7G6b2Kj/KtFOX63j8mRHHoIlE'

    The resulting string is not urs-safe.

    For more information about SEGUID, see:

    Babnigg, G., & Giometti, C. S. (2006). A database of unique protein
    sequence identifiers for proteome studies. Proteomics, 6(16), 4514–4522.
    https://doi.org/10.1002/pmic.200600032
    """
    m = hashlib.sha1()
    m.update(seq.encode().upper())
    hsh = base64.encodebytes(m.digest())
    return hsh.decode().rstrip("\n").rstrip("=")


def useguid(seq: str) -> str:
    """Url-safe SEGUID checksum for a sequence (string).

    This is the SEGUID checksum where the
    '+' and '/' characters of standard Base64 encoding are replaced by
    '-' and '_'.

    Examples
    --------
    >>> useguid("aaa")
    'YG7G6b2Kj_KtFOX63j8mRHHoIlE'
    """
    return seguid(seq).replace("+", "-").replace("/", "_")


def lseguid_blunt(seq: str) -> str:
    """Linear SEGUID checksum (lSEGUID) for a dsDNA sequence.

    lSEGUID is a checksum for a double stranded (dsDNA) DNA molecule
    represented by the seq argument. For the string `gtatgc`, there are two
    representations:

        5-gtatgc-3           5-gcatac-3
          ||||||               ||||||
        3-catacg-5           3-cgtatg-5

    Usually, only the top strand is used when referring to this molecule,
    This means that `gtatgc` and `gcatac` are equally valid representations.

    The lSEGUID algorithm selects the lexicographically smallest of the top
    and inferred bottom strands and applied the uSEGUID checksum.

    This means that:

        lseguid_blunt("gtatgc") == lseguid_blunt("gcatac")

    Examples
    --------
    >>> lseguid_blunt("gtatgc")
    'RAgd7GiTGrnLcI2VQ55u-lZiGsw'
    >>> lseguid_blunt("gcatac")
    'RAgd7GiTGrnLcI2VQ55u-lZiGsw'
    """
    return useguid(min(seq.upper(), str(rc(seq)).upper()))


def lseguid_sticky(watson: str, crick: str, overhang: int) -> str:
    """Linear SEGUID (lSEGUID) checksum for staggered linear dsDNA.

    Calculates the lSEGUID checksum for a ds DNA sequence
    described by two strings (watson and crick) representing the two
    complementary DNA strands and an integer describing the stagger
    between the two strands in the 5' end.

    The overhang is defined as the amount of 3' overhang in the 5'
    side of the molecule. A molecule with 5' overhang has a negative
    overhang value.

    See examples below:

        dsDNA       overhang

          nnn...    2
        nnnnn...

         nnnn...    1
        nnnnn...

        nnnnn...    0
        nnnnn...

        nnnnn...   -1
         nnnn...

        nnnnn...   -2
          nnn...


        dsDNA       overhang  lSEGUID

         TATGCC     1         Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU
        catacg

         gcatac     1         Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU
        CCGTAT


    Examples
    --------
    >>> lseguid_sticky("TATGCC", "gcatac", 1)
    'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
    >>> lseguid_sticky("gcatac", "TATGCC", 1)
    'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'

    """
    watson = watson.upper()
    crick = crick.upper()
    lw = len(watson)
    lc = len(crick)
    if overhang == 0 and lw == lc:
        return lseguid_blunt(watson)
    else:
        w, c, o = min(((watson, crick, overhang),
                       (crick, watson, lw - lc + overhang)))

    return useguid(f"{o*chr(32)}{w}{chr(10)}{-o*chr(32)}{c[::-1]}")


def smallest_rotation(s):
    """Smallest rotation of a string using a fast algorithm.

    This implementation uses python bindings to the libdivsufsort
    https://github.com/y-256/libdivsufsort C-algorithm by Yuta Mori.
    The python code (https://github.com/louisabraham/pydivsufsort)
    is distributed as a python package installable from PyPi:

         install pydivsufsort

    Examples
    --------
    >>> smallest_rotation("taaa")
    'aaat'
    """
    from pydivsufsort import min_rotation
    k = min_rotation(s)
    return s[k:] + s[:k]


def smallest_rotation_py(s):
    """Smallest rotation of a string, pure Python.

    Algorithm described in Pierre Duval, Jean. 1983. Factorizing Words
    over an Ordered Alphabet. Journal of Algorithms & Computational Technology
    4 (4) (December 1): 363–381. and Algorithms on strings and sequences based
    on Lyndon words, David Eppstein 2011.
    https://gist.github.com/dvberkel/1950267

    This is a pure python implementation, considerably slower than the
    smallest_rotation function above. This should only be used if a pure python
    implementation is requred.

    Examples
    --------
    >>> smallest_rotation_py("taaa")
    'aaat'

    """
    from array import array as _array
    prev, rep = None, 0
    ds = _array("u", 2 * s)
    lens, lends = len(s), len(ds)
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
                return "".join(w * rep)


def cseguid(seq: str, fun=smallest_rotation) -> str:
    """Circular SEGUID (cSEGUID) checksum for circular dsDNA.

    The cSEGUID is the uSEGUID checksum calculated for the lexicographically
    smallest string rotation of a dsDNA sequence. Only defined for circular
    sequences. The fun argument has to take a string as an argument and
    return another string.

    Examples
    --------
    >>> cseguid("attt")
    'oopV-6158nHJqedi8lsshIfcqYA'
    >>> cseguid("ttta")
    'oopV-6158nHJqedi8lsshIfcqYA'
    """
    return useguid(min(fun(seq.upper()), fun(str(rc(seq)).upper())))


def tuple_from_representation(rpr: str) -> tuple:
    """
    Generate a tuple from dsDNA text representation.

    This function can generate a tuple (watson, crick, overhang)
    from a dsDNA figure such as the one depicted below:

        5'-TATGCC-3'
           |||||
       3'-catacg-5'

    or:

        TATGCC
        |||||
       catacg

    or:

        TATGCC
       catacg

    The figures above all represent the same dsDNA molecule.

    The generated tuple can be used directly as an argument for the
    lseguid_sticky function.

    Examples
    --------
    >>> s = \"""
    ...         5'-TATGCC-3'
    ...            |||||
    ...        3'-catacg-5'  \"""
    >>> tuple_from_representation(s)
    ('TATGCC', 'gcatac', 1)

    """
    allowed = ambiguous_dna.upper() + ambiguous_dna.lower() + chr(10) + chr(32)

    cleaned_rpr = "".join(c if c in allowed else chr(32) for c in rpr)

    cleaned_rpr = "".join(ln for ln in
                          cleaned_rpr.splitlines(keepends=True) if ln.strip())

    cleaned_rpr = dedent(cleaned_rpr)

    w, c = [x.rstrip() for x in cleaned_rpr.splitlines()]

    overhang = len(w) - len(w.lstrip()) - (len(c) - len(c.lstrip()))

    return w.strip(), c.strip()[::-1], overhang


def dsseguid(watson: str,
             crick: str = "",
             overhang: int = 0,
             circular: bool = False,
             fun=smallest_rotation) -> str:
    """Double stranded SEGUID checksum for linear or circular dsDNA.

    Calculates the dsSEGUID checksum for a ds DNA sequence described by two
    strings (watson and crick) representing the two complementary DNA strands.

    The crick strand is an optional argument, if not given it is inferred from
    the watson strand.

    For linear DNA fragments, an optional overhang is defined as the amount
    of 3' overhang in the 5' side (left side ) of the molecule.
    A molecule with 5' overhang has a negative overhang value.

    See examples below:

        dsDNA       overhang

          nnn...    2
        nnnnn...

         nnnn...    1
        nnnnn...

        nnnnn...    0
        nnnnn...

        nnnnn...   -1
         nnnn...

        nnnnn...   -2
          nnn...

        dsDNA       overhang  SEGUID

         TATGCC     1         Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU
         |||||
        catacg

         gcatac     1         Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU
         |||||
        CCGTAT

        TATGCC      0         jZwQDKGSBpAF1wOh9icAv13hC2c
        ||||||
        atacgg

        ggcata      0         jZwQDKGSBpAF1wOh9icAv13hC2c
        ||||||
        CCGTAT


    If the optional boolean circular argument is set to True, the given
    sequence is interpreted as a circular double stranded DNA sequence.
    The overhang is assumed to be zero for circular sequences.

        dsDNA    dsseguid

         ATTT    oopV-6158nHJqedi8lsshIfcqYA
         ||||
         taaa

         aaat    oopV-6158nHJqedi8lsshIfcqYA
         ||||
         TTTA


    For a circular sequence dsSEGUID is the uSEGUID checksum calculated for
    the lexicographically smallest of the minimum rotation of the watson and
    the minimum rotation of the crick strand. The optional fun argument sets
    a function for the caclulation of minimum rotation for a string. The fun
    argument is expected to take a string as an argument and return another
    string.

    The dsSEGUID uses the url-safe Base64 encoding.

    Examples
    --------
    >>> dsseguid("TATGCC", "gcatac", 1)
    'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
    >>> dsseguid("gcatac", "TATGCC", 1)
    'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
    >>> dsseguid("TATGCC")
    'jZwQDKGSBpAF1wOh9icAv13hC2c'
    >>> dsseguid("ggcata")
    'jZwQDKGSBpAF1wOh9icAv13hC2c'
    >>> dsseguid("attt", circular=True)
    'oopV-6158nHJqedi8lsshIfcqYA'
    >>> dsseguid("ttta", circular=True)
    'oopV-6158nHJqedi8lsshIfcqYA'
    >>> dsseguid("ttta", overhang=1, circular=True)
    Traceback (most recent call last):
      File "/home/bjorn/mambaforge/envs/bjorn311/lib/python3.11/doctest.py", line 1351, in __run
        exec(compile(example.source, filename, "single",
      File "<doctest seguids.dsseguid[6]>", line 1, in <module>
        dsseguid("ttta", overhang=1, circular=True)
      File "/home/bjorn/Desktop/ProjectCheckSum/SEGUID-code-snippets/Python/seguids.py", line 430, in dsseguid
        raise ValueError("Circular sequence with overhang != 0.")
    ValueError: Circular sequence with overhang != 0.
    """

    watson = watson.upper()
    crick = crick.upper() if crick else rc(watson)

    if circular and overhang:
        raise ValueError("Circular sequence with overhang != 0.")

    if circular:
        msg = min(fun(watson), fun(crick))
    else:
        w, c, o = min(((watson, crick, overhang),
                       (crick, watson, len(watson) - len(crick) + overhang)))
        msg = f"{o*chr(32)}{w}{chr(10)}{-o*chr(32)}{c[::-1]}"
    return useguid(msg)

# https://github.com/ajalt/python-sha1


if __name__ == "__main__":

    import timeit

    dna500 = """\
    CAGTGAAATCAGAACCCATGAGGGCGGACGAGTCATATCC
    GGTATTAGAGATTTATACAGTCTGGACACCTAGCGAACCG
    ACTTGAACCACCAGGATTGAAGACGAAACCTTAGAGTATA
    GTAATGCCGTACGTGTCGGGGCCCACGCATCTAGGACAGG
    ATCGCATGATGGTGGTTTTAGTTGCCGTTGTACCGGATTT
    CTTAGTAGTATAAGCATGAGGATAAGTGAAACCGGGTGAA
    GGTGGTTTGTGTGAGTGCCTAATAGTCCGACTCCCCGAGG
    GGAGTAGGCACTGCCTTCAGCGTTCAGTTATTGAGCACGT
    CCGCCCGGCGAAAGATGGCTTTGAGCTCCACTGACAGCCA
    GGGACCGCGTGCATGAGGCTAGAGCAGAGTCGTTGACAGT
    GAGATTAGATTGATCATTTTTATCTGAAACGGCAGCATAC
    CGACAGTTGTTCTCAAGCAAAGTGGTCTTGCCTAGATTCA
    ATATTGCCCACAATCAGCTC""".replace("\n", "")

    print("pure Python : ", end="")
    print(timeit.timeit("cseguid(dna500, fun=smallest_rotation_py)",
                        globals=globals(),
                        number=1000))
    print("pydivsufsort: ", end="")
    print(timeit.timeit("cseguid(dna500)",
                        globals=globals(),
                        number=1000))






    # import base64, hashlib
    # data = b"aaa"
    # sha512_digest = hashlib.sha512(data).digest()
    # sha512t24u = base64.urlsafe_b64encode(sha512_digest[:24]).decode(" ascii")
    # sha512t24u
