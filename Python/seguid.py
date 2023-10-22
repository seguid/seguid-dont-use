#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The seguid module provides ten functions for calculations of checksums of
biological sequences. Som auxillary functions are also provided.

"""

import hashlib
import base64
from textwrap import dedent

# https://github.com/ajalt/python-sha1

# Definition of Complementary IUPAC Ambigous DNA Symbols
COMPLEMENT_TABLE = str.maketrans("ABCDGHKMSTVWNabcdghkmstvwn",
                                 "TVGHCDMKSABWNtvghcdmksabwn")


def rc(sequence: str,
       table: dict = COMPLEMENT_TABLE,
       strict: bool = False):
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

    The optional boolean strict argument can optionally check for
    characters not defined in the translation table and raise a
    ValueError if such are encountered in the sequence.

    This implementation

    Examples
    --------
    >>> rc("Gat")
    'atC'
    >>> rc("GTT", strict=True)
    'AAC'
    >>> from seguid import rc
    >>> rc("GTZ", strict=True)
    Traceback (most recent call last):
        ...
    ValueError: Character(s) Z not permitted.
    >>>
    """
    if strict:
        not_in_table = set(c for c in sequence
                           if c not in (chr(k) for k in table.keys()))
        if not_in_table:
            raise ValueError("Character(s) "
                             f"{' '.join(not_in_table)} not permitted.")
    return sequence.translate(COMPLEMENT_TABLE)[::-1]


def smallest_rotation(s):
    """Smallest rotation of a string using a fast C-algorithm.

    This implementation uses the pydivsufsort python extension for the
    libdivsufsort algorithm by Yuta Mori, written in C.

    The python extension (https://github.com/louisabraham/pydivsufsort)
    is available as a python package installable from PyPi:
    ::

         pip install pydivsufsort

    The C-code for libdivsufsort can be found here:
    https://github.com/y-256/libdivsufsort

    Note that this function is case-sensitive.

    Examples
    --------
    >>> smallest_rotation("taaa")
    'aaat'
    >>> smallest_rotation("abaabaaabaababaaabaaababaab")
    'aaabaaababaababaabaaabaabab'
    >>> smallest_rotation("abaabaaabaababaaabaaaBabaab")
    'Babaababaabaaabaababaaabaaa'
    """
    from pydivsufsort import min_rotation
    k = min_rotation(s)
    return s[k:] + s[:k]


def smallest_rotation_py(s):
    """Smallest rotation of a string, pure Python.

    Algorithm described in:

    Pierre Duval, Jean. 1983. Factorizing Words
    over an Ordered Alphabet. Journal of Algorithms & Computational Technology
    4 (4) (December 1): 363–381. and Algorithms on strings and sequences based
    on Lyndon words, David Eppstein 2011.
    https://gist.github.com/dvberkel/1950267

    This is a pure python implementation, considerably slower than the
    smallest_rotation function above.

    This should only be used if a pure python implementation is required.

    Note that this function is case-sensitive.

    Examples
    --------
    >>> smallest_rotation_py("taaa")
    'aaat'
    >>> smallest_rotation_py("abaabaaabaababaaabaaababaab")
    'aaabaaababaababaabaaabaabab'
    >>> smallest_rotation_py("abaabaaabaababaaabaaaBabaab")
    'Babaababaabaaabaababaaabaaa'
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


def tuple_from_representation(rpr: str,
                              allowed: str = ("ABCDGHKMSTVWNabcdghkmstvwn"
                                              + chr(10) + chr(32))) -> tuple:
    """Generate a tuple from dsDNA text representation.

    This function can generate a tuple (watson, crick, overhang)
    from a dsDNA figure such as the ones depicted below. The resulting
    tuple can be used as an argument for the lSEGUID_sticky or nseguid
    functions. See these functions for the definition of watson, crick and
    overhang.
    ::


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


    The three figures above represent the same dsDNA molecule.


    Examples
    --------
    >>> s = \"""
    ...         5'-TATGCC-3'
    ...            |||||
    ...        3'-catacg-5'  \"""
    >>> tuple_from_representation(s)
    ('TATGCC', 'gcatac', 1)
    """
    cleaned_rpr = "".join(c if c in allowed else chr(32) for c in rpr)

    cleaned_rpr = "".join(ln for ln in
                          cleaned_rpr.splitlines(keepends=True) if ln.strip())

    cleaned_rpr = dedent(cleaned_rpr)

    w, c = [x.rstrip() for x in cleaned_rpr.splitlines()]

    overhang = len(w) - len(w.lstrip()) - (len(c) - len(c.lstrip()))

    return w.strip(), c.strip()[::-1], overhang

def seguid(seq: str,
           encoding=base64.standard_b64encode) -> str:
    """Return the SEGUID (string) for a sequence (string).

    Given a nucleotide or amino-acid sequence `seq` returns
    a string containing the SEquence Globally Unique IDentifier (SEGUID).

    The optional ´encoding` argument expects a function accepting a
    byte string an returning another byte string with the encoded string.
    Several such functions are available from the standard library:

    https://docs.python.org/3/library/base64.html

    The SEGUID is defined as the Base64 encoded sha1 checksum calculated for
    the sequence in upercase with the trailing "=" character removed. This
    means that upper or lower case symbols in `seq` do not affect the result.

    The resulting string is not url-safe as the Base64 encoding it potentially
    produces / and + characters, carrying special meaning in a Uniform Resource
    Locator (url). It can also not be used as an identifier or variable name in
    programming languanges such as Python.

    This implementation follows the original SEGUID definition by
    Babnigg et al. 2006. For more information:

    Babnigg, G., & Giometti, C. S. (2006). A database of unique protein
    sequence identifiers for proteome studies. Proteomics, 6(16), 4514–4522.
    https://doi.org/10.1002/pmic.200600032

    Examples
    --------
    >>> seguid("AT")
    'Ax/RG6hzSrMEEWoCO1IWMGska+4'
    >>> seguid("at")
    'Ax/RG6hzSrMEEWoCO1IWMGska+4'
    """
    m = hashlib.sha1()
    m.update(seq.encode("ASCII").upper())
    hs = encoding(m.digest())
    return hs.decode("ASCII").rstrip("=")


def useguid(seq: str) -> str:
    """Url-safe SEGUID checksum for a sequence (string).

    This is the SEGUID checksum where the '+' and '/' characters of standard
    Base64 encoding are replaced by '-' and '_', respectively following
    the standard in RFC 4648 section 5.

    Examples
    --------
    >>> useguid("AT")
    'Ax_RG6hzSrMEEWoCO1IWMGska-4'
    >>> seguid("AT")
    'Ax/RG6hzSrMEEWoCO1IWMGska+4'
    """
    return seguid(seq, encoding=base64.urlsafe_b64encode)


def lseguid_blunt(seq: str) -> str:
    r"""Linear SEGUID checksum (lSEGUID) for a blunt dsDNA sequence.

    lSEGUID is a checksum for a blunt double stranded (dsDNA) DNA molecule
    defined by *one* of the DNA strands in 5'-3' (`seq`). The complementary
    strand is deduced from the given strand.

    For example, the seq argument `GTATGC` represents the blunt dsDNA
    molecule depicted in the two figures below:
    ::

            5-GTATGC-3           5-gcatac-3
              ||||||               ||||||
            3-catacg-5           3-CGTATG-5


    The lSEGUID function constructs a string made from the lexicographically
    smallest of the top or bottom strands, a line break (ASCII 10) and the
    other complementary strand in reverse.

    The resulting string is used as an argument for the useguid function and
    the resulting checksum is returned.

    For the example above, gcatac < GTATGC, so the string used for the
    calculation will be:
    ::

        "gcatac" + chr(10) + "CGTATG"


    This means that lseguid_blunt("GTATGC") == lseguid_blunt("gcatac")


    Examples
    --------
    >>> lseguid_blunt("gtatgc")
    'b0Xa5pLe4LNd5T8fhGWHicCI_f4'
    >>> lseguid_blunt("GTATGC")
    'b0Xa5pLe4LNd5T8fhGWHicCI_f4'
    >>> useguid("gcatac" + chr(10) + "CGTATG")
    'b0Xa5pLe4LNd5T8fhGWHicCI_f4'
    """
    return nseguid(watson=seq, crick="", overhang=0)


def lseguid_sticky(watson: str, crick: str, overhang: int) -> str:
    """Linear SEGUID (lSEGUID) checksum for staggered linear dsDNA.

    Calculates the lSEGUID checksum for a dsDNA sequence defined by two
    strings (watson and crick) representing two complementary DNA strands
    and an integer describing the stagger between the two strands in the 5'
    (left) end of the molecule.

    The overhang is defined as the amount of 3' overhang in the 5'
    side of the molecule. A molecule with 5' overhang has a negative
    overhang value.

    See examples below:
    ::


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


    The lSEGUID_sticky function first selects the lexicographically smallest
    of the top or bottom strands.

    For positive overhang, the top strand is is left padded with the number
    of whitespace characters (ASCII 32) indicated by the overhang value.
    For negative overhang the reverse of the bottom strand is similarly padded.

    The two string are joined, separated by a line break (ASCII 10).

    The resulting string is used as an argument for the useguid function
    and the result is returned.
    ::

        dsDNA       overhang  lseguid_sticky

         TATGCC     1         Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU
         |||||
        catacg

         gcatac     1         Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU
         |||||
        CCGTAT


    For the linear dsDNA sequence defined by Watson = "TATGCC", Crick ="gcatac"
    and overhang = 1 (see figures above), The "gcatac" strand is selected as
    "gcatac" < "TATGCC".

    Overhang is positive, so the first strand is padded to " gcatac".

    A string is constructed like so:
    ::

        " gcatac" + chr(10) + "CCGTAT"


    and used as and argument for the useguid function.

    For overhang == 0, the same result is produced as would the lseguid_blunt
    function.


    Examples
    --------
    >>> lseguid_sticky("TATGCC", "gcatac", 1)
    'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
    >>> lseguid_sticky("gcatac", "TATGCC", 1)
    'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
    >>> useguid(" gcatac" + chr(10) + "CCGTAT")
    'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
    >>> lseguid_sticky("gcatac", "GTATGC", 0)
    'b0Xa5pLe4LNd5T8fhGWHicCI_f4'
    """
    return nseguid(watson=watson, crick=crick, overhang=overhang)


def cseguid(seq: str, srfunc=smallest_rotation) -> str:
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
    return nseguid(watson=seq, circular=True, srfunc=srfunc)


def nseguid(watson: str,
            crick: str = "",
            overhang: int = 0,
            circular: bool = False,
            ds: bool = True,
            srfunc: callable = smallest_rotation,
            rc: callable = rc,
            encoding: callable = base64.urlsafe_b64encode,
            cksumfunc: callable = seguid) -> str:
    """SEGUID checksums for linear or circular single or double stranded DNA.


    aaaaa

    Double stranded, circular DNA (dcSEGUID)
    ========================================
    Most bacterial plasmids are good examples of this kind of molecule.


    Single stranded, circular DNA  (scSEGUID)
    =========================================
    Some bacteriphages such as the M13 phage (Genbank NC_003287) is an example
    of this kind of molecule. This molecule has no complementary strand, so
    not a part of the checksum calculation. It does have n equally valid
    rotations where n is equal to the number of bases in the DNA string.
    The SEGUID for a single stranded, circular DNA (scSEGUID) is defined as
    the sha1 digest of the smallest rotation of the DNA string in uppercase.



    scSEGUID()



    CtfRtztP5kiC7OSqXWAIuiJfpok


    https://www.ncbi.nlm.nih.gov/nuccore/NC_003287.2

    Double stranded, linear DNA (dlSEGUID)
    ======================================

    Single stranded, linear DNA (slSEGUID)
    ======================================


    This function can calculate checksums or cSEGUID checksums depending on the
    arguments given.

    The crick strand is an optional argument, if not given it is inferred from
    the watson strand.

    If overhang is not 0, the crick argument is required as the complementary
    strand can not safely by deduced.

    If the optional boolean circular argument is set to False, the given
    sequence is interpreted as a linear double stranded DNA sequence.

    The lSEGUID checksum is calculated as described for the lSEGUID_sticky
    function.

    If the optional boolean circular argument is set to True, the given
    sequence is interpreted as a circular double stranded DNA sequence.
    The overhang is assumed to be zero for circular sequences.

    For a circular sequence nseguid is the uSEGUID checksum calculated for
    the lexicographically smallest of the minimum rotation of the watson and
    the minimum rotation of the crick strand. The optional fun argument sets
    a function for the caclulation of minimum rotation for a string. The fun
    argument is expected to take a string as an argument and return another
    string.

    The rc function is used by default to deduce the complementary sequence if
    needed. Another rc function can be used by setting the rc argument.

    The nseguid function uses the url-safe Base64 encoding by default, but
    this can be changed using the ´encoding´ argument.

    ttta
    aaat

    Examples
    --------
    >>> nseguid("TATGCC", "gcatac", 1)
    'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
    >>> nseguid("gcatac", "TATGCC", 1)
    'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
    >>> nseguid("TATGCC")
    'jZwQDKGSBpAF1wOh9icAv13hC2c'
    >>> nseguid("ggcata")
    'jZwQDKGSBpAF1wOh9icAv13hC2c'
    >>> nseguid("attt", circular=True)
    'oopV-6158nHJqedi8lsshIfcqYA'
    >>> nseguid("ttta", circular=True)
    'oopV-6158nHJqedi8lsshIfcqYA'
    >>> nseguid("ttta", ds=True)
    '3mcb7TbmWRazzfCKk5iohMo7REg'
    >>> nseguid("ttta", ds=False)
    '8zCvKwyQAEsbPtC4yTV-pY0H93Q'
    >>> nseguid("ttta", "TAAA", ds=False)  #
    Traceback (most recent call last):
        ...
    ValueError: ssDNA can only have a watson strand.
    >>> nseguid("ttta", overhang=1, circular=True)
    Traceback (most recent call last):
        ...
    ValueError: Circular sequence with overhang != 0.
    >>> nseguid("ttta", overhang=1)
    Traceback (most recent call last):
        ...
    ValueError: Overhang != 0 requires a crick strand.
    """
    if not ds and (watson and crick):
        raise ValueError("ssDNA can only have a watson strand.")
    if circular and overhang:
        raise ValueError("Circular sequence with overhang != 0.")
    if overhang and not crick:
        raise ValueError("Overhang != 0 requires a crick strand.")

    watson = watson.upper()

    w, c, o = "", "", 0

    if ds:

        crick = crick.upper() if crick else rc(watson)

        if circular:

            w = min(srfunc(watson), srfunc(crick))

        else:

            w, c, o = min(((watson, crick, overhang),
                           (crick, watson, len(watson) - len(crick) + overhang)))
    else:
        if circular:

            w = srfunc(watson)

        else:

            w = watson

    msg = f"{o*chr(32)}{w}{chr(10)}{-o*chr(32)}{c[::-1]}".rstrip()
    print(msg)
    return cksumfunc(msg,
                     encoding=encoding)


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
    print(timeit.timeit("cseguid(dna500, srfunc=smallest_rotation_py)",
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


    # for watson, crick in [("",""),("a",""),("","b"),("a","b")]:
    #     r = bool(watson and crick)
    #     print(r, watson, crick)
