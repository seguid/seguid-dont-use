from seguid.asserts import *
from seguid.tables import COMPLEMENT_TABLE

def rotate(seq: str, amount: int = 0) -> str:
    """Rotates a circular, DNA sequence a certain amount.

    Rotates sequence 'seq', 'amount' number of symbols to the right.
    A rotation 'amount' is the same as a rotation 'amount + n * len(seq)'
    for any integer 'n'.

    Returns the rotated sequence as a string of length 'len(seq)'.

    Examples
    --------
    >>> seq = "ABCDEFGH"
    >>> rotate(seq, 0)
    'ABCDEFGH'

    >>> rotate(seq, +1)
    'BCDEFGHA'

    >>> rotate(seq, +7)
    'HABCDEFG'

    >>> rotate(seq, -1)
    'HABCDEFG'

    >>> rotate(seq, +8)
    'ABCDEFGH'
    """
    assert isinstance(seq, str),    "Argument 'seq' must be an string"
    assert isinstance(amount, int), "Argument 'amount' must be an integer"

    ## Nothing to rotate?
    if len(seq) == 0:
        return seq

    amount = amount % len(seq)

    ## Rotate?
    if amount > 0:
        seq = seq[amount:] + seq[:amount]

    return seq


def complementary(seq: str, table: dict = COMPLEMENT_TABLE) -> str:
    """Complement of a DNA sequence.
    """
    ## Validate 'table':
    assert_table(table)

    ## Validate 'seq':
    assert_in_alphabet(seq, alphabet = set(table.keys()))

    return seq.translate(table)


def reverse(seq) -> str:
    """Reverses a DNA sequence
    """
    assert isinstance(seq, str),    "Argument 'seq' must be an string"
    return seq[::-1]


def rc(seq: str, table: dict = COMPLEMENT_TABLE) -> str:
    """Reverse complement of sequence.

    Returns the reverse complement for a DNA strand.

    The default complement table accepts GATC only.

    The tables module defines and alternative table containing the
    ambiguous codes suggested by IUPAC.

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
    return reverse(complementary(seq, table = table))
