from asserts import *
from tables import *

def rotate(seq: str, amount: int = 0) -> str:
    """Rotates a circular, genomic sequence a certain amount.

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
    """Complement of a genomic sequence.
    """
    ## Validate 'table':
    assert_table(table)
    
    ## Validate 'seq':
    assert_in_alphabet(seq, alphabet = set(table.keys()))
    
    return seq.translate(table)


def reverse(seq) -> str:
    """Reverses a genomic sequence
    """
    assert isinstance(seq, str),    "Argument 'seq' must be an string"
    return seq[::-1]


def rc(seq: str, table: dict = COMPLEMENT_TABLE) -> str:
    """Reverse complement of a genomic sequence.
    """
    return reverse(complementary(seq, table = table))
