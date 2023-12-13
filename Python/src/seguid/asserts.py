def assert_in_alphabet(seq: str,
                       alphabet: set):

    assert isinstance(seq, str), "Argument 'seq' must be an string"
    assert isinstance(alphabet, set), "Argument 'alphabet' must be a set"
    assert len(alphabet) > 0, "Argument 'alphabet' must not be empty"
    # breakpoint()
    # Nothing to do?
    if len(seq) == 0:
        return

    first = list(alphabet)[0]
    if isinstance(first, int):
        unknown = set(
            c for c in seq if c not in (chr(k) for k in alphabet)
        )
    elif isinstance(first, str):
        unknown = set(
            c for c in seq if c not in (k for k in alphabet)
        )
    else:
        raise ValueError("Unknown type of the elements in 'alphabet'")

    if unknown:
        missing = ' '.join(unknown)
        raise ValueError(
            "Detected symbols " f"{missing} in not in the 'alphabet'"
        )


def assert_table(table: dict):
    assert isinstance(table, dict), "Argument 'table' must be a dict"
    keys, values = table.items()

    ## Assert that the set of values are also in the set of keys
    unknown = set(
        chr(v) for v in values if v not in (k for k in keys)
    )

    if unknown:
        missing = ' '.join(unknown)
        raise ValueError(
            "Detected values (" f"{missing}) in 'table' that are not in the keys"
        )
