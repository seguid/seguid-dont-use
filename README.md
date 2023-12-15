[![Python checks](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/check-python.yml/badge.svg)](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/check-python.yml)
[![Python Code Coverage](https://codecov.io/gh/MetabolicEngineeringGroupCBMA/seguid/graph/badge.svg)](https://codecov.io/gh/MetabolicEngineeringGroupCBMA/seguid)
[![R checks](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/R-CMD-check.yaml)


# seguid

Functions for calculating uSEGUID, lSEGUID and cSEGUID checksums for
biological sequences.


## R


## C++



## Rust



## Python

[![Documentation Status](https://readthedocs.org/projects/seguid/badge/?version=latest)](https://seguid.readthedocs.io/en/latest/?badge=latest)

Install with:

```sh
$ python -m pip install --user seguid
```

use:

```python
Python 3.8.18 | packaged by conda-forge | (default, Oct 10 2023, 15:44:36)
[GCC 12.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from seguid import dsseguid
>>> dsseguid("TATGCC", "gcatac", 1)
'Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU'
>>>
```

Run tests with pytest without arguments in the Python directory;

```sh
$ cd Python/
$ pytest
```



