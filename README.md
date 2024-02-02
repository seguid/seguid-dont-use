[![Python checks](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/check-python.yml/badge.svg)](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/check-python.yml)
[![R checks](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/check-r.yml/badge.svg)](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/check-r.yml)


# SEGUID v2: Checksums Circular, Linear, Single- and Double-Stranded Sequences

This repository hosts packaged implementations in [Python](python/)
and [R](r/) of SEGUID v2 together with the original SEGUID algorithm.

## Quick examples

Python:

```python
>>> from seguid import *

>>> lsseguid("AT")
'lsseguid-Ax_RG6hzSrMEEWoCO1IWMGska-4'

>>> cdseguid("AT", "AT")
'cdseguid-AWD-dt5-TEua8RbOWfnctJIu9nA'
```

```sh
$ python -m seguid --type="lsseguid" <<< "AT"
lsseguid-Ax_RG6hzSrMEEWoCO1IWMGska-4

$ python -m seguid --type="cdseguid" <<< $'AT\nTA'
cdseguid-AWD-dt5-TEua8RbOWfnctJIu9nA
```


R:

```r
> library(seguid)

> lsseguid("AT")
[1] "lsseguid-Ax_RG6hzSrMEEWoCO1IWMGska-4"

> cdseguid("AT", "AT")
[1] "cdseguid-AWD-dt5-TEua8RbOWfnctJIu9nA"
```

```sh
$ Rscript -e seguid::lsseguid <<< "AT"
lsseguid-Ax_RG6hzSrMEEWoCO1IWMGska-4

$ Rscript -e seguid::cdseguid <<< $'AT\nTA'
cdseguid-AWD-dt5-TEua8RbOWfnctJIu9nA
```
