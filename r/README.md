[![R checks](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/check-r.yml/badge.svg)](https://github.com/MetabolicEngineeringGroupCBMA/seguid/actions/workflows/check-r.yml)

# SEGUID v2: Checksums Circular, Linear, Single- and Double-Stranded Sequences

This R package, **seguid**, implements SEGUID v2 together with the
original SEGUID algorithm.


## Example

```r
> library(seguid)

> lsseguid("AT")
[1] "lsseguid-Ax_RG6hzSrMEEWoCO1IWMGska-4"

> lsseguid("AT")
[1] "lsseguid-Ax_RG6hzSrMEEWoCO1IWMGska-4"

> csseguid("AT")
[1] "csseguid-Ax_RG6hzSrMEEWoCO1IWMGska-4"

> csseguid("TA")
[1] "csseguid-Ax_RG6hzSrMEEWoCO1IWMGska-4"

> cdseguid("AT", "AT")
[1] "cdseguid-AWD-dt5-TEua8RbOWfnctJIu9nA"

> cdseguid("TA", "TA")
[1] "cdseguid-AWD-dt5-TEua8RbOWfnctJIu9nA"
```


## Documentation

See `help(package = "seguid")` in R.


## Installation

To install the **seguid** R package, use:

```sh
$ remotes::install_github("MetabolicEngineeringGroupCBMA/seguid")
```
