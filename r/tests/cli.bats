#!/usr/bin/env bats

setup() {
    load "${BATS_SUPPORT_HOME:?}/load.bash"
    load "${BATS_ASSERT_HOME:?}/load.bash"
}


## --------------------------------------------------------
## Help and version
## --------------------------------------------------------
@test "Rscript -e seguid::seguid --args --version" {
    run Rscript -e seguid::seguid --args --version
    assert_success
    ## Assert numeric x.y.z format
    assert_output --regexp "\b[[:digit:]]+([.][[:digit:]]+)$"
}

@test "Rscript -e seguid::seguid --args --help" {
    run Rscript -e seguid::seguid --args --help
    assert_success
    assert_output --partial "seguid"
    assert_output --partial "--version"
    assert_output --partial "--help"
    assert_output --partial "Usage:"
    assert_output --partial "Options:"
    assert_output --partial "Examples:"
}



## --------------------------------------------------------
## Succesfully cases
## --------------------------------------------------------
@test "Rscript -e seguid::seguid <<< \"ACGT\"" {
    run Rscript -e seguid::seguid <<< "ACGT"
    assert_success
    assert_output "seguid:IQiZThf2zKn/I1KtqStlEdsHYDQ"
}


@test "Rscript -e seguid::seguid --type=seguid <<< \"ACGT\"" {
    run Rscript -e seguid::seguid --type=seguid <<< "ACGT"
    assert_success
    assert_output "seguid:IQiZThf2zKn/I1KtqStlEdsHYDQ"
}


@test "Rscript -e seguid::seguid --type=slseguid <<< \"ACGT\"" {
    run Rscript -e seguid::seguid --type=slseguid <<< "ACGT"
    assert_success
    assert_output "slseguid:IQiZThf2zKn_I1KtqStlEdsHYDQ"
}


@test "Rscript -e seguid::seguid --type=scseguid <<< \"ACGT\"" {
    run Rscript -e seguid::seguid --type=scseguid <<< "ACGT"
    assert_success
    assert_output "scseguid:IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::seguid --type=scseguid <<< \"CGTA\" (rotation invariant)" {
    run Rscript -e seguid::seguid --type=scseguid <<< "CGTA"
    assert_success
    assert_output "scseguid:IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::seguid --type=scseguid <<< \"GTAC\" (rotation invariant)" {
    run Rscript -e seguid::seguid --type=scseguid <<< "GTAC"
    assert_success
    assert_output "scseguid:IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::seguid --type=dlseguid <<< \$'ACGT\\nTGCA'" {
    run Rscript -e seguid::seguid --type=dlseguid <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dlseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::seguid --type=dcseguid <<< \$'ACGT\\nTGCA'" {
    run Rscript -e seguid::seguid --type=dcseguid <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::seguid --type=dcseguid <<< \$'CGTA\\nGCAT' (rotation invariant)" {
    run Rscript -e seguid::seguid --type=dcseguid <<< $'CGTA\nGCAT'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::seguid --type=dcseguid <<< \$'GTAC\\nCATG' (rotation invariant)" {
    run Rscript -e seguid::seguid --type=dcseguid <<< $'GTAC\nCATG'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::seguid --type=dlseguid <<< \$'ACGT\nTGCA'" {
    run Rscript -e seguid::seguid --type=dlseguid <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dlseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::seguid --type=dcseguid <<< \$'ACGT\nTGCA'" {
    run Rscript -e seguid::seguid --type=dcseguid <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::seguid --type=dcseguid <<< \$'CGTA\nGCAT' (rotation invariant)" {
    run Rscript -e seguid::seguid --type=dcseguid <<< $'CGTA\nGCAT'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::seguid --type=dlseguid <<< \$'-CGT\\nTGCA'" {
    run Rscript -e seguid::seguid --type=dlseguid <<< $'-CGT\nTGCA'
    assert_success
    assert_output "dlseguid:MpPe6pJoya3CoRh3BAw2qgEOcKI"
}

@test "Rscript -e seguid::seguid --type=dlseguid <<< \$'-CGT\nTGC-'" {
    run Rscript -e seguid::seguid --type=dlseguid <<< $'-CGT\nTGC-'
    assert_success
    assert_output "dlseguid:a_o4Ga8vQrhlvI_zkjUg0uu6obA"
}


## --------------------------------------------------------
## Corner cases (single-symbol input)
## --------------------------------------------------------
@test "Rscript -e seguid::seguid --type=seguid <<< \"\" (single-symbol input)" {
    run Rscript -e seguid::seguid --type=seguid <<< "A"
    assert_success
    assert_output "seguid:bc1M4j2I4u6VaLpUbAB8Y9kTHBs"
}

@test "Rscript -e seguid::seguid --type=slseguid <<< \"\" (single-symbol input)" {
    run Rscript -e seguid::seguid --type=slseguid <<< "A"
    assert_success
    assert_output "slseguid:bc1M4j2I4u6VaLpUbAB8Y9kTHBs"
}

@test "Rscript -e seguid::seguid --type=scseguid <<< \"\" (single-symbol input)" {
    run Rscript -e seguid::seguid --type=scseguid <<< "A"
    assert_success
    assert_output "scseguid:bc1M4j2I4u6VaLpUbAB8Y9kTHBs"
}

@test "Rscript -e seguid::seguid --type=dlseguid <<< \$'A\nT' (single-symbol input)" {
    run Rscript -e seguid::seguid --type=dlseguid <<< $'A\nT'
    assert_success
    assert_output "dlseguid:S4AfmFCoHYVrWNQ_d7-lVVF2t20"
}

@test "Rscript -e seguid::seguid --type=dcseguid <<< \$'A\nT' (single-symbol input)" {
    run Rscript -e seguid::seguid --type=dcseguid <<< $'A\nT'
    assert_success
    assert_output "dcseguid:S4AfmFCoHYVrWNQ_d7-lVVF2t20"
}


## --------------------------------------------------------
## Corner cases (empty input)
## --------------------------------------------------------
@test "Rscript -e seguid::seguid --type=seguid <<< \"\" (empty input)" {
    run Rscript -e seguid::seguid --type=seguid <<< ""
    assert_success
    assert_output "seguid:2jmj7l5rSw0yVb/vlWAYkK/YBwk"
}

@test "Rscript -e seguid::seguid --type=slseguid <<< \"\" (empty input)" {
    run Rscript -e seguid::seguid --type=slseguid <<< ""
    assert_success
    assert_output "slseguid:2jmj7l5rSw0yVb_vlWAYkK_YBwk"
}

@test "Rscript -e seguid::seguid --type=scseguid <<< \"\" (empty input)" {
    run Rscript -e seguid::seguid --type=scseguid <<< ""
    assert_success
    assert_output "scseguid:2jmj7l5rSw0yVb_vlWAYkK_YBwk"
}

@test "Rscript -e seguid::seguid --type=dlseguid <<< \"\" (empty input)" {
    skip "dlseguid() does not support empty input for now"
    run Rscript -e seguid::seguid --type=dlseguid <<< ""
    assert_success
}

@test "Rscript -e seguid::seguid --type=dcseguid <<< \"\" (empty input)" {
    skip "dcseguid() does not support empty input for now"
    run Rscript -e seguid::seguid --type=dcseguid <<< ""
    assert_success
}


## --------------------------------------------------------
## Failing cases
## --------------------------------------------------------
@test "Rscript -e seguid::seguid <<< \"aCGT\" gives error (invalid character)" {
    run Rscript -e seguid::seguid <<< "aCGT"
    assert_failure
}

@test "Rscript -e seguid::seguid <<< \" ACGT\" gives error (invalid character)" {
    run Rscript -e seguid::seguid <<< " ACGT"
    assert_failure
}

@test "Rscript -e seguid::seguid --type=dlseguid <<< \$' ACGT\\nTGCA' gives error (invalid character)" {
    run Rscript -e seguid::seguid --type=dlseguid <<< $' ACGT\nTGCA '
    assert_failure
}

@test "Rscript -e seguid::seguid --type=dlseguid <<< \$'ACGT\\nTGC' gives error (unbalanced lengths)" {
    run Rscript -e seguid::seguid --type=dlseguid <<< $'ACGT\nTGC'
    assert_failure
}

@test "Rscript -e seguid::seguid --type=dlseguid <<< \$'ACGT\\nTGCC' gives error (incompatible sequences)" {
    run Rscript -e seguid::seguid --type=dlseguid <<< $'ACGT\nTGCC'
    assert_failure
}


## --------------------------------------------------------
## The follow CLI calls are specific to R
## --------------------------------------------------------
@test "Rscript -e seguid::slseguid <<< \"ACGT\"" {
    run Rscript -e seguid::slseguid <<< "ACGT"
    assert_success
    assert_output "slseguid:IQiZThf2zKn_I1KtqStlEdsHYDQ"
}


@test "Rscript -e seguid::scseguid <<< \"ACGT\"" {
    run Rscript -e seguid::scseguid <<< "ACGT"
    assert_success
    assert_output "scseguid:IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::scseguid <<< \"CGTA\" (rotation invariant)" {
    run Rscript -e seguid::scseguid <<< "CGTA"
    assert_success
    assert_output "scseguid:IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::scseguid <<< \"GTAC\" (rotation invariant)" {
    run Rscript -e seguid::scseguid <<< "GTAC"
    assert_success
    assert_output "scseguid:IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::dlseguid <<< \$'ACGT\\nTGCA'" {
    run Rscript -e seguid::dlseguid <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dlseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid <<< \$'ACGT\\nTGCA'" {
    run Rscript -e seguid::dcseguid <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid <<< \$'CGTA\\nGCAT' (rotation invariant)" {
    run Rscript -e seguid::dcseguid <<< $'CGTA\nGCAT'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid <<< \$'GTAC\\nCATG' (rotation invariant)" {
    run Rscript -e seguid::dcseguid <<< $'GTAC\nCATG'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dlseguid <<< \$'ACGT\nTGCA'" {
    run Rscript -e seguid::dlseguid <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dlseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid <<< \$'ACGT\nTGCA'" {
    run Rscript -e seguid::dcseguid <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid <<< \$'CGTA\nGCAT' (rotation invariant)" {
    run Rscript -e seguid::dcseguid <<< $'CGTA\nGCAT'
    assert_success
    assert_output "dcseguid:ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dlseguid <<< \$'-CGT\\nTGCA'" {
    run Rscript -e seguid::dlseguid <<< $'-CGT\nTGCA'
    assert_success
    assert_output "dlseguid:MpPe6pJoya3CoRh3BAw2qgEOcKI"
}

@test "Rscript -e seguid::dlseguid <<< \$'-CGT\nTGC-'" {
    run Rscript -e seguid::dlseguid <<< $'-CGT\nTGC-'
    assert_success
    assert_output "dlseguid:a_o4Ga8vQrhlvI_zkjUg0uu6obA"
}
