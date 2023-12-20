#!/usr/bin/env bats

setup() {
    load "${BATS_SUPPORT_HOME:?}/load.bash"
    load "${BATS_ASSERT_HOME:?}/load.bash"
}

@test "Rscript -e seguid::seguid --args --version" {
    run Rscript -e seguid::seguid --args --version
    assert_success
}

@test "Rscript -e seguid::seguid --args --help" {
    run Rscript -e seguid::seguid --args --help
    assert_success
}

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
