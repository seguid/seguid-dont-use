#!/usr/bin/env bats

setup() {
    load "${BATS_SUPPORT_HOME:?}/load.bash"
    load "${BATS_ASSERT_HOME:?}/load.bash"
}


## --------------------------------------------------------
## The follow CLI calls are specific to R
## --------------------------------------------------------
@test "Rscript -e seguid::slseguid --args <<< \"ACGT\"" {
    run Rscript -e seguid::slseguid --args <<< "ACGT"
    assert_success
    assert_output "slseguid-IQiZThf2zKn_I1KtqStlEdsHYDQ"
}


@test "Rscript -e seguid::scseguid --args <<< \"ACGT\"" {
    run Rscript -e seguid::scseguid --args <<< "ACGT"
    assert_success
    assert_output "scseguid-IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::scseguid --args <<< \"CGTA\" (rotation invariant)" {
    run Rscript -e seguid::scseguid --args <<< "CGTA"
    assert_success
    assert_output "scseguid-IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::scseguid --args <<< \"GTAC\" (rotation invariant)" {
    run Rscript -e seguid::scseguid --args <<< "GTAC"
    assert_success
    assert_output "scseguid-IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::dlseguid --args <<< \$'ACGT\\nTGCA'" {
    run Rscript -e seguid::dlseguid --args <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dlseguid-ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid --args <<< \$'ACGT\\nTGCA'" {
    run Rscript -e seguid::dcseguid --args <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dcseguid-ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid --args <<< \$'CGTA\\nGCAT' (rotation invariant)" {
    run Rscript -e seguid::dcseguid --args <<< $'CGTA\nGCAT'
    assert_success
    assert_output "dcseguid-ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid --args <<< \$'GTAC\\nCATG' (rotation invariant)" {
    run Rscript -e seguid::dcseguid --args <<< $'GTAC\nCATG'
    assert_success
    assert_output "dcseguid-ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dlseguid --args <<< \$'ACGT\nTGCA'" {
    run Rscript -e seguid::dlseguid --args <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dlseguid-ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid --args <<< \$'ACGT\nTGCA'" {
    run Rscript -e seguid::dcseguid --args <<< $'ACGT\nTGCA'
    assert_success
    assert_output "dcseguid-ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dcseguid --args <<< \$'CGTA\nGCAT' (rotation invariant)" {
    run Rscript -e seguid::dcseguid --args <<< $'CGTA\nGCAT'
    assert_success
    assert_output "dcseguid-ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::dlseguid --args <<< \$'-CGT\\nTGCA'" {
    run Rscript -e seguid::dlseguid --args <<< $'-CGT\nTGCA'
    assert_success
    assert_output "dlseguid-MpPe6pJoya3CoRh3BAw2qgEOcKI"
}

@test "Rscript -e seguid::dlseguid --args <<< \$'-CGT\nTGC-'" {
    run Rscript -e seguid::dlseguid --args <<< $'-CGT\nTGC-'
    assert_success
    assert_output "dlseguid-a_o4Ga8vQrhlvI_zkjUg0uu6obA"
}
