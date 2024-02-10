#!/usr/bin/env bats

setup() {
    load "${BATS_SUPPORT_HOME:?}/load.bash"
    load "${BATS_ASSERT_HOME:?}/load.bash"
}


## --------------------------------------------------------
## The follow CLI calls are specific to R
## --------------------------------------------------------
@test "Rscript -e seguid::lsseguid --args <<< \"ACGT\"" {
    run Rscript -e seguid::lsseguid --args <<< "ACGT"
    assert_success
    assert_output "lsseguid=IQiZThf2zKn_I1KtqStlEdsHYDQ"
}


@test "Rscript -e seguid::csseguid --args <<< \"ACGT\"" {
    run Rscript -e seguid::csseguid --args <<< "ACGT"
    assert_success
    assert_output "csseguid=IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::csseguid --args <<< \"CGTA\" (rotation invariant)" {
    run Rscript -e seguid::csseguid --args <<< "CGTA"
    assert_success
    assert_output "csseguid=IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::csseguid --args <<< \"GTAC\" (rotation invariant)" {
    run Rscript -e seguid::csseguid --args <<< "GTAC"
    assert_success
    assert_output "csseguid=IQiZThf2zKn_I1KtqStlEdsHYDQ"
}

@test "Rscript -e seguid::ldseguid --args <<< \$'ACGT\\nTGCA'" {
    run Rscript -e seguid::ldseguid --args <<< $'ACGT\nTGCA'
    assert_success
    assert_output "ldseguid=ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::cdseguid --args <<< \$'ACGT\\nTGCA'" {
    run Rscript -e seguid::cdseguid --args <<< $'ACGT\nTGCA'
    assert_success
    assert_output "cdseguid=ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::cdseguid --args <<< \$'CGTA\\nGCAT' (rotation invariant)" {
    run Rscript -e seguid::cdseguid --args <<< $'CGTA\nGCAT'
    assert_success
    assert_output "cdseguid=ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::cdseguid --args <<< \$'GTAC\\nCATG' (rotation invariant)" {
    run Rscript -e seguid::cdseguid --args <<< $'GTAC\nCATG'
    assert_success
    assert_output "cdseguid=ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::ldseguid --args <<< \$'ACGT\nTGCA'" {
    run Rscript -e seguid::ldseguid --args <<< $'ACGT\nTGCA'
    assert_success
    assert_output "ldseguid=ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::cdseguid --args <<< \$'ACGT\nTGCA'" {
    run Rscript -e seguid::cdseguid --args <<< $'ACGT\nTGCA'
    assert_success
    assert_output "cdseguid=ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::cdseguid --args <<< \$'CGTA\nGCAT' (rotation invariant)" {
    run Rscript -e seguid::cdseguid --args <<< $'CGTA\nGCAT'
    assert_success
    assert_output "cdseguid=ZubWOQ_QWYz_5fy9qVNfCGXZhag"
}

@test "Rscript -e seguid::ldseguid --args <<< \$'-CGT\\nTGCA'" {
    run Rscript -e seguid::ldseguid --args <<< $'-CGT\nTGCA'
    assert_success
    assert_output "ldseguid=ONPHQCrPDPDbypL85mg8vXNQGPw"
}

@test "Rscript -e seguid::ldseguid --args <<< \$'-CGT\nTGC-'" {
    run Rscript -e seguid::ldseguid --args <<< $'-CGT\nTGC-'
    assert_success
    assert_output "ldseguid=a_o4Ga8vQrhlvI_zkjUg0uu6obA"
}
