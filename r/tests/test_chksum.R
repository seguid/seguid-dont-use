library(seguid)

rc <- seguid:::rc
reverse <- seguid:::reverse

assert_error <- function(expr, envir = parent.frame()) {
  expr <- substitute(expr)
  res <- tryCatch(eval(expr, envir = envir), error = identity)
  if (!inherits(res, "error")) {
    stop("Call did not result in an error: ", deparse(expr))
  }
}


stopifnot(  seguid("AT") ==   "seguid=Ax/RG6hzSrMEEWoCO1IWMGska+4")
stopifnot(lsseguid("AT") == "lsseguid=Ax_RG6hzSrMEEWoCO1IWMGska-4")

NP_313053_1 <- paste0(
  "MKALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSG",
  "ASRGIRLLQEEEEGLPLVGRVAAGEPLLAQQHIEGHYQVDPSLFKPNADFLLRVSGMSMKD",
  "IGIMDGDLLAVHKTQDVRNGQVVVARIDDEVTVKRLKKQGNKVELLPENSEFKPIVVDLRQ",
  "QSFTIEGLAVGVIRNGDWL"
)
#stopifnot(seguid(NP_313053_1, alphabet="{IUPAC}") == "seguid=2c4yjE+JqjvzYF1d0OmUh8pCpz8")


m13dna <- readLines("test_data/M13.txt")
truth <- "csseguid=aAjgnsF9cPI6cu8IQ81sYnstVzU"
stopifnot(csseguid(m13dna) == truth)


dlDNA <- "AT"
truth <- "ldseguid=AWD-dt5-TEua8RbOWfnctJIu9nA"
stopifnot(ldseguid(dlDNA, rc(dlDNA), 0) == truth)

dlDNAb <- list("AT", rc("AT"), 0)
stopifnot(do.call(ldseguid, args = dlDNAb) == truth)

dlDNA2 <- list("AT", "TA", 1)
truth <- "ldseguid=JwB2eUmZkCNjyWAv471JeUbiSDM"
stopifnot(do.call(ldseguid, args = dlDNA2) == truth)

dlDNA3 <- list("TA", "AT", -1)
truth <- "ldseguid=XBcVadfQevTW_lklW4rdqw5udQ8"
stopifnot(do.call(ldseguid, args = dlDNA3) == truth)

dlDNA4 <- list("CTATAG", "AT", -2)
truth <- "ldseguid=_E05Xeo7KnLxrjsqDdpXNw_AIDE"
stopifnot(do.call(ldseguid, args = dlDNA4) == truth)

dlDNA5 <- list("AT", "CTATAG", 2)
truth <- "ldseguid=np3hncfQvOh8rZ8Co1Ts_02NXg4"
stopifnot(do.call(ldseguid, args = dlDNA5) == truth)

truth <- "cdseguid=tYeHZYwxQGDHTqGDcrebERag0AU"
stopifnot(cdseguid("ACGTT", "AACGT") == truth)
stopifnot(cdseguid("AACGT", "ACGTT") == truth)

pUC19dna <- readLines("test_data/pUC19.txt", warn = FALSE)
truth <- "cdseguid=zhw8Yrxfo3FO5DDccx4PamBVPCQ"
stopifnot(cdseguid(pUC19dna, rc(pUC19dna)) == truth)
bfr <- readLines("test_data/pUC19_minimal_rotation_watson_linebreak_crick.txt", warn = FALSE)
w <- bfr[1]
c <- bfr[2]
stopifnot(ldseguid(w, reverse(c), 0) == gsub("cdseguid=", "ldseguid=", truth))


## Empty input is considered an error
assert_error(seguid::seguid(""))
assert_error(seguid::lsseguid(""))
assert_error(seguid::csseguid(""))
assert_error(seguid::ldseguid("", "", overhang = 0))
assert_error(seguid::cdseguid("", ""))



## Use checksums as filenames
seq <- "GATTACA"
## Comment:
## The   SEGUID check is seguid=tp2jzeCM2e3W4yxtrrx09CMKa/8
## The slSEGUID check is seguid=tp2jzeCM2e3W4yxtrrx09CMKa_8
td <- tempdir()
filename <- seguid::lsseguid(seq)
pathname <- file.path(td, filename)
cat(seq, file = pathname)
stopifnot(utils::file_test("-f", pathname))
content <- readLines(pathname, warn = FALSE)
stopifnot(identical(content, seq))
file.remove(pathname)
unlink(td)
