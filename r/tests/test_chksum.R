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


stopifnot(  seguid("AT") ==   "seguid-Ax/RG6hzSrMEEWoCO1IWMGska+4")
stopifnot(slseguid("AT") == "slseguid-Ax_RG6hzSrMEEWoCO1IWMGska-4")

NP_313053_1 <- paste0(
  "MKALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSG",
  "ASRGIRLLQEEEEGLPLVGRVAAGEPLLAQQHIEGHYQVDPSLFKPNADFLLRVSGMSMKD",
  "IGIMDGDLLAVHKTQDVRNGQVVVARIDDEVTVKRLKKQGNKVELLPENSEFKPIVVDLRQ",
  "QSFTIEGLAVGVIRNGDWL"
)
#stopifnot(seguid(NP_313053_1, table="{IUPAC}") == "seguid-2c4yjE+JqjvzYF1d0OmUh8pCpz8")


m13dna <- readLines("test_data/M13.txt")
truth <- "scseguid-aAjgnsF9cPI6cu8IQ81sYnstVzU"
stopifnot(scseguid(m13dna) == truth)


dlDNA <- "AT"
truth <- "dlseguid-AWD-dt5-TEua8RbOWfnctJIu9nA"
stopifnot(dlseguid(dlDNA, rc(dlDNA), 0) == truth)

dlDNAb <- list("AT", rc("AT"), 0)
stopifnot(do.call(dlseguid, args = dlDNAb) == truth)

dlDNA2 <- list("AT", "TA", 1)
truth <- "dlseguid-JwB2eUmZkCNjyWAv471JeUbiSDM"
stopifnot(do.call(dlseguid, args = dlDNA2) == truth)

dlDNA3 <- list("TA", "AT", -1)
truth <- "dlseguid-bv0UOR12eWrBeaAx79PNZvveviU"
stopifnot(do.call(dlseguid, args = dlDNA3) == truth)

dlDNA4 <- list("CTATAG", "AT", -2)
truth <- "dlseguid-np3hncfQvOh8rZ8Co1Ts_02NXg4"
stopifnot(do.call(dlseguid, args = dlDNA4) == truth)

dlDNA5 <- list("AT", "CTATAG", 2)
truth <- "dlseguid-np3hncfQvOh8rZ8Co1Ts_02NXg4"
stopifnot(do.call(dlseguid, args = dlDNA5) == truth)

truth <- "dcseguid-tYeHZYwxQGDHTqGDcrebERag0AU"
stopifnot(dcseguid("ACGTT", "AACGT") == truth)
stopifnot(dcseguid("AACGT", "ACGTT") == truth)

pUC19dna <- readLines("test_data/pUC19.txt", warn = FALSE)
truth <- "dcseguid-zhw8Yrxfo3FO5DDccx4PamBVPCQ"
stopifnot(dcseguid(pUC19dna, rc(pUC19dna)) == truth)
bfr <- readLines("test_data/pUC19_minimal_rotation_watson_linebreak_crick.txt", warn = FALSE)
w <- bfr[1]
c <- bfr[2]
stopifnot(dlseguid(w, reverse(c), 0) == gsub("dcseguid-", "dlseguid-", truth))


## Empty input is considered an error
assert_error(seguid::seguid(""))
assert_error(seguid::slseguid(""))
assert_error(seguid::scseguid(""))
assert_error(seguid::dlseguid("", "", overhang = 0))
assert_error(seguid::dcseguid("", ""))
