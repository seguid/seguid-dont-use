library(seguid)

TABLE_IUPAC_PROTEIN <- seguid:::TABLE_IUPAC_PROTEIN
rc <- seguid:::rc
reverse <- seguid:::reverse

stopifnot(  seguid("AT") ==   "seguid:Ax/RG6hzSrMEEWoCO1IWMGska+4")
stopifnot(slseguid("AT") == "slseguid:Ax_RG6hzSrMEEWoCO1IWMGska-4")

NP_313053_1 <- paste0(
  "MKALTARQQEVFDLIRDHISQTGMPPTRAEIAQRLGFRSPNAAEEHLKALARKGVIEIVSG",
  "ASRGIRLLQEEEEGLPLVGRVAAGEPLLAQQHIEGHYQVDPSLFKPNADFLLRVSGMSMKD",
  "IGIMDGDLLAVHKTQDVRNGQVVVARIDDEVTVKRLKKQGNKVELLPENSEFKPIVVDLRQ",
  "QSFTIEGLAVGVIRNGDWL"
)
#stopifnot(seguid(NP_313053_1, table=TABLE_IUPAC_PROTEIN) == "seguid:2c4yjE+JqjvzYF1d0OmUh8pCpz8")


m13dna <- readLines("test_data/M13.txt")
truth <- "scseguid:aAjgnsF9cPI6cu8IQ81sYnstVzU"
stopifnot(scseguid(m13dna) == truth)


dlDNA <- "AT"
truth <- "dlseguid:AWD-dt5-TEua8RbOWfnctJIu9nA"
stopifnot(dlseguid(dlDNA, rc(dlDNA), 0) == truth)

dlDNAb <- list("AT", rc("AT"), 0)
stopifnot(do.call(dlseguid, args = dlDNAb) == truth)

dlDNA2 <- list("AT", "TA", 1)
truth <- "dlseguid:JwB2eUmZkCNjyWAv471JeUbiSDM"
stopifnot(do.call(dlseguid, args = dlDNA2) == truth)

dlDNA3 <- list("TA", "AT", -1)
truth <- "dlseguid:bv0UOR12eWrBeaAx79PNZvveviU"
stopifnot(do.call(dlseguid, args = dlDNA3) == truth)

dlDNA4 <- list("CTATAG", "AT", -2)
truth <- "dlseguid:np3hncfQvOh8rZ8Co1Ts_02NXg4"
stopifnot(do.call(dlseguid, args = dlDNA4) == truth)

dlDNA5 <- list("AT", "CTATAG", 2)
truth <- "dlseguid:np3hncfQvOh8rZ8Co1Ts_02NXg4"
stopifnot(do.call(dlseguid, args = dlDNA5) == truth)

pUC19dna <- readLines("test_data/puc19.txt", warn = FALSE)
truth <- "dcseguid:zhw8Yrxfo3FO5DDccx4PamBVPCQ"
stopifnot(dcseguid(pUC19dna, rc(pUC19dna)) == truth)
bfr <- readLines("test_data/pUC19msg.txt", warn = FALSE)
w <- bfr[1]
c <- bfr[2]
stopifnot(dlseguid(w, reverse(c), 0) == gsub("dcseguid:", "dlseguid:", truth))
