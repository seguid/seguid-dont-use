library(seguid)

x <- "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"
truth <- "cl5ukSUdlvZeBaBLEUhxisdRaL8"

## AD HOC: Because 'x' has both upper and lower-case symbols, which
## shouldn't really be allowed. The R package rejects them. /HB 2023-10-21
x <- toupper(x)

stopifnot(
  seguid(x) == truth,
  seguid(tolower(x)) == seguid(toupper(x)),
  
  seguid("ACGTACGTACGT") == "If6HIvcnRSQDVNiAoefAzySc6i4",
  seguid("acgtACGTacgt") == "If6HIvcnRSQDVNiAoefAzySc6i4"
)
