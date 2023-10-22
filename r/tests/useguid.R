library(seguid)

x <- "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"
truth <- "cl5ukSUdlvZeBaBLEUhxisdRaL8"

## AD HOC: Because 'x' has both upper and lower-case symbols, which
## shouldn't really be allowed. The R package rejects them. /HB 2023-10-21
x <- toupper(x)

stopifnot(
  useguid(x) == truth,
  useguid(tolower(x)) == seguid(toupper(x)),

  useguid("aaa") == "YG7G6b2Kj_KtFOX63j8mRHHoIlE",
  useguid("aAA") == "YG7G6b2Kj_KtFOX63j8mRHHoIlE"
)


