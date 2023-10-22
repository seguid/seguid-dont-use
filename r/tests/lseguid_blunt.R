library(seguid)

reverse_sequence <- seguid:::reverse_sequence
rc <- seguid:::rc
space <- seguid:::space

x <- "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"
truth <- "bHrqalTJ793oAigMQ5_qCttJRTk"
##truth <- "R7VPKJXvozX-xPk0wFeNxbZd_dM"

## AD HOC: Because 'x' has both upper and lower-case symbols, which
## shouldn't really be allowed. The R package rejects them. /HB 2023-10-21
x <- toupper(x)

stopifnot( 
  lseguid_blunt(x) == truth,
  lseguid_blunt(x) == useguid(paste0(c(rc(x), space(10L), reverse_sequence(x))))
)
