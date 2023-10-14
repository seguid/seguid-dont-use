library(seguid)

x <- "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"
stopifnot(lseguid_blunt(x) == "bHrqalTJ793oAigMQ5_qCttJRTk")
