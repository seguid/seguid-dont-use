rc <- seguid:::rc
complementary <- seguid:::complementary
rotate <- seguid:::rotate
min_rotation <- seguid:::min_rotation
rotate_to_min <- seguid:::rotate_to_min

table <- c(G = "C", A = "T", T = "A", C = "G")

seq <- "AACCGGTT"
seq_complementary <- "TTGGCCAA"
stopifnot(complementary(seq, table = table) == seq_complementary)
stopifnot(complementary(seq_complementary, table = table) == seq)
stopifnot(complementary(complementary(seq, table = table), table = table) == seq)

seq <- "AACCGGTTxx"
res <- tryCatch({
  complementary(seq, table = table)
}, error = identity)
stopifnot(inherits(res, "error"))

watson <- "ACGTAACCGGTT"
crick <- "AACCGGTTACGT"
stopifnot(rc(watson, table = table) == crick)
stopifnot(rc(crick, table = table) == watson)
stopifnot(rc(rc(watson, table = table), table = table) == watson)

seq <- "ACGTAACCGGTT"
n <- nchar(seq)
stopifnot(rotate(seq,   0) ==        seq     )
stopifnot(rotate(seq,   n) == rotate(seq,  0))
stopifnot(rotate(seq, 2*n) == rotate(seq,  0))
stopifnot(rotate(seq, n-1) == rotate(seq, -1))

stopifnot(rotate("", 0) == "")
stopifnot(rotate("", 1) == "")

## Rotate on the complementary strand
stopifnot(complementary(rotate(complementary(seq, table = table), +1), table = table) == rotate(seq, +1))

watson <- "ACGTAACCGGTT"
crick <- "AACCGGTTACGT"

## Rotate on Watson, is the opposite rotation on Crick
stopifnot(rc(watson, table = table)                 == crick             )
stopifnot(rc(rotate(crick, -1), table = table)      == rotate(watson, +1))
stopifnot(rc(rotate(rc(watson, table = table), -1), table = table) == rotate(watson, +1))

stopifnot(rc("GAT", table = table) == "ATC")
stopifnot(rc("GTT", table = table) == "AAC")

res <- tryCatch({
  rc("GTZ", table = table)
}, error = identity)
stopifnot(inherits(res, "error"))


stopifnot(
  min_rotation("Aa") == 0
)

seq <- "TAAA"
amount <- min_rotation(seq)
stopifnot(amount == 1)
stopifnot(rotate(seq, amount = amount) == "AAAT")

seq <- "ACAACAAACAACACAAACAAACACAAC"
amount <- min_rotation(seq)
stopifnot(amount == 14)
stopifnot(rotate(seq, amount = amount) == "AAACAAACACAACACAACAAACAACAC")


stopifnot(rotate_to_min("taaa") == "aaat")
stopifnot(rotate_to_min("abaabaaabaababaaabaaababaab") == "aaabaaababaababaabaaabaabab")
stopifnot(rotate_to_min("abaabaaabaababaaabaaaBabaab") == "Babaababaabaaabaababaaabaaa")
stopifnot(rotate_to_min("TAAA") == "AAAT")
stopifnot(rotate_to_min("ACAACAAACAACACAAACAAACACAAC") == "AAACAAACACAACACAACAAACAACAC")
