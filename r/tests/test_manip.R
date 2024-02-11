rc <- seguid:::rc
complementary <- seguid:::complementary
rotate <- seguid:::rotate
min_rotation <- seguid:::min_rotation
rotate_to_min <- seguid:::rotate_to_min

alphabet <- c(G = "C", A = "T", T = "A", C = "G")

seq <- "AACCGGTT"
seq_complementary <- "TTGGCCAA"
stopifnot(complementary(seq, alphabet = alphabet) == seq_complementary)
stopifnot(complementary(seq_complementary, alphabet = alphabet) == seq)
stopifnot(complementary(complementary(seq, alphabet = alphabet), alphabet = alphabet) == seq)

seq <- "AACCGGTTxx"
res <- tryCatch({
  complementary(seq, alphabet = alphabet)
}, error = identity)
stopifnot(inherits(res, "error"))

watson <- "ACGTAACCGGTT"
crick <- "AACCGGTTACGT"
stopifnot(rc(watson, alphabet = alphabet) == crick)
stopifnot(rc(crick, alphabet = alphabet) == watson)
stopifnot(rc(rc(watson, alphabet = alphabet), alphabet = alphabet) == watson)

seq <- "ACGTAACCGGTT"
n <- nchar(seq)
stopifnot(rotate(seq,   0) ==        seq     )
stopifnot(rotate(seq,   n) == rotate(seq,  0))
stopifnot(rotate(seq, 2*n) == rotate(seq,  0))
stopifnot(rotate(seq, n-1) == rotate(seq, -1))

stopifnot(rotate("", 0) == "")
stopifnot(rotate("", 1) == "")

## Rotate on the complementary strand
stopifnot(complementary(rotate(complementary(seq, alphabet = alphabet), +1), alphabet = alphabet) == rotate(seq, +1))

watson <- "ACGTAACCGGTT"
crick <- "AACCGGTTACGT"

## Rotate on Watson, is the opposite rotation on Crick
stopifnot(rc(watson, alphabet = alphabet)                 == crick             )
stopifnot(rc(rotate(crick, -1), alphabet = alphabet)      == rotate(watson, +1))
stopifnot(rc(rotate(rc(watson, alphabet = alphabet), -1), alphabet = alphabet) == rotate(watson, +1))

stopifnot(rc("GAT", alphabet = alphabet) == "ATC")
stopifnot(rc("GTT", alphabet = alphabet) == "AAC")

res <- tryCatch({
  rc("GTZ", alphabet = alphabet)
}, error = identity)
stopifnot(inherits(res, "error"))


stopifnot(
  min_rotation("Aa") == 0
)

stopifnot(
  min_rotation("A") == 0
)

stopifnot(
  min_rotation("") == 0
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
