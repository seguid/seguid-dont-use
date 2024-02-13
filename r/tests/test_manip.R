rotate <- seguid:::rotate
min_rotation <- seguid:::min_rotation
rotate_to_min <- seguid:::rotate_to_min

alphabet <- c(G = "C", A = "T", T = "A", C = "G")

seq <- "ACGTAACCGGTT"
n <- nchar(seq)
stopifnot(rotate(seq,   0) ==        seq     )
stopifnot(rotate(seq,   n) == rotate(seq,  0))
stopifnot(rotate(seq, 2*n) == rotate(seq,  0))
stopifnot(rotate(seq, n-1) == rotate(seq, -1))

stopifnot(rotate("", 0) == "")
stopifnot(rotate("", 1) == "")

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
