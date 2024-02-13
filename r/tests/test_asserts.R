assert_valid_alphabet <- seguid:::assert_valid_alphabet
assert_in_alphabet <- seguid:::assert_in_alphabet
assert_alphabet <- seguid:::assert_alphabet

seq <- "ABCDEFGH"
alphabet <- c("A", "C", "G", "T")

assert_valid_alphabet(alphabet)
res <- tryCatch({
  assert_valid_alphabet("!")
}, error = identity)
stopifnot(inherits(res, "error"))

assert_in_alphabet("ACGT", alphabet = alphabet)
assert_in_alphabet("AAAA", alphabet = alphabet)
assert_in_alphabet("", alphabet = alphabet)

res <- tryCatch({
  assert_in_alphabet("x", alphabet = alphabet)
}, error = identity)
stopifnot(inherits(res, "error"))

res <- tryCatch({
  assert_in_alphabet("ACGTx", alphabet = alphabet)
}, error = identity)
stopifnot(inherits(res, "error"))



alphabet <- c(G = "C", A = "T", T = "A", C = "G")
assert_alphabet(alphabet)

res <- tryCatch({
  alphabet <- c(G = "C", A = "T", T = "A", x = "G")
  assert_alphabet(alphabet)
}, error = identity)
stopifnot(inherits(res, "error"))
