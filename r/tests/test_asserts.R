assert_anneal <- seguid:::assert_anneal
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



tuples <- list(
   list("AT", "TA", 1),
   list("CTATAG", "AT", -2),
   list("AT", "CTATAG", 2),
   list("AT", "AT", 0)
)

alphabet <- c(G = "C", A = "T", T = "A", C = "G")
for (tuple in tuples) {
  watson <- tuple[[1]]
  crick <- tuple[[2]]
  overhang <- tuple[[3]]
  assert_anneal(watson, crick, overhang  = overhang, alphabet = alphabet)
}

tuples = list(
  list("AT", "CG", 1),
  list("CTATAG", "AT", -3),
  list("AT", "CTATAG", 1)
)

alphabet <- c(G = "C", A = "T", T = "A", C = "G")
for (tuple in tuples) {
  watson <- tuple[[1]]
  crick <- tuple[[2]]
  overhang <- tuple[[3]]
  res <- tryCatch({
    assert_anneal(watson, crick, overhang = overhang, alphabet = alphabet)
  }, error = identity)
  stopifnot(inherits(res, "error"))
}

res <- tryCatch({
  assert_anneal("AT", "AT", overhang = 4, alphabet = alphabet)
}, error = identity)
stopifnot(inherits(res, "error"))
