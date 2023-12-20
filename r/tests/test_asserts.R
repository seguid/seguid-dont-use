assert_anneal <- seguid:::assert_anneal
assert_alphabet <- seguid:::assert_alphabet
assert_in_alphabet <- seguid:::assert_in_alphabet
assert_table <- seguid:::assert_table

seq <- "ABCDEFGH"
alphabet <- c("A", "C", "G", "T")

assert_alphabet(alphabet)
res <- tryCatch({
  assert_alphabet("!")
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



table <- c(G = "C", A = "T", T = "A", C = "G")
assert_table(table)

res <- tryCatch({
  table <- c(G = "C", A = "T", T = "A", x = "G")
  assert_table(table)
}, error = identity)
stopifnot(inherits(res, "error"))



tuples <- list(
   list("AT", "TA", 1),
   list("CTATAG", "AT", -2),
   list("AT", "CTATAG", 2),
   list("AT", "AT", 0)
)

table <- c(G = "C", A = "T", T = "A", C = "G")
for (tuple in tuples) {
  watson <- tuple[[1]]
  crick <- tuple[[2]]
  overhang <- tuple[[3]]
  assert_anneal(watson, crick, overhang  = overhang, table = table)
}

tuples = list(
  list("AT", "CG", 1),
  list("CTATAG", "AT", -3),
  list("AT", "CTATAG", 1)
)

table <- c(G = "C", A = "T", T = "A", C = "G")
for (tuple in tuples) {
  watson <- tuple[[1]]
  crick <- tuple[[2]]
  overhang <- tuple[[3]]
  res <- tryCatch({
    assert_anneal(watson, crick, overhang = overhang, table = table)
  }, error = identity)
  stopifnot(inherits(res, "error"))
}

res <- tryCatch({
  assert_anneal("AT", "AT", overhang = 4, table = table)
}, error = identity)
stopifnot(inherits(res, "error"))
