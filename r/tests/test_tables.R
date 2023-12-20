get_table <- seguid:::get_table

table <- get_table("dna")
table <- get_table("rna")
table <- get_table("iupac")
table <- get_table("protein")

res <- tryCatch({
  table <- get_table("unknown")
}, error = identity)
stopifnot(inherits(res, "error"))
