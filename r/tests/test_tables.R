get_table <- seguid:::get_table

table <- get_table("{DNA}")
table <- get_table("{RNA}")
table <- get_table("{IUPAC}")
table <- get_table("{protein}")

res <- tryCatch({
  table <- get_table("unknown")
}, error = identity)
stopifnot(inherits(res, "error"))
