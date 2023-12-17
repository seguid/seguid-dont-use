assert_in_alphabet <- function(seq, alphabet) {
  stopifnot(length(seq) == 1, is.character(seq), !is.na(seq))
  stopifnot(length(alphabet) > 0, is.character(alphabet), !anyNA(alphabet))

  ## Nothing to do?
  if (nchar(seq) == 0) {
    return(seq)
  }

  seq <- strsplit(seq, split = "", fixed = TRUE)[[1]]
  unknown <- setdiff(seq, alphabet)
  if (length(unknown) > 0) {
    missing <- paste(unknown, collapse = " ")
    stop(sprintf("Detected symbols %s in not in the 'alphabet'", missing))
  }
}

assert_table <- function(table) {
  stopifnot(is.character(table), !anyNA(table), is.character(names(table)))
  unknown <- setdiff(table, names(table))
  if (length(unknown) > 0) {
    missing <- paste(unknown, collapse = " ")
    stop(sprintf("Detected values (%s) in 'table' that are not in the names", missing))
  }
}

assert_anneal <- function(watson, crick, overhang, table = COMPLEMENT_TABLE_DNA) {
  assert_table(table)
  assert_in_alphabet(watson, alphabet = names(table))
  assert_in_alphabet(crick, alphabet = names(table))
  stopifnot(length(overhang) == 1, is.numeric(overhang), !is.na(overhang))
  stopifnot(-nchar(watson) < overhang, overhang < nchar(crick))

  crick_rc <- rc(crick, table = table)
  
  up <- substr(watson,   start = max(-overhang, 0) + 1, stop = nchar(crick)  - overhang)
  dn <- substr(crick_rc, start = max(+overhang, 0) + 1, stop = nchar(watson) + overhang)

  if (up != dn) {
    stop("Mismatched basepairs.")
  }
}
