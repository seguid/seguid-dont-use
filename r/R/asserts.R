escape_symbols <- function(s) {
  for (kk in seq_along(s)) {
    ch <- s[kk]
    if (grepl("^[[:space:]]$", ch)) {
      ch <- gsub('(^"|"$)', "", deparse(ch))
      s[kk] <- ch;
    }
  }
  s
}

assert_alphabet <- function(alphabet) {
  stopifnot(length(alphabet) > 0, is.character(alphabet), !anyNA(alphabet), all(nchar(alphabet) == 1))
  
  ## Allow only for A-Z, a-Z, '-', '\n'
  unknown <- setdiff(alphabet, c(LETTERS, letters, "-", "\n"))
  if (length(unknown) > 0) {
    unknown <- escape_symbols(unique(sort(unknown)))
    stop(sprintf("Non-supported symbols in alphabet: [n=%d] %s", length(unknown), paste(sQuote(unknown), collapse = ", ")))
  }  
}  

assert_in_alphabet <- function(seq, alphabet) {
  stopifnot(length(seq) == 1, is.character(seq), !is.na(seq))
  assert_alphabet(alphabet)

  ## Nothing to do?
  if (nchar(seq) == 0) {
    return(seq)
  }

  seq <- strsplit(seq, split = "", fixed = TRUE)[[1]]
  unknown <- setdiff(seq, alphabet)
  if (length(unknown) > 0) {
    unknown <- escape_symbols(unique(sort(unknown)))
    stop(sprintf("Sequence symbols not in alphabet: [n=%d] %s", length(unknown), paste(sQuote(unknown), collapse = ", ")))
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


assert_checksum <- function(checksum) {
  stopifnot(length(checksum) == 1L, !is.na(checksum))
  pattern <- "^(|sl|sc|ds|dc)seguid:"
  stopifnot(grepl(pattern, checksum))
  checksum <- sub(pattern, "", checksum)
  stopifnot(nchar(checksum) == 27)
}
