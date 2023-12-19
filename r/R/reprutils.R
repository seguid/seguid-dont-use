## Emulate Python's 'cleandoc()' in the 'inspect' module
cleandoc_split <- function(s) {
  stopifnot(length(s) == 1L, is.character(s), !is.na(s))
  
  ## Split up into lines
  s <- strsplit(s, split = "\n", fixed = TRUE)[[1]]
  
  ## Drop empty lines
  s <- s[!grepl("^[[:space:]]*$", s)]

  ## Trim right
  s <- gsub("[[:space:]]*$", "", s)

  ## Trim left; common amount of whitespace
  padding <- gsub("[^[:space:]].*", "", s)
  n <- min(nchar(padding))
  if (n > 0) {
    for (kk in seq_along(s)) {
      s[kk] <- substr(s[kk], start = n + 1L, stop = nchar(s[kk]))
    }
  }

  s
}

tuple_from_repr <- function(rpr, table = COMPLEMENT_TABLE_DNA, space = "-") {
  stopifnot(length(space) == 1, is.character(space), !is.na(space))

  ## From 'from string import whitespace; set(whitespace)' in Python
  whitespace <- c('\x0b', '\n', '\t', ' ', '\x0c', '\r')

  assert_in_alphabet(rpr, alphabet = c(names(table), space, whitespace))

  rpr <- cleandoc_split(rpr)
  watson <- rpr[[1]]
  crickrv <- rpr[[2]]
  stopifnot(nchar(watson) == nchar(crickrv))
  prefix <- sprintf("^[%s]+", space)
  suffix <- sprintf("[%s]+$", space)
  stopifnot(!(grepl(prefix, watson) && grepl(prefix, crickrv)))
  stopifnot(!(grepl(suffix, watson) && grepl(suffix, crickrv)))
  assert_in_alphabet(watson, alphabet = c(names(table), space))
  assert_in_alphabet(crickrv, alphabet = c(names(table), space))

  lstrip <- function(s, symbol) {
    pattern <- sprintf("^[%s]+", symbol)
    sub(pattern, "", s)
  }

  strip <- function(s, symbol) {
    pattern <- sprintf("(^[%s]+|[%s]+$)", symbol, symbol)
    gsub(pattern, "", s)
  }

  overhang <- (nchar(watson)  - nchar(lstrip(watson , space))) -
              (nchar(crickrv) - nchar(lstrip(crickrv, space)))

  watson <- strip(watson, space)
  crick <- reverse(strip(crickrv, space))
  
  assert_in_alphabet(watson, alphabet = c(names(table), space))
  assert_in_alphabet(crick , alphabet = c(names(table), space))
  
  table2 <- c(table, space)
  names(table2)[length(table2)] <- space
  assert_table(table2)
  
  assert_anneal(watson, crick, table = table2, overhang = overhang)
  
  list(watson = watson, crick = crick, overhang = overhang)
}


pad <- function(n, symbol = " ") {
  stopifnot(length(n) == 1, is.numeric(n), !is.na(n))
  stopifnot(length(symbol) == 1, is.character(symbol), !is.na(symbol))
  if (n > 0) paste(rep(symbol, times = n), collapse = "") else ""
}


repr_from_tuple <- function(watson, crick, overhang, table = COMPLEMENT_TABLE_DNA, space = "-") {
  assert_anneal(watson, crick, overhang = overhang, table = table)

  roverhang <- overhang + nchar(watson) - nchar(crick)
  rpr <- paste0(
    paste0(pad(+overhang, symbol = space), watson,         pad(-roverhang, symbol = space)),
    "\n",
    paste0(pad(-overhang, symbol = space), reverse(crick), pad(+roverhang, symbol = space))
  )
  assert_in_alphabet(rpr, alphabet = c(names(table), space, "\n"))
  rpr
}

  
