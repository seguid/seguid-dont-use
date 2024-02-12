## Emulate Python's 'cleandoc()' in the 'inspect' module
cleandoc <- function(s) {
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

  paste(s, collapse = "\n")
}

tuple_from_repr <- function(rpr, alphabet = COMPLEMENT_ALPHABET_DNA, space = "-") {
  stopifnot(length(space) == 1, is.character(space), !is.na(space))
  rpr <- cleandoc(rpr)
  assert_in_alphabet(rpr, alphabet = c(names(alphabet), space, "\n"))
  rpr <- strsplit(rpr, split = "\n", fixed = TRUE)[[1]]
  watson <- rpr[[1]]
  crickrv <- rpr[[2]]
  stopifnot(nchar(watson) == nchar(crickrv))
  prefix <- sprintf("^[%s]+", space)
  suffix <- sprintf("[%s]+$", space)
  stopifnot(!(grepl(prefix, watson) && grepl(prefix, crickrv)))
  stopifnot(!(grepl(suffix, watson) && grepl(suffix, crickrv)))
  assert_in_alphabet(watson, alphabet = c(names(alphabet), space))
  assert_in_alphabet(crickrv, alphabet = c(names(alphabet), space))

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
  
  assert_in_alphabet(watson, alphabet = c(names(alphabet), space))
  assert_in_alphabet(crick , alphabet = c(names(alphabet), space))
  
  alphabet2 <- c(alphabet, space)
  names(alphabet2)[length(alphabet2)] <- space
  assert_alphabet(alphabet2)
  
  assert_anneal(watson, crick, alphabet = alphabet2, overhang = overhang)
  
  list(watson = watson, crick = crick, overhang = overhang)
}


pad <- function(n, symbol = " ") {
  stopifnot(length(n) == 1, is.numeric(n), !is.na(n))
  stopifnot(length(symbol) == 1, is.character(symbol), !is.na(symbol))
  if (n > 0) paste(rep(symbol, times = n), collapse = "") else ""
}


repr_from_tuple <- function(watson, crick, overhang, alphabet = COMPLEMENT_ALPHABET_DNA, space = "-") {
  assert_anneal(watson, crick, overhang = overhang, alphabet = alphabet)

  roverhang <- overhang + nchar(watson) - nchar(crick)
  rpr <- paste0(
    paste0(pad(+overhang, symbol = space), watson,         pad(-roverhang, symbol = space)),
    "\n",
    paste0(pad(-overhang, symbol = space), reverse(crick), pad(+roverhang, symbol = space))
  )
  assert_in_alphabet(rpr, alphabet = c(names(alphabet), space, "\n"))
  rpr
}

  

watson_crick_from_tuple <- function(watson, crick, overhang, space = "-") {
  ## Nothing to do?
  if (overhang == 0L) {
    return(list(watson, crick))
  }
  
  if (overhang > 0) {
    watson_pad_left  <- pad(+overhang, symbol = "-")
    watson_pad_right <- pad(-overhang+nchar(crick)-nchar(watson), symbol = "-")
    crick_pad_left   <- pad(+overhang+nchar(watson)-nchar(crick), symbol = "-")
    crick_pad_right  <- pad(-overhang, symbol = "-")
  } else if (overhang < 0) {
    crick_pad_left   <- pad(-overhang, symbol = "-")
    crick_pad_right  <- pad(+overhang+nchar(watson)-nchar(crick), symbol = "-")
    watson_pad_left  <- pad(-overhang+nchar(crick)-nchar(watson), symbol = "-")
    watson_pad_right <- pad(+overhang, symbol = "-")
  }
  
  watson <- paste0(watson_pad_left, watson, watson_pad_right)
  crick  <- paste0(crick_pad_left, crick, crick_pad_right)

  return(list(watson, crick))
}


dsseq_to_tuple <- function(watson, crick, overhang) {
  ## Staggeredness specified via dash symbols ('-')
  if (!grepl("-", watson, fixed = TRUE) && !grepl("-", crick, fixed = TRUE)) {
    return(list(watson = watson, crick = crick, overhang = overhang))
  }

  pattern <- "^([-]*)([^-]*)([-]*)$"
  stopifnot(grepl(pattern, watson), grepl(pattern, crick), nchar(watson) == nchar(crick))
  watson_trim <- sub(pattern, "\\2", watson)
  watson_pad_left <- sub(pattern, "\\1", watson)
  watson_pad_right <- sub(pattern, "\\3", watson)
  crick_trim <- sub(pattern, "\\2", crick)
  crick_pad_left <- sub(pattern, "\\1", crick)
  crick_pad_right <- sub(pattern, "\\3", crick)
  
  if (nchar(watson_pad_left) > 0) {
    stopifnot(nchar(crick_pad_right) == 0)
    overhang <- nchar(watson_pad_left)
  } else if (nchar(crick_pad_left) > 0) {
    stopifnot(nchar(watson_pad_right) == 0)
    overhang <- -nchar(crick_pad_left)
  }

  watson <- watson_trim
  crick <- crick_trim
  
  list(watson = watson, crick = crick, overhang = overhang)
} ## dsseq_to_tuple()
