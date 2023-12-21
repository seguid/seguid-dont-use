rotate <- function(seq, amount = 0) {
  stopifnot(length(seq) == 1, is.character(seq), !is.na(seq))
  stopifnot(length(amount) == 1, is.numeric(amount), !is.na(amount))

  ## Nothing to rotate?
  if (nchar(seq) == 0) {
    return(seq)
  }

  amount <- amount %% nchar(seq)

  ## Rotate?
  if (amount > 0) {
    seq <- paste(
      substr(seq, start = amount + 1, stop = nchar(seq)),
      substr(seq, start = 1         , stop = amount    ),
      sep = ""
    )
  }

  seq
}

complementary <- function(seq, table = COMPLEMENT_TABLE_DNA) {
  ## Validate 'table':
  assert_table(table)

  ## Validate 'seq':
  assert_in_alphabet(seq, alphabet = names(table))

  seq <- strsplit(seq, split = "", fixed = TRUE)[[1]]
  seq_c <- table[seq]
  paste(seq_c, collapse = "")
}

reverse <- function(seq) {
  seq <- strsplit(seq, split = "", fixed = TRUE)[[1]]
  seq <- rev(seq)
  paste(seq, collapse = "")
} 

# Reverse complement
rc <- function(seq, table = COMPLEMENT_TABLE_DNA) {
  reverse(complementary(seq, table = table))
}


min_rotation <- function(s) {
  ## Make sure to collate in the 'C' locale
  old_locale <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old_locale))
  Sys.setlocale("LC_COLLATE", "C")
  
  ## Turn string into character vector
  s <- strsplit(s, split = "", fixed = TRUE)[[1]]

  prev <- NULL
  rep <- 0

  ds <- c(s, s)
  lens <- length(s)
  lends <- lens * 2
  old <- 1
  k <- 1
  w <- character(0)
  while (k <= lends) {
    i <- k
    j <- k + 1
#    str(list(k = k, i = i, j = j, ds_i = ds[i], ds_j = ds[j], old = old, prev = prev, w = w))
    
    while (j <= lends && ds[i] <= ds[j]) {
      if (ds[i] == ds[j]) {
        i <- i + 1
      } else {
        i <- k
      }
      j <- j + 1
    }
#    str(list(i = i, j = j, ds_i = ds[i], ds_j = ds[j]))
    
    while (k < i + 1) {
      k <- k + (j - i)
      prev <- w
      w <- ds[old:(k-1)]
#      str(list(k = k, i = i, j = j, old = old, prev = prev, w = w))
      
      old <- k
      if (identical(w, prev)) {
        rep <- rep + 1
      } else {
        prev <- w
        rep <- 1
      }
      if (length(w) * rep == lens) {
#        message("Found: ", old - i)
        return(old - i)
      }
    }
  } ## while (k < lends)
  
  0L
}


rotate_to_min <- function(s) {
  ## Assert that upper-case letters are ordered before lower-case letters
  stopifnot(min_rotation("Aa") == 0)
  
  amount <- min_rotation(s)
  rotate(s, amount = amount)
}
