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


min_rotation_py <- function(s, table = COMPLEMENT_TABLE_DNA) {
}
