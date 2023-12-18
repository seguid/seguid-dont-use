#' @importFrom base64enc base64encode
b64encode <- function(s) {
  s <- base64encode(s)
  s <- sub("[=]$", "", s)
  s <- sub("[\n]+$", "", s)
  s
}

b64encode_urlsafe <- function(s) {
  s <- b64encode(s)
  s <- gsub("+", "-", s, fixed = TRUE)
  s <- gsub("/", "_", s, fixed = TRUE)
  s
}


#' @import digest digest
.seguid <- function(seq, table = COMPLEMENT_TABLE_DNA, encoding = b64encode, prefix = "seguid:") {
    assert_table(table)
    assert_in_alphabet(seq, alphabet = names(table))
    stopifnot(is.function(encoding))
    stopifnot(length(prefix) == 1, is.character(prefix), !is.na(prefix))

    checksum <- digest(seq, algo = "sha1", serialize = FALSE, raw = TRUE)
    checksum <- encoding(checksum)
    checksum <- sub("[\n]+$", "", checksum)
    checksum <- sub("[=]$", "", checksum)
    paste(prefix, checksum, sep = "")
}


#' The SEGUID of nucleotide and amino-acid sequences
#'
#' @param seq A character string.
#'
#' @param table A named character string.
#'
#' @param encoding A encoding function.
#'
#' @param prefix A character string prepended the checksum.
#'
#' @return
#' A character string.
#'
#' @examples
#' seguid("ACGTACGTACGT")
#'
#' @references
#' 1. Babnigg G, Giometti CS. A database of unique protein sequence
#'    identifiers for proteome studies. Proteomics.
#'    2006 Aug;6(16):4514-22. \doi{10.1002/pmic.200600032}.
#'
#' @importFrom base64enc base64encode
#' @importFrom digest digest
#' @export
seguid <- function(seq, table = COMPLEMENT_TABLE_DNA, encoding = b64encode, prefix = "seguid:") {
  .seguid(seq, table = table, encoding = encoding, prefix = prefix)
}


#' @rdname seguid
#' @export
slseguid <- function(seq, table = COMPLEMENT_TABLE_DNA, encoding = b64encode_urlsafe, prefix = "slseguid:") {
  .seguid(seq, table = table, encoding = encoding, prefix = prefix)
}


#' @param watson,crick (character string) Two complementary DNA strands.
#'
#' @param overhang (integer) Amount of 3' overhang in the 5' side of
#' the molecule. A molecule with 5' overhang has a negative value.
#' 
#' @rdname seguid
#' @export
dlseguid <- function(watson, crick, overhang, table = COMPLEMENT_TABLE_DNA, encoding = b64encode_urlsafe, prefix = "dlseguid:") {
  assert_anneal(watson, crick, overhang = overhang, table = table)

  lw <- nchar(watson)
  lc <- nchar(crick)
  df <- data.frame(A = c(watson, crick), B = c(crick, watson), C = c(overhang, lw - lc + overhang))
  o <- order(df$A, df$B, df$C, decreasing = FALSE)
  df <- df[o[1],]
  w <- df$A
  c <- df$B
  o <- df$C

  msg <- repr_from_tuple(watson = w, crick = c, overhang = o, table = table, space = "-")

  table2 <- c(table, "-" = "-", "\n" = "\n")
  slseguid(msg, table = table2, encoding = encoding, prefix = prefix)
}



#' @param min_rotation A function.
#' 
#' @rdname seguid
#' @export
scseguid <- function(seq, table = COMPLEMENT_TABLE_DNA, min_rotation = min_rotation_R, prefix = "scseguid:") {
  start <- min_rotation_R(seq)
  slseguid(rotate(seq, amount = start), table = table, prefix = prefix)
}

#' @rdname seguid
#' @export
dcseguid <- function(watson, crick, table = COMPLEMENT_TABLE_DNA, min_rotation = min_rotation_R, prefix = "dcseguid:") {
  stopifnot(nchar(watson) == nchar(crick))
  assert_anneal(watson, crick, overhang = 0, table = table)

  x <- min_rotation(watson)
  y <- min_rotation(crick)
  minwatson <- rotate(watson, amount = x)
  mincrick <- rotate(crick, amount = y)

  ln <- nchar(watson)
  if (minwatson < mincrick) {
    w <- minwatson
    c <- rotate(crick, amount = ln - x)
  } else {
    w <- mincrick
    c <- rotate(watson, amount = ln - y)
  }

  dlseguid(w, c, overhang = 0, table = table, prefix = prefix)
}
