#' @importFrom base64enc base64encode
b64encode <- function(s) {
  base64encode(s) 
}

b64encode_urlsafe <- function(s) {
  s <- b64encode(s)
  s <- gsub("+", "-", s, fixed = TRUE)
  s <- gsub("/", "_", s, fixed = TRUE)
  s
}

#' @import digest digest
sha1_b64encode <- function(seq) {
  checksum <- digest(seq, algo = "sha1", serialize = FALSE, raw = TRUE)
  checksum <- b64encode(checksum)

  ## Drop newlines (just in case)
  checksum <- sub("[\n]+$", "", checksum)

  ## SHA-1 (160 bits = 20 bytes = 40 hexadecimal character) needs
  ## at most 160/log2(64) = 26.6667 = 27 symbols. Base64 pads to
  ## multiples of 4 symbols, i.e. 28 symbols.  Thus, the last
  ## symbol is always a pad symbol, when using SHA-1. This is
  ## why we drop the last symbol.
  checksum <- sub("[=]$", "", checksum)

  checksum
}    

#' @import digest digest
sha1_b64encode_urlsafe <- function(seq) {
  checksum <- digest(seq, algo = "sha1", serialize = FALSE, raw = TRUE)
  checksum <- b64encode_urlsafe(checksum)

  ## Drop newlines (just in case)
  checksum <- sub("[\n]+$", "", checksum)

  ## SHA-1 (160 bits = 20 bytes = 40 hexadecimal character) needs
  ## at most 160/log2(64) = 26.6667 = 27 symbols. Base64 pads to
  ## multiples of 4 symbols, i.e. 28 symbols.  Thus, the last
  ## symbol is always a pad symbol, when using SHA-1. This is
  ## why we drop the last symbol.
  checksum <- sub("[=]$", "", checksum)

  checksum
}    


with_prefix <- function(s, prefix) {
  paste0(prefix, sub("^(|(l|c)(s|d))*seguid-", "", s))
}

.seguid <- function(seq, table, encoding, prefix = "") {
    assert_table(table)
    assert_in_alphabet(seq, alphabet = names(table))
    stopifnot(is.function(encoding))
    stopifnot(length(prefix) == 1, is.character(prefix), !is.na(prefix))

    checksum <- encoding(seq)

    checksum <- paste0(prefix, checksum)
    assert_checksum(checksum, prefix = prefix)
    checksum
}


#' SEGUID checksum for protein and linear single-stranded DNA
#'
#' @param seq (character string) The sequence for which the checksum
#' should be calculated.  The sequence may only comprise of characters
#' in the alphabet specified by the `table` argument.
#'
#' @param table (character string) The type of sequence used.
#' If `"{DNA}"` (default), then the input is a DNA sequence.
#' If `"iupac"`, then the input is a DNA sequence specified with
#' IUPAC ambigous DNA symbols (3).
#' If `"{RNA}"`, then the input is an RNA sequence.
#' If `"{protein}"`, then the input is an amino-acid sequence.
#'
#' @return
#' `seguid()` returns a character string composed of the prefix `seguid-`
#' followed by a _base64_ encoding (2) ("Base 64 Encoding").
#' A base64 encoding is always 27-character long, and may comprise
#' non-URL-safe characters.
#'
#' @details
#' The Sequence Globally Unique Identifiers (SEGUID) (1) is defined as the
#' Base64-encoded SHA-1 checksum (4) calculated for the sequence in
#' uppercase with the trailing pad character `=` removed.
#' In contrast to the original implementation (1), this function returns
#' the SEGUID checksum prefixed with `seguid-`.
#'
#' @section Base64 checksums are neither filename nor URL safe:
#' The Base64 checksum is not guaranteed to comprise symbols that can
#' safely be used as-is in Uniform Resource Locator (URL). Specifically,
#' it may consist of forward slashes (`/`) and plus symbols (`+`), which
#' are characters that carry special meaning in a URL.
#' For the same reason, a Base64 checksum can be guaranteed to be used
#' as a file or directory name, because it may have a forward slash.
#'
#' @examples
#' ## Linear single-stranded DNA:
#' ## GATTACA
#' seguid("GATTACA")
#' #> seguid-tp2jzeCM2e3W4yxtrrx09CMKa/8
#'
#' ## Linear single-stranded DNA
#' ## GATTACA
#' lsseguid("GATTACA")
#' #> lsseguid-tp2jzeCM2e3W4yxtrrx09CMKa_8
#'
#' ## Circular single-stranded DNA
#' ## GATTACA = ATTACAG = ... = AGATTAC
#' csseguid("GATTACA")
#' #> csseguid-mtrvbtuwr6_MoBxvtm4BEpv-jKQ
#'
#' ## Linear double-stranded DNA
#' ## GATTACA
#' ## CTAATGT
#' seguid::ldseguid("GATTACA", "TGTAATC", overhang = 0)
#' #> ldseguid-XscjVNyZarYrROVgGXUCleJcMC
#'
#' ## Circular double-stranded DNA
#' ## GATTACA = ATTACAG = ... = AGATTAC
#' ## CTAATGT = TAATGTC = ... = TCTAATG
#' seguid::cdseguid("GATTACA", "TGTAATC")
#' #> cdseguid-zCuq031K3_-40pArbl-Y4N9RLnA
#'
#' @references
#' 1. Babnigg, G., Giometti, CS. A database of unique protein sequence
#'    identifiers for proteome studies. Proteomics.
#'    2006 Aug;6(16):4514-22. \doi{10.1002/pmic.200600032}.
#' 2. Josefsson, S., The Base16, Base32, and Base64 Data Encodings,
#'    RFC 4648, \doi{10.17487/RFC4648}, October 2006,
#'    <https://www.rfc-editor.org/info/rfc4648>.
#' 3. Wikpedia article 'Nucleic acid notation', December 2023.
#'    <https://en.wikipedia.org/wiki/Nucleic_acid_notation>.
#' 4. Wikipedia article 'SHA-1' (Secure Hash Algorithm 1), December 2023.
#'    <https://en.wikipedia.org/wiki/SHA-1>.
#'
#' @importFrom base64enc base64encode
#' @importFrom digest digest
#' @export
seguid <- function(seq, table = "{DNA}") {
  if (nchar(seq) == 0) {
    stop("A sequence must not be empty")
  }
  
  table2 <- get_table(table)
  with_prefix(.seguid(seq, table = table2, encoding = sha1_b64encode), prefix = "seguid-")
}


#' @return
#' `lsseguid()` returns a character string composed of the prefix `lsseguid-`
#' followed by a _base64url_ encoding (2) ("Base 64 Encoding with URL and
#' Filename Safe Alphabet").
#'
#' @section Base64url checksums are filename and URL safe:
#' The base64url encoding is the base64 encoding with non-URL-safe characters
#' substituted with URL-safe ones. Specifically, the plus symbol (`+`) is
#' replaced by the minus symbol (`-`), and the forward slash (`/`) is
#' replaced by the underscore symbol (`_`).
#'
#' @rdname seguid
#' @export
lsseguid <- function(seq, table = "{DNA}") {
  if (nchar(seq) == 0) {
    stop("A sequence must not be empty")
  }
  
  table2 <- get_table(table)
  with_prefix(.seguid(seq, table = table2, encoding = sha1_b64encode_urlsafe), prefix = "lsseguid-")
}


#' @return
#' `csseguid()` returns a character string composed of the prefix `csseguid-`
#' followed by a _base64url_ encoding.
#'
#' @rdname seguid
#' @export
csseguid <- function(seq, table = "{DNA}") {
  if (nchar(seq) == 0) {
    stop("A sequence must not be empty")
  }
  
  with_prefix(lsseguid(rotate_to_min(seq), table = table), prefix = "csseguid-")
}


#' @param watson,crick (character string) Two reverse-complementary DNA strands.
#'
#' @param overhang (integer) Amount of 3' overhang in the 5' side of
#' the molecule. A molecule with 5' overhang has a negative value.
#' 
#' @rdname seguid
#' @export
ldseguid <- function(watson, crick, overhang, table = "{DNA}") {
  if (nchar(watson) == 0 || nchar(crick) == 0) {
    stop("A sequence must not be empty")
  }
  
  table2 <- get_table(table)
  assert_anneal(watson, crick, overhang = overhang, table = table2)

  lw <- nchar(watson)
  lc <- nchar(crick)
  df <- data.frame(A = c(watson, crick), B = c(crick, watson), C = c(overhang, lw - lc + overhang))
  o <- order(df$A, df$B, df$C, decreasing = FALSE)
  df <- df[o[1],]
  w <- df$A
  c <- df$B
  o <- df$C

  msg <- repr_from_tuple(watson = w, crick = c, overhang = o, table = table2, space = "-")

  table2 <- paste0(table, "+[-\n]")
  with_prefix(lsseguid(msg, table = table2), prefix = "ldseguid-")
}


#' @return
#' `cdseguid()` returns a character string composed of the prefix `cdseguid-`
#' followed by a _base64url_ encoding.
#'
#' @rdname seguid
#' @export
cdseguid <- function(watson, crick, table = "{DNA}") {
  if (nchar(watson) == 0 || nchar(crick) == 0) {
    stop("A sequence must not be empty")
  }
  
  stopifnot(nchar(watson) == nchar(crick))
  table2 <- get_table(table)
  assert_anneal(watson, crick, overhang = 0, table = table2)

  watson_min <- rotate_to_min(watson)
  crick_min <- rotate_to_min(crick)

  ## Keep the "minimum" of the two variants
  if (is_seq_less_than(watson_min, crick_min)) {
      w <- watson_min
  } else {
      w <- crick_min
  }

  with_prefix(ldseguid(w, rc(w, table = table2), overhang = 0, table = table), prefix = "cdseguid-")
}
