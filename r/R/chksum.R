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


with_prefix <- function(s, prefix, form = c("long", "short", "both")) {
  form <- match.arg(form)
  
  checksum <- sub("^(|(l|c)(s|d))*seguid=", "", s)
  assert_checksum(checksum, prefix = "")

  if (form == "both") form <- c("short", "long")

  res <- character(0L)
  for (ff in form) {
    if (ff == "long") {
      res <- c(res, paste0(prefix, checksum))
    } else if (ff == "short") {
      res <- c(res, substr(checksum, start = 1L, stop = 6L))
    }
  }
  
  res
}

.seguid <- function(seq, alphabet, encoding, prefix = "") {
    assert_alphabet(alphabet)
    assert_in_alphabet(seq, alphabet = names(alphabet))
    stopifnot(is.function(encoding))
    stopifnot(length(prefix) == 1, is.character(prefix), !is.na(prefix))

    checksum <- encoding(seq)

    checksum <- paste0(prefix, checksum)
    assert_checksum(checksum, prefix = prefix)
    checksum
}


#' SEGUID checksums for linear, circular, single- and double-stranded sequences
#'
#' @param seq (character string) The sequence for which the checksum
#' should be calculated.  The sequence may only comprise of symbols
#' in the alphabet specified by the `alphabet` argument.
#'
#' @param alphabet (character string) The type of sequence used.
#' If `"{DNA}"` (default), then the input is a DNA sequence.
#' If `"{IUPAC}"`, then the input is a DNA sequence specified with
#' IUPAC ambigous DNA symbols (3).
#' If `"{RNA}"`, then the input is an RNA sequence.
#' If `"{protein}"`, then the input is an amino-acid sequence.
#' A custom alphabet may also be used.
#' A non-complementary alphabet is specified as a comma-separated
#' set of single symbols, e.g. `"X,Y,Z"`.
#' A complementary alphabet is specified as a comma-separated
#' set of paired symbols, e.g. `"AT,TA,CG,GC"`.
#' It is also possible to extend a pre-defined alphabet, e.g.
#' `"{DNA},XY,YX"`.
#'
#' @param form (character string) How the checksum is presented.
#' If `"long"` (default), the full-length checksum is outputted.
#' If `"short"`, the short, six-digit checksum is outputted.
#' If `"both"`, both the short and the long checksums are outputted.
#'
#' @return
#' The SEGUID functions return a single character string, if `form` is
#' either `"long"` or `"short"`. If `form` is `"both"`, then a character
#' vector of length two is return, where the first component holds the
#' "short" checksum and the second the "long" checksum.
#' The long checksum, without the prefix, is string with 27 characters.
#' The short checksum, without the prefix, is the first six characters
#' of the long checksum.
#' All checksums are prefixed with a label indicating which SEGUID
#' method was used.
#' Except for `seguid()`, which uses _base64_ encoding, all functions
#' produce checksums using the _base64url_ encoding ("Base 64 Encoding
#' with URL and Filename Safe Alphabet").
#'
#' `seguid()` calculates the SEGUID v1 checksum for a linear,
#' single-stranded sequence. 
#'
#' @section Base64 and Base64url encodings:
#' The base64url encoding is the base64 encoding with non-URL-safe characters
#' substituted with URL-safe ones. Specifically, the plus symbol (`+`) is
#' replaced by the minus symbol (`-`), and the forward slash (`/`) is
#' replaced by the underscore symbol (`_`).
#'
#' The Base64 checksum, which is used for the original SEGUID checksum,
#' is not guaranteed to comprise symbols that can
#' safely be used as-is in Uniform Resource Locator (URL). Specifically,
#' it may consist of forward slashes (`/`) and plus symbols (`+`), which
#' are characters that carry special meaning in a URL.
#' For the same reason, a Base64 checksum cannot safely be used
#' as a file or directory name, because it may have a forward slash.
#'
#' The checksum returned is always 27-character long. This is because the
#" SHA-1 hash (4) is 160-bit long (20 bytes), which result in the encoded
#' represention always end with a padding character (`=`) so that the length
#' is a multiple of four character. We relax this requirement, by dropping
#' the padding character.
#'
#' @example incl/seguid.R
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
seguid <- function(seq, alphabet = "{DNA}", form = c("long", "short", "both")) {
  if (nchar(seq) == 0) {
    stop("A sequence must not be empty")
  }
  
  alphabet2 <- get_alphabet(alphabet)
  with_prefix(.seguid(seq, alphabet = alphabet2, encoding = sha1_b64encode), prefix = "seguid=", form = form)
}


#' @return
#' `lsseguid()` calculates the SEGUID v2 checksum for a linear,
#' single-stranded sequence.
#'
#' @rdname seguid
#' @export
lsseguid <- function(seq, alphabet = "{DNA}", form = c("long", "short", "both")) {
  if (nchar(seq) == 0) {
    stop("A sequence must not be empty")
  }
  
  alphabet2 <- get_alphabet(alphabet)
  with_prefix(.seguid(seq, alphabet = alphabet2, encoding = sha1_b64encode_urlsafe), prefix = "lsseguid=", form = form)
}


#' @return
#' `csseguid()` calculates the SEGUID v2 checksum for a circular,
#' single-stranded sequence.
#'
#' @rdname seguid
#' @export
csseguid <- function(seq, alphabet = "{DNA}", form = c("long", "short", "both")) {
  if (nchar(seq) == 0) {
    stop("A sequence must not be empty")
  }
  
  with_prefix(lsseguid(rotate_to_min(seq), alphabet = alphabet), prefix = "csseguid=", form = form)
}


#' @param watson,crick (character strings) Two reverse-complementary DNA
#' sequences. Both sequences should be specified in the 5'-to-3' direction.
#'
#' @return
#' `ldseguid()` calculates the SEGUID v2 checksum for a linear,
#' double-stranded sequence.
#'
#' @rdname seguid
#' @export
ldseguid <- function(watson, crick, alphabet = "{DNA}", form = c("long", "short", "both")) {
  ## Make sure to collate in the 'C' locale
  old_locale <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old_locale))
  Sys.setlocale("LC_COLLATE", "C")

  if (nchar(watson) == 0 || nchar(crick) == 0) {
    stop("A sequence must not be empty")
  }

  alphabet2 <- paste0(alphabet, "+[-\n]")
  assert_complementary(watson, crick, alphabet = alphabet2)

  rcrick <- reverse(crick)
  if (is_seq_less_than(watson, rcrick)) {
    spec <- paste(watson, rcrick, sep = "\n")
  } else {
    spec <- paste(rcrick, watson, sep = "\n")
  }
  with_prefix(lsseguid(spec, alphabet = alphabet2), prefix = "ldseguid=", form = form)
}


#' @return
#' `cdseguid()` calculates the SEGUID v2 checksum for a circular,
#' double-stranded sequence.
#'
#' @rdname seguid
#' @export
cdseguid <- function(watson, crick, alphabet = "{DNA}", form = c("long", "short", "both")) {
  if (nchar(watson) == 0 || nchar(crick) == 0) {
    stop("A sequence must not be empty")
  }
  
  stopifnot(nchar(watson) == nchar(crick))
  
  assert_complementary(watson, crick, alphabet = alphabet)

  watson_min <- rotate_to_min(watson)
  crick_min <- rotate_to_min(crick)

  ## Keep the "minimum" of the two variants
  alphabet2 <- get_alphabet(alphabet)
  if (is_seq_less_than(watson_min, crick_min)) {
      w <- watson_min
      c <- rc(w, alphabet = alphabet2)  ## FIXME [#74]
  } else {
      w <- crick_min
      c <- rc(w, alphabet = alphabet2)  ## FIXME [#74]
  }

  with_prefix(ldseguid(watson = w, crick = c, alphabet = alphabet), prefix = "cdseguid=", form = form)
}
