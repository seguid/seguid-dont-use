#' The SEGUID of nucleotide and amino-acid sequences
#'
#' @param seq A character string.
#'
#' @return
#' A character string.
#'
#' @examples
#' seguid("ACGTACGTACGT")
#' seguid("acgtACGTacgt")
#' seguid("ACGTACGTACGT") == seguid("acgtACGTacgt")
#'
#' @references
#' 1. Babnigg G, Giometti CS. A database of unique protein sequence
#'    identifiers for proteome studies. Proteomics.
#'    2006 Aug;6(16):4514-22. \doi{10.1002/pmic.200600032}.
#'
#' @importFrom base64enc base64encode
#' @importFrom digest digest
#' @export
seguid <- function(seq) {
  seq <- toupper(seq)
  checksum <- digest::digest(seq, algo = "sha1", serialize = FALSE, raw = TRUE)
  checksum <- base64encode(checksum)
  checksum <- sub("[\n]+$", "", checksum)
  checksum <- sub("[=]$", "", checksum)
  checksum
}


#' @examples
#' useguid("aaa")
#'
#' @rdname seguid
#' @export
useguid <- function(seq) {
  checksum <- seguid(seq)
  checksum <- gsub("+", "-", checksum, fixed = TRUE)
  checksum <- gsub("/", "_", checksum, fixed = TRUE)
  checksum
}


#' @rdname seguid
#' @export
lseguid <- function(seq) {
}


#' @examples
#' cseguid("attt")
#' cseguid("ttta")
#' cseguid("attt") == cseguid("ttta")
#' 
#' @rdname seguid
#' @export
cseguid <- function(seq) {
}



# Internal functions
smallest_rotation <- function(s) {
}

rc <- function(sequence) {
}
        