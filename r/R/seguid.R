#' The SEGUID of nucleotide and amino-acid sequences
#'
#' @param seq A character string.
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
#' @export
seguid <- function(seq) {
}


#' @rdname seguid
#' @export
lseguid <- function(seq) {
}

#' @examples
#' useguid("aaa")
#'
#' @rdname seguid
#' @export
useguid <- function(seq) {
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
        