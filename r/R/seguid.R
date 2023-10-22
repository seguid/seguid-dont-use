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

  if (!is_dna_sequence(seq) && !is_rna_sequence(seq) && !is_amino_acid_sequence(seq)) {
     stop("Sequence is neither a DNA sequence, an RNA sequence, nor an amino-acid sequence: ", sQuote(seq))
  }

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
  seq <- toupper(seq)
  if (!is_dna_sequence(seq) && !is_rna_sequence(seq)) {
     stop("Sequence is neither a DNA sequence nor an RNA sequence: ", sQuote(seq))
  }

  checksum <- seguid(seq)
  checksum <- gsub("+", "-", checksum, fixed = TRUE)
  checksum <- gsub("/", "_", checksum, fixed = TRUE)
  checksum
}

#' @examples
#' lseguid_blunt("ACGGGT")
#'
#' @rdname seguid
#' @export
lseguid_blunt <- function(seq) {
  seq <- toupper(seq)
  if (!is_dna_sequence(seq) && !is_rna_sequence(seq)) {
     stop("Sequence is neither a DNA sequence nor an RNA sequence: ", sQuote(seq))
  }
  useguid(min(seq, rc(seq)))
}


#' @param watson,crick (character string) Two complementary DNA strands.
#'
#' @param overhang (integer) Amount of 3' overhang in the 5' side of
#' the molecule. A molecule with 5' overhang has a negative value.
#' 
#' @rdname seguid
#' @export
lseguid_sticky <- function(watson, crick, overhang) {
  watson <- toupper(watson)
  crick <- toupper(crick)

  if (!is_dna_sequence(watson) && !is_rna_sequence(watson)) {
     stop("The 'watson' sequence is neither a DNA sequence nor an RNA sequence: ", sQuote(watson))
  }

  if (!is_dna_sequence(crick) && !is_rna_sequence(crick)) {
     stop("The 'crick' sequence is neither a DNA sequence nor an RNA sequence: ", sQuote(crick))
  }

  lw <- nchar(watson)
  lc <- nchar(crick)
  
  if (overhang == 0L && lw == lc) {
    lseguid_blunt(watson)
  } else {
    w <- min(watson, crick)
    c <- min(crick, watson)
    o <- min(overhang, lw - lc + overhang)
    spaces_w <- if (o > 0L) rep(" ", times =  o) else ""
    spaces_c <- if (o < 0L) rep(" ", times = -o) else ""
    seq_w <- paste(spaces_w,     w , sep = "")
    seq_c <- paste(spaces_c, rev(c), sep = "")
    seq <- paste(seq_w, seq_c, sep = "\n")
    useguid(seq)
  }
}


#' @examples
#' cseguid("attt")
#' cseguid("ttta")
#' cseguid("attt") == cseguid("ttta")
#' 
#' @rdname seguid
#' @export
cseguid <- function(seq) {
  seq <- toupper(seq)
  if (!is_dna_sequence(seq) && !is_rna_sequence(seq)) {
     stop("Sequence is neither a DNA sequence nor an RNA sequence: ", sQuote(seq))
  }
}



# ---------------------------------------------------------
# Internal functions
# ---------------------------------------------------------
dna_complement <- c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
rna_complement <- c("A" = "U", "C" = "G", "G" = "C", "U" = "A")
amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

ambiguous_dna_complement <- c(
    dna_complement,
    rna_complement,
    "M" = "K", "R" = "Y",
    "W" = "W", "S" = "S", "Y" = "R",
    "K" = "M", "V" = "B", "H" = "D",
    "D" = "H", "B" = "V", "X" = "X",
    "N" = "N"
)


is_dna_sequence <- function(seq) {
  seq <- strsplit(seq, split = "", fixed = TRUE)[[1]]
  all(seq %in% dna_complement)
}

is_rna_sequence <- function(seq) {
  seq <- strsplit(seq, split = "", fixed = TRUE)[[1]]
  all(seq %in% rna_complement)
}

is_amino_acid_sequence <- function(seq) {
  seq <- strsplit(seq, split = "", fixed = TRUE)[[1]]
  all(seq %in% amino_acids)
}


space <- function(n) {
  paste(rep(" ", times = n), collapse = "")
}

reverse_sequence <- function(sequence) {
  sequence <- strsplit(sequence, split = "", fixed = TRUE)[[1]]
  paste(sequence, collapse = "")
} 

# Reverse complement
reverse_complement <- function(sequence) {
  if (!is_dna_sequence(sequence) && !is_rna_sequence(sequence)) {
     stop("Sequence is neither a DNA sequence nor an RNA sequence: ", sQuote(sequence))
  }
  sequence <- toupper(sequence)
  sequence <- strsplit(sequence, split = "", fixed = TRUE)[[1]]
  sequence_c <- ambiguous_dna_complement[sequence]
  stopifnot(!anyNA(sequence_c))
  sequence_cr <- rev(sequence_c)
  paste(sequence_cr, collapse = "")
}

rc <- reverse_complement


tuple_from_representation <- function(bfr) {
  ## Split up into lines
  bfr <- strsplit(bfr, split = "\n", fixed = TRUE)[[1]]
  
  ## Drop empty lines
  bfr <- bfr[!grepl("^[[:space:]]*$", bfr)]

  ## Drop lines with |||
  bfr <- bfr[!grepl("^[[:space:]]*[|]+[[:space:]]*$", bfr)]

  stopifnot(length(bfr) == 2L)

  ## Drop 5' and 3' prefixes and suffixes
  bfr <- gsub("([35]'-|-[35]')", "", bfr)
  
  ## Trim right
  bfr <- gsub("[[:space:]]*$", "", bfr)
  
  pattern <- "^([[:space:]]*)([^[:space:]]+)$"
  ls <- gsub(pattern, "\\1", bfr)
  nls <- nchar(ls) - min(nchar(ls))
  ms <- gsub(pattern, "\\2", bfr)

  stopifnot(
    is_dna_sequence(toupper(ms[1])) || is_rna_sequence(toupper(ms[1])),
    is_dna_sequence(toupper(ms[2])) || is_rna_sequence(toupper(ms[2]))
  )
  
  list(watson = ms[1], cricket = ms[2], overhang = nls[1])
}  


# Smallest rotation of a string.
#
# Algorithm described in:
# 1. Pierre Duval , Jean. 1983. Factorizing Words over an Ordered
#    Alphabet. Journal of Algorithms & Computational Technology
#    4 (4) (December 1): 363–381
# 2. Algorithms on strings and sequences based on Lyndon words,
#    David Eppstein 2011 <https://gist.github.com/dvberkel/1950267>
#
# Examples:
# smallest_rotation("taaa") == "aaat"
#
smallest_rotation <- function(s) {
#    from array import array as _array
#    prev <- ""
#    rep <- 0
#    ds <- c("u", paste0(s, s))
#    lens <- nchar(s)
#    lends <- nchar(ds)
#    old <- 0
#    k <- 0
#    w <- ""
#    
#    while (k < lends) {
#        i <- k
#        j <- k + 1
#        
#        while (j < lends && ds[i] <= ds[j]) {
#            i <- (ds[i] == ds[j]) and i + 1 or k
#            j <- j + 1
#        }
#        
#        while (k < i + 1) {
#            k <- k + (j - i)
#            prev <- w
#            w <- ds[old:k]
#            old <- k
#            
#            if (w == prev) {
#                rep <- rep + 1
#            } else {
#                prev <- w
#                rep <- 1
#            }
#
#            ## Done?
#            if (nchar(w) * rep == lens) {
#                return(paste(rep(w, times = rep), collapse = ""))
#            }
#        }
#    }
}

