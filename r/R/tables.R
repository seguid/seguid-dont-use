# Definition of Complementary DNA Symbols
COMPLEMENT_TABLE_DNA <- c(G="C",
                          A="T",
                          C="G",
                          T="A")


# Definition of Complementary RNA Symbols
COMPLEMENT_TABLE_RNA <- c(G="C",
                          A="U",
                          C="G",
                          U="A")


# Definition of Complementary IUPAC Ambigous DNA Symbols
COMPLEMENT_TABLE_IUPAC <- c(COMPLEMENT_TABLE_DNA, c(B="V",
                                                    D="H",
                                                    H="D",
                                                    K="M",
                                                    M="K",
                                                    S="S",
                                                    V="B",
                                                    W="W",
                                                    N="N"))


TABLE_IUPAC_PROTEIN <- c(A="",
                         C="",
                         D="",
                         E="",
                         F="",
                         G="",
                         H="",
                         I="",
                         K="",
                         L="",
                         M="",
                         N="",
                         P="",
                         Q="",
                         R="",
                         S="",
                         T="",
                         V="",
                         W="",
                         Y="")


get_table <- function(name) {
  stopifnot(length(name) == 1, is.character(name), !is.na(name))
  
  ## Extras? Example: "dna+[-]+[\n]" and "dna+[-\n]"
  extras <- NULL
  pattern <- "(.*)[+][[](.*)[]]$"
  while (grepl(pattern, name)) {
    name0 <- name
    extra <- sub(pattern, "\\2", name)
    name <- sub(pattern, "\\1", name)
    extra <- strsplit(extra, split = "", fixed = TRUE)[[1]]
    extras <- c(extras, extra)
  }
  
  if (name == "dna") {
    table <- COMPLEMENT_TABLE_DNA
  } else if (name == "rna") {
    table <- COMPLEMENT_TABLE_RNA
  } else if (name == "iupac") {
    table <- COMPLEMENT_TABLE_IUPAC
  } else if (name == "protein") {
    table <- TABLE_IUPAC_PROTEIN
  } else {
    stop("Unknown table: ", sQuote(name))
  }

  ## Add extras?
  if (length(extras) > 0) {
    names(extras) <- extras
    table <- c(table, extras)
  }

  table
}
