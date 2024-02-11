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


make_table <- function(definition) {
  stopifnot(length(definition) == 1, is.character(definition), !is.na(definition))
  table <- strsplit(definition, split = ",", fixed = TRUE)[[1]]
  n <- nchar(table)[1]
  stopifnot(nchar(table) == n, n %in% 1:2)
  
  table <- strsplit(table, split = "", fixed = TRUE)
  keys <- vapply(table, FUN = function(x) x[1], FUN.VALUE = NA_character_)
  assert_valid_alphabet(keys)
  
  if (n == 1L) {
    values <- rep("", times = length(keys))
    names(values) <- keys
  } else if (n == 2L) {
    values <- vapply(table, FUN = function(x) x[2], FUN.VALUE = NA_character_)
    assert_valid_alphabet(values)
    names(values) <- keys
    assert_table(values)
  }
  values
}

get_table <- function(spec) {
  stopifnot(length(spec) == 1, is.character(spec), !is.na(spec))
  
  ## Extras? Example: "{DNA}+[-]+[\n]" and "{DNA}+[-\n]"
  extras <- NULL
  pattern <- "(.*)[+][[](.*)[]]$"
  while (grepl(pattern, spec)) {
    spec0 <- spec
    extra <- sub(pattern, "\\2", spec)
    spec <- sub(pattern, "\\1", spec)
    extra <- strsplit(extra, split = "", fixed = TRUE)[[1]]
    extras <- c(extras, extra)
  }

  parts <- strsplit(spec, split = ",", fixed = TRUE)[[1]]
  for (kk in seq_along(parts)) {
    part <- parts[kk]
    if (grepl("^[{][[:alpha:]][[:alnum:]]+[}]$", part)) {
      if (part == "{DNA}") {
        table <- COMPLEMENT_TABLE_DNA
      } else if (part == "{RNA}") {
        table <- COMPLEMENT_TABLE_RNA
      } else if (part == "{IUPAC}") {
        table <- COMPLEMENT_TABLE_IUPAC
      } else if (part == "{protein}") {
        table <- TABLE_IUPAC_PROTEIN
      } else {
        stop("Unknown table: ", sQuote(part))
      }
      part <- paste(sprintf("%s%s", names(table), table), collapse = ",")
      parts[kk] <- part
    }
  }
  parts <- paste(parts, collapse = ",")
  table <- make_table(parts)

  ## Add extras?
  if (length(extras) > 0) {
    names(extras) <- extras
    table <- c(table, extras)
  }

  table
}
