is_staggered <- function(watson, crick) {
  (grepl("-", watson, fixed = TRUE) || !grepl("-", crick, fixed = TRUE))
}


escape_sequence_spec <- function(spec) {
  spec <- gsub("\n", "\\\\n", spec, fixed = TRUE)
  spec
}

reverse <- function(seq) {
  seq <- strsplit(seq, split = "", fixed = TRUE)[[1]]
  seq <- rev(seq)
  paste(seq, collapse = "")
} 

parse_sequence_string <- function(spec) {
  stopifnot(length(spec) == 1L, is.character(spec), !is.na(spec), nzchar(spec))

  ## Single- or double-stranded sequence?
  pattern <- "^([[:alnum:]]+)$"
  if (grepl(pattern, spec)) {
    return(list(type = "ss", specification = spec))
  }

  ## "<watson>;<crick>"?
  pattern <- "^([[:alnum:]-]+);([[:alnum:]-]+)$"
  if (grepl(pattern, spec)) {
    watson <- sub(pattern, "\\1", spec)
    crick  <- sub(pattern, "\\2", spec)
    if (nchar(watson) != nchar(crick)) {
      stop(sprintf("Double-strand sequence string specifies two strands of different lengths (%d != %d): %s", nchar(watson), nchar(crick), sQuote(spec)))
    }
    if (is_staggered(watson, crick)) {
      rcrick <- reverse(crick)
      ## Watson and reverse Crick must not be staggered on the same side
      if (grepl("^-", watson) && grepl("^-", rcrick)) {
        stop(sprintf("Please trim the staggering. Watson and Crick are both staggered at the beginning of the double-stranded sequence: %s", sQuote(spec)))
      } else if (grepl("-$", watson) && grepl("-$", rcrick)) {
        stop(sprintf("Please trim the staggering. Watson and Crick are both staggered at the end of the double-stranded sequence: %s", sQuote(spec)))
      }
    }
    return(list(type = "ds", specification = spec, watson = watson, crick = crick))
  }
  
  ## "<watson>\n<rev crick>"?
  pattern <- "^([[:alnum:]-]+)\n([[:alnum:]-]+)$"
  if (grepl(pattern, spec)) {
    watson <- sub(pattern, "\\1", spec)
    rcrick  <- sub(pattern, "\\2", spec)
    if (nchar(watson) != nchar(rcrick)) {
      stop(sprintf("Double-strand sequence string specifies two strands of different lengths (%d != %d): %s", nchar(watson), nchar(rcrick), sQuote(escape_sequence_spec(spec))))
    }
    crick <- reverse(rcrick)
    if (is_staggered(watson, crick)) {
      ## Watson and reverse Crick must not be staggered on the same side
      if (grepl("^-", watson) && grepl("^-", rcrick)) {
        stop(sprintf("Please trim the staggering. Watson and Crick are both staggered at the beginning of the double-stranded sequence: %s", sQuote(escape_sequence_spec(spec))))
      } else if (grepl("-$", watson) && grepl("-$", rcrick)) {
        stop(sprintf("Please trim the staggering. Watson and Crick are both staggered at the end of the double-stranded sequence: %s", sQuote(escape_sequence_spec(spec))))
      }
    }
    return(list(type = "ds", specification = spec, watson = watson, crick = crick))
  }

  stop(sprintf("Syntax error in sequence string: %s", sQuote(spec)))
}
