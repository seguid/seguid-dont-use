cli_help_string <- '
{{ package }}: {{ title }}

Usage:
 Rscript -e seguid::seguid [options] <<< "<sequence>"

Options:  
 --help             Display the full help page with examples
 --version          Output version of this software

Examples:
Rscript -e seguid::seguid --version
Rscript -e seguid::seguid --help

echo "ACGT" | Rscript -e seguid::lsseguid
Rscript -e seguid::lsseguid <<< "ACGT"
Rscript -e seguid::ldseguid <<< $\'-CGT\nTGCA\'
Rscript -e seguid::cdseguid <<< $\'ACGT\nTGCA\'

Version: {{ version }}
Copyright: Henrik Bengtsson (2023-2024)
License: MIT
'

#' @importFrom utils capture.output file_test str
cli_call_fcn <- function(..., alphabet = "{DNA}", file = NULL, debug = FALSE, fcn) {
  if (is.character(fcn)) {
    fcn <- get(fcn, mode = "function", envir = getNamespace(.packageName), inherits = FALSE)
  }
  stopifnot(length(debug) == 1, is.logical(debug), !is.na(debug))

  alphabet2 <- get_alphabet(alphabet)
  
  seq <- NULL
  
  args <- list(...)
  nargs <- length(args)
  if (nargs > 0) {
    names <- names(args)
    if (is.null(names)) names <- rep("", times = nargs) 
    ## At most one unnamed CLI option
    unnamed <- which(nchar(names) == 0)
    n <- length(unnamed)
    if (n > 2) {
      unknown <- args[seq_len(n - 2)]
      stop("Unknown option(s): ", paste(sQuote(unknown), collapse = ", "))
    }
    if (n >= 1) {
      if (!is.null(file)) {
        stop("Option --file=<filename> already specified")
      }
      seq <- unlist(args[unnamed])
      seq <- paste(seq, collapse = "\n")
    }
  }

  if (is.null(seq)) {
    ## Get sequence input from file?
    if (!is.null(file)) {
      if (!file_test("-f", file)) {
        stop("No such file: ", sQuote(file))
      }
      seq <- readLines(file)
      seq <- paste(seq, collapse = "\n")
    } else {
      ## Get arguments from the standard input
      seq <- readLines("stdin")
      seq <- paste(seq, collapse = "\n")
    }
  }

  if (debug) {
    message(sprintf("Sequence data:\n%s", seq))
    message(sprintf("Arguments:\n%s", paste(capture.output(str(args)), collapse = "\n")))
  }

  ## Parse sequence string. This will throw an error if not meeting the minimal specifications
  argnames <- names(formals(fcn))
  if (is.element("crick", argnames)) {
    nseq <- length(strsplit(seq, split = "\n", fixed = TRUE)[[1]])
    if (nseq == 1) {
      stop(sprintf("Specified sequence string contains a single sequence, but expected two: %s", seq))
    } else {
      seq_spec <- parse_sequence_string(seq)
      args2 <- list(watson = seq_spec[["watson"]], crick = seq_spec[["crick"]])
    }
    if (debug) {
      msg <- sprintf("Double-stranded sequence pair:\nwatson=%s\ncrick=%s", sQuote(args2[[1]]), sQuote(args2[[2]]))
    }
    if (debug) message(msg)
  } else {
    args2 <- list(seq)
  }
  args <- c(args2, args, alphabet = alphabet)
  if (debug) {
    message(sprintf("Arguments:\n%s", paste(capture.output(str(args)), collapse = "\n")))
  }

  res <- do.call(fcn, args = args)
  if (debug) {
    message(sprintf("Result:\n%s", paste(capture.output(str(res)), collapse = "\n")))
  }

  res <- paste(res, collapse = " ")
  cat(res, "\n", sep = "")
}


class(seguid) <- c("cli_function", class(seguid))
attr(seguid, "cli") <- function(..., type = c("seguid", "lsseguid", "csseguid", "ldseguid", "cdseguid")) {
  type <- match.arg(type)
  cli_call_fcn(..., fcn = type)
}

class(lsseguid) <- c("cli_function", class(lsseguid))
attr(lsseguid, "cli") <- function(...) {
  cli_call_fcn(..., fcn = lsseguid)
}

class(csseguid) <- c("cli_function", class(csseguid))
attr(csseguid, "cli") <- function(...) {
  cli_call_fcn(..., fcn = csseguid)
}

class(ldseguid) <- c("cli_function", class(ldseguid))
attr(ldseguid, "cli") <- function(...) {
  cli_call_fcn(..., fcn = ldseguid)
}

class(cdseguid) <- c("cli_function", class(cdseguid))
attr(cdseguid, "cli") <- function(...) {
  cli_call_fcn(..., fcn = cdseguid)
}
