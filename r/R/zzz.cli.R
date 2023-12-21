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

echo "ACGT" | Rscript -e seguid::slseguid
Rscript -e seguid::slseguid <<< "ACGT"
Rscript -e seguid::dlseguid <<< $\'-CGT\nTGCA\'
Rscript -e seguid::dcseguid <<< $\'ACGT\nTGCA\'

Version: {{ version }}
Copyright: Henrik Bengtsson (2023)
License: MIT
'

#' @importFrom utils file_test
cli_call_fcn <- function(..., table = "dna", file = NULL, debug = FALSE, fcn) {
  if (is.character(fcn)) {
    fcn <- get(fcn, mode = "function", envir = getNamespace(.packageName), inherits = FALSE)
  }
  stopifnot(length(debug) == 1, is.logical(debug), !is.na(debug))

  table2 <- get_table(table)
  
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
    message(sprintf("Arguments:\n%s", paste(utils::capture.output(str(args)), collapse = "\n")))
  }

  argnames <- names(formals(fcn))
  if (is.element("crick", argnames)) {
    nseq <- length(strsplit(seq, split = "\n", fixed = TRUE)[[1]])
    if (nseq == 1) {
      args2 <- list(watson = seq, crick = rc(seq, table = table2))
    } else {
      args2 <- tuple_from_repr(seq, table = table2)
    }
    if (debug) {
      msg <- sprintf("Sequence tuple:\nwatson=%s\ncrick=%s", args2[[1]], args2[[2]])
    }
    if (!is.element("overhang", argnames)) {
      args2 <- args2[-3]
    } else if (debug) {
      msg <- sprintf("%s\noverhang=%d", msg, args2[[3]])
    }
    if (debug) message(msg)
  } else {
    args2 <- list(seq)
  }
  args <- c(args2, args, table = table)

  res <- do.call(fcn, args = args)
  
  cat(res)
  cat("\n")
}


class(seguid) <- c("cli_function", class(seguid))
attr(seguid, "cli") <- function(..., type = c("seguid", "slseguid", "scseguid", "dlseguid", "dcseguid")) {
  type <- match.arg(type)
  cli_call_fcn(..., fcn = type)
}

class(slseguid) <- c("cli_function", class(slseguid))
attr(slseguid, "cli") <- function(...) {
  cli_call_fcn(..., fcn = slseguid)
}

class(scseguid) <- c("cli_function", class(scseguid))
attr(scseguid, "cli") <- function(...) {
  cli_call_fcn(..., fcn = scseguid)
}

class(dlseguid) <- c("cli_function", class(dlseguid))
attr(dlseguid, "cli") <- function(...) {
  cli_call_fcn(..., fcn = dlseguid)
}

class(dcseguid) <- c("cli_function", class(dcseguid))
attr(dcseguid, "cli") <- function(...) {
  cli_call_fcn(..., fcn = dcseguid)
}
