cli_help_string <- '
{{ package }}: {{ title }}

Usage:
 Rscript -e seguid::seguid [options]

Options:  
 --help             Display the full help page with examples
 --version          Output version of this software

Examples:
Rscript -e port4me::port4me --version
Rscript -e port4me::port4me --help

Version: {{ version }}
Copyright: Henrik Bengtsson (2023)
License: MIT
'

#' @importFrom utils file_test
cli_call_fcn <- function(..., file = NULL, debug = FALSE, fcn) {
  if (is.character(fcn)) {
    fcn <- get(fcn, mode = "function", envir = getNamespace(.packageName), inherits = FALSE)
  }
  stopifnot(length(debug) == 1, is.logical(debug), !is.na(debug))
  
  seq <- NULL
  
  args <- list(...)
  nargs <- length(args)
  if (nargs > 0) {
    names <- names(args)
    if (is.null(names)) names <- rep("", times = nargs) 
    ## At most one unnamed CLI option
    unnamed <- which(nchar(names) == 0)
    if (length(unnamed) > 1) {
      unknown <- unnamed[-length(unnamed)]
      stop("Unknown option: ", paste(sQuote(unknown), collapse = ", "))
    }
    seq <- args[[unnamed]]
    if (!is.null(file)) {
      stop("Option --file=<filename> already specified")
    }
  }

  if (is.null(seq)) {
    ## Get sequence input from file?
    if (!is.null(file)) {
      if (!file_test("-f", file)) {
        stop("No such file: ", sQuote(file))
      }
      seq <- readLines(file)
    } else {
      ## Get arguments from the standard input
      seq <- readLines("stdin")
    }
  }
  seq <- paste(seq, collapse = "\n")

  if (debug) {
    message(sprintf("Sequence data:\n%s", seq))
  }
  
  argnames <- names(formals(fcn))
  if (is.element("crick", argnames)) {
    args <- tuple_from_repr(seq)
    if (debug) {
      msg <- sprintf("Sequence tuple:\nwatson=%s\ncrick=%s", args[[1]], args[[2]])
    }
    if (!is.element("overhang", argnames)) {
      args <- args[-3]
    } else if (debug) {
      msg <- sprintf("%s\noverhang=%d", msg, args[[3]])
    }
    if (debug) message(msg)
  } else {
    args <- list(seq)
  }

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


