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

cli_call_fcn <- function(fcn, ...) {
  if (is.character(fcn)) {
    fcn <- get(fcn, mode = "function", envir = getNamespace(.packageName), inherits = FALSE)
  }

  ## Get arguments
  seq <- readLines('stdin')

  argnames <- names(formals(fcn))
  if (is.element("crick", argnames)) {
    seq <- paste(seq, collapse = "\n")
    args <- tuple_from_repr(seq)
    if (!is.element("overhang", argnames)) {
      args <- args[-3]
    }
  } else {
    args <- list(seq)
  }

  res <- do.call(fcn, args = args)
  
  cat(res)
  cat("\n")
}


class(seguid) <- c("cli_function", class(seguid))
attr(seguid, "cli") <- function(type = c("seguid", "slseguid", "scseguid", "dlseguid", "dcseguid"), ...) {
  type <- match.arg(type)
  cli_call_fcn(type, ...)
}

class(slseguid) <- c("cli_function", class(slseguid))
attr(slseguid, "cli") <- function(...) {
  cli_call_fcn(slseguid, ...)
}

class(scseguid) <- c("cli_function", class(scseguid))
attr(scseguid, "cli") <- function(...) {
  cli_call_fcn(scseguid, ...)
}

class(dlseguid) <- c("cli_function", class(dlseguid))
attr(dlseguid, "cli") <- function(...) {
  cli_call_fcn(dlseguid, ...)
}

class(dcseguid) <- c("cli_function", class(dcseguid))
attr(dcseguid, "cli") <- function(...) {
  cli_call_fcn(dcseguid, ...)
}


