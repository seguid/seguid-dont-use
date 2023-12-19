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

main <- function(type = c("seguid", "slseguid", "scseguid", "dlseguid", "dcseguid"), ...) {
  type <- match.arg(type)
  fcn <- get(type, mode = "function", envir = getNamespace(.packageName), inherits = FALSE)

  ## Get arguments
  seq <- readLines('stdin')

  if (grepl("^d", type)) {
    seq <- paste(seq, collapse = "\n")
    args <- tuple_from_repr(seq)
    if (!is.element("overhang", names(formals(fcn)))) {
      args <- args[-3]
    }
  } else {
    args <- seq
  }

  res <- do.call(fcn, args = args)
  
  cat(res)
  cat("\n")
}

class(seguid) <- c("cli_function", class(seguid))
attr(seguid, "cli") <- main
