parse_cli_args <- function() {
  ## Parse command-line arguments
  cli_args <- commandArgs(trailingOnly = TRUE)

  more_opts <- TRUE
  args <- list()
  while (length(cli_args) > 0) {
    arg <- cli_args[1]
    if (more_opts) {
      if (grepl(pattern <- "^--args$", arg)) {
        ## Ignore --args
      } else if (grepl(pattern <- "^--$", arg)) {
        more_opts <- FALSE
      } else if (grepl(pattern <- "^--([[:alnum:]]+)=(.*)$", arg)) {
        name <- gsub(pattern, "\\1", arg)
        value <- gsub(pattern, "\\2", arg)
        if (grepl("^[+-]?[[:digit:]]+$", value)) {
          value_int <- suppressWarnings(as.integer(value))
          if (!is.na(value_int)) value <- value_int
        }
        args[[name]] <- value
      } else if (grepl(pattern <- "^--([[:alnum:]]+)$", arg)) {
        name <- gsub(pattern, "\\1", arg)
        args[[name]] <- I(TRUE)
      } else if (grepl(pattern <- "^--", arg)) {
        stop(sprintf("Unknown command-line argument: %s", arg))
      } else {
        args[[length(args) + 1L]] <- arg
      }
    } else {
      args[[length(args) + 1L]] <- arg
    }
    cli_args <- cli_args[-1]
  }

  args
}


#' @importFrom utils packageDescription packageVersion
#' @export
print.cli_function <- function(x, ..., envir = parent.frame()) {
  if (interactive()) {
    attr(x, "cli") <- NULL
    class(x) <- setdiff(class(x), "cli_function")
    return(NextMethod())
  }

  args <- parse_cli_args()

  if (isTRUE(args$version)) {
    cat(as.character(packageVersion(.packageName)), "\n", sep = "")
  } else if (isTRUE(args$help)) {
    msg <- cli_help_string
    msg <- sub("{{ package }}", .packageName, msg, fixed = TRUE)
    msg <- sub("{{ title }}", packageDescription(.packageName)[["Title"]], msg, fixed = TRUE)
    msg <- sub("{{ version }}", packageVersion(.packageName), msg, fixed = TRUE)
    cat(msg)
  } else {
    ## Is there a custom function?
    fcn <- attr(x, "cli", exact = TRUE)
    if (is.null(fcn)) fcn <- x

    res <- withVisible(do.call(fcn, args = args, envir = envir))

    # Should the result be printed?
    if (res$visible) {
      value <- res$value
      if (is.integer(value)) {
        cat(sprintf("%i\n", value), collapse = "", sep = "")
      } else if (is.logical(value)) {
        quit(save = "no", status = as.integer(!value))
      }
    }
  }
  invisible()
}
