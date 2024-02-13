parse_sequence_string <- seguid:::parse_sequence_string


spec <- "TATGCC\nATACGG"
truth <- list(type = "ds", specification = spec, watson = "TATGCC", crick = "GGCATA")
res <- parse_sequence_string(spec)
str(res)
stopifnot(all.equal(res, truth))

spec <- "-TATGCC\nCATACG-"
truth <- list(type = "ds", specification = spec, watson = "-TATGCC", crick = "-GCATAC")
res <- parse_sequence_string(spec)
str(res)
stopifnot(all.equal(res, truth))

spec <- "TATGCC--\n-TACGGGG"
truth <- list(type = "ds", specification = spec, watson = "TATGCC--", crick = "GGGGCAT-")
res <- parse_sequence_string(spec)
str(res)
stopifnot(all.equal(res, truth))

spec <- "TATGCC--\n--ACGGGG"
truth <- list(type = "ds", specification = spec, watson = "TATGCC--", crick = "GGGGCA--")
res <- parse_sequence_string(spec)
str(res)
stopifnot(all.equal(res, truth))


spec <- "- TGCC\nATACGG"
res <- tryCatch({
  parse_sequence_string(spec)
}, error = identity)
print(res)
stopifnot(inherits(res, "error"))

spec <- " -TGCC\nATACGG"
res <- tryCatch({
  parse_sequence_string(spec)
}, error = identity)
print(res)
stopifnot(inherits(res, "error"))

spec <- "---TGCC\nATACGG"
res <- tryCatch({
  parse_sequence_string(spec)
}, error = identity)
print(res)
stopifnot(inherits(res, "error"))

spec <- "---TGCC\n-ATACGG"
res <- tryCatch({
  parse_sequence_string(spec)
}, error = identity)
print(res)
stopifnot(inherits(res, "error"))

spec <- "TATGCC-\nATACGG-"
res <- tryCatch({
  parse_sequence_string(spec)
}, error = identity)
print(res)
stopifnot(inherits(res, "error"))
