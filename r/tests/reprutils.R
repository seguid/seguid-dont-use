tuple_from_repr <- seguid:::tuple_from_repr
# repr_from_tuple <- seguid:::repr_from_tuple

assert_equal_tuple <- function(tuple, truth) {
  res <- all.equal(tuple, truth, check.attributes = FALSE)
  if (!isTRUE(res)) {
    str(list(tuple = tuple, truth = truth))
    stopifnot(res)
  }
}


truth <- list("TATGCC", "GCATAC", 1L)
rpr <-
"             \
      -TATGCC \
      CATACG- \
"
tuple <- tuple_from_repr(rpr)
assert_equal_tuple(tuple, truth)


truth <- list("TATGCC", "GGGGCAT", -1L)
rpr <-
"               \
       TATGCC-- \
       -TACGGGG \
"
tuple <- tuple_from_repr(rpr)
assert_equal_tuple(tuple, truth)


rpr <-
" # Source code comments like this one are not allowed \
        TATGCC \
        ATACGG \
"
res <- tryCatch({
  tuple_from_repr(rpr)
}, error = identity)
stopifnot(inherits(res, "error"))


rpr <-
"               \
       TATGCC-- \
       --ACGGGG \
"
tuple <- tuple_from_repr(rpr)
assert_equal_tuple(tuple, list("TATGCC", "GGGGCA", -2L))


rpr <-
"             \
       TATGCC \
       ATACGG \
"
tuple <- tuple_from_repr(rpr)
assert_equal_tuple(tuple_from_repr(rpr), list('TATGCC', 'GGCATA', 0L))


rpr <-
"             \
       - TGCC \
       ATACGG \
"
res <- tryCatch({
  tuple_from_repr(rpr)
}, error = identity)
stopifnot(inherits(res, "error"))

        
rpr <-
"             \
        -TGCC \
       ATACGG \
"
res <- tryCatch({
  tuple_from_repr(rpr)
}, error = identity)
stopifnot(inherits(res, "error"))

        
rpr <-
"             \
      ---TGCC \
       ATACGG \
"
res <- tryCatch({
  tuple_from_repr(rpr)
}, error = identity)
stopifnot(inherits(res, "error"))


rpr <-
"             \
      ---TGCC \
      -ATACGG \
"
res <- tryCatch({
  tuple_from_repr(rpr)
}, error = identity)
stopifnot(inherits(res, "error"))


rpr <-
"              \
       TATGCC- \
       ATACGG- \
"
res <- tryCatch({
  tuple_from_repr(rpr)
}, error = identity)
stopifnot(inherits(res, "error"))


#assert_equal_tuple(repr_from_tuple(*("TATGCC", "GGGGCA", -2)) == "TATGCC--\n--ACGGGG"
