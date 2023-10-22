library(seguid)

tuple_from_representation <- seguid:::tuple_from_representation

rprs <- list()

rprs[[1]] <- "   \
                 \
    5'-TATGCC-3' \
       |||||     \
   3'-catacg-5'  \
                 \
"

rprs[[2]] <- "   \
                 \
    5'-TATGCC-3' \
       |||||     \
   3'-catacg-5'  \
                 \
"

rprs[[3]] <- "   \
       TATGCC    \
       |||||     \
      catacg     \
"

rprs[[4]] <- "   \
   TATGCC        \
  catacg         \
                 \
"

tuples <- lapply(rprs, FUN = tuple_from_representation)
void <- mapply(rprs, tuples, FUN = function(rpr, tuple) {
  cat(rpr, sep = "")
  str(tuple)
  stopifnot(identical(tuple, tuples[[1]]))
})


