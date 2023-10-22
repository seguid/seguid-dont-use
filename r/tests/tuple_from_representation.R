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

truth <- list(watson = "TATGCC", crick = "gcatac", overhang = 1L)

tuples <- lapply(rprs, FUN = function(rpr) {
  cat(rpr, sep = "")
  tuple <- tuple_from_representation(rpr)
  str(tuple)
  stopifnot(identical(tuple, truth))
})


