library(seguid)

tuple_from_representation <- seguid:::tuple_from_representation

# ----------------------------------------------------------------
# Test 1
# ----------------------------------------------------------------
truth <- list(watson = "TATGCC", crick = "gcatac", overhang = 1L)

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

tuples <- lapply(rprs, FUN = function(rpr) {
  cat(rpr, sep = "")
  tuple <- tuple_from_representation(rpr)
  str(tuple)
  stopifnot(identical(tuple, truth))
})




# ----------------------------------------------------------------
# Test 2
# ----------------------------------------------------------------
cat("--------------------------------\n")
truth <- list(watson = "TATGCC", crick = "ggcata", overhang = 0L)

rprs <- list()

rprs[[1]] <- "   \
                 \
    5'-TATGCC-3' \
       ||||||    \
    3'-atacgg-5' \
                 \
"

rprs[[2]] <- "   \
                 \
    5'-TATGCC-3' \
       ||||||    \
    3'-atacgg-5' \
                 \
"

rprs[[3]] <- "   \
       TATGCC    \
       ||||||    \
       atacgg    \
"

rprs[[4]] <- "   \
   TATGCC        \
   atacgg        \
                 \
"

tuples <- lapply(rprs, FUN = function(rpr) {
  cat(rpr, sep = "")
  tuple <- tuple_from_representation(rpr)
  str(tuple)
  stopifnot(identical(tuple, truth))
})



# ----------------------------------------------------------------
# Test 3
# ----------------------------------------------------------------
cat("--------------------------------\n")
truth <- list(watson = "TATGCC", crick = "ggggca", overhang = -1L)

rprs <- list()

rprs[[1]] <- "   \
                 \
  5'-TATGCC-3'   \
      |||||      \
   3'-acgggg-5'  \
                 \
"

rprs[[2]] <- "   \
                 \
  5'-TATGCC-3'   \
      |||||      \
   3'-acgggg-5'  \
                 \
"

rprs[[3]] <- "   \
     TATGCC      \
      |||||      \
      acgggg     \
"

rprs[[4]] <- "   \
     TATGCC      \
      acgggg     \
                 \
"

tuples <- lapply(rprs, FUN = function(rpr) {
  cat(rpr, sep = "")
  tuple <- tuple_from_representation(rpr)
  str(tuple)
  stopifnot(identical(tuple, truth))
})
