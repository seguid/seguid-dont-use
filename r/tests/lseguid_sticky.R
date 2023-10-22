library(seguid)
tuple_from_representation <- seguid:::tuple_from_representation

# ------------------------------------------------------------
# Case 1
# ------------------------------------------------------------
rpr <- "\
   TATGCC \
  catacg  \
"
truth <- "Jv9Z9JJ0IYnG-dTPBGwhDyAqnmU"
tuple <- tuple_from_representation(rpr)
stopifnot(
  do.call(lseguid_sticky, tuple) == truth,
  useguid(toupper(" gcatac\nCCGTAT")) == truth
)


# ------------------------------------------------------------
# Case 2
# ------------------------------------------------------------
rpr <- "\
  gTATGC  \
  catacg  \
"
truth <- "b0Xa5pLe4LNd5T8fhGWHicCI_f4"
tuple <- tuple_from_representation(rpr)
stopifnot(
  do.call(lseguid_sticky, tuple) == truth,
  lseguid_blunt("gcatac") == truth
)

# ------------------------------------------------------------
# Case 3
# ------------------------------------------------------------
rpr <- "\
  gTATGC   \
   atacgc  \
"
truth <- "3PZNLU2tPDYs78AJP3mg5w1uEw4"
tuple <- tuple_from_representation(rpr)
stopifnot(
  do.call(lseguid_sticky, tuple) == truth,
  useguid("cgcata\n CGTATg") == truth
)

