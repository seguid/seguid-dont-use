get_table <- seguid:::get_table
make_table <- seguid:::make_table

table <- get_table("{DNA}")
table <- get_table("{RNA}")
table <- get_table("{IUPAC}")
table <- get_table("{protein}")

## Unknown table
res <- tryCatch({
  table <- get_table("unknown")
}, error = identity)
stopifnot(inherits(res, "error"))


truth <- get_table("{DNA}")
table <- make_table("AT,TA,CG,GC")
stopifnot(identical(sort(table), sort(truth)))

truth <- get_table("{RNA}")
table <- make_table("AU,UA,CG,GC")
stopifnot(identical(sort(table), sort(truth)))

truth <- get_table("{IUPAC}")
table <- make_table("AT,BV,CG,DH,GC,HD,KM,MK,NN,SS,TA,VB,WW")
stopifnot(identical(sort(table), sort(truth)))

truth <- get_table("{protein}")
table <- make_table("A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y")
stopifnot(identical(sort(table), sort(truth)))

## Missing complementary component
res <- tryCatch({
  table <- make_table("AT")
}, error = identity)
stopifnot(inherits(res, "error"))

## Invalid components; should be either one or two characters
res <- tryCatch({
  table <- make_table("ATT")
}, error = identity)
stopifnot(inherits(res, "error"))

## Incompatible components
res <- tryCatch({
  table <- make_table("AT,T")
}, error = identity)
stopifnot(inherits(res, "error"))

## Duplicates
res <- tryCatch({
  table <- make_table("AT,AT,TA,TA")
}, error = identity)
stopifnot(inherits(res, "error"))


truth <- c(get_table("{DNA}"), a = "u", u = "a")
table <- get_table("{DNA},au,ua")
stopifnot(identical(sort(table), sort(truth)))
