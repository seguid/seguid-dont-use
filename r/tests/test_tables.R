get_alphabet <- seguid:::get_alphabet
make_alphabet <- seguid:::make_alphabet

alphabet <- get_alphabet("{DNA}")
alphabet <- get_alphabet("{RNA}")
alphabet <- get_alphabet("{IUPAC}")
alphabet <- get_alphabet("{protein}")

## Unknown alphabet
res <- tryCatch({
  alphabet <- get_alphabet("unknown")
}, error = identity)
stopifnot(inherits(res, "error"))


truth <- get_alphabet("{DNA}")
alphabet <- make_alphabet("AT,TA,CG,GC")
stopifnot(identical(sort(alphabet), sort(truth)))

truth <- get_alphabet("{RNA}")
alphabet <- make_alphabet("AU,UA,CG,GC")
stopifnot(identical(sort(alphabet), sort(truth)))

truth <- get_alphabet("{IUPAC}")
alphabet <- make_alphabet("AT,BV,CG,DH,GC,HD,KM,MK,NN,SS,TA,VB,WW")
stopifnot(identical(sort(alphabet), sort(truth)))

truth <- get_alphabet("{protein}")
alphabet <- make_alphabet("A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y")
stopifnot(identical(sort(alphabet), sort(truth)))

## Missing complementary component
res <- tryCatch({
  alphabet <- make_alphabet("AT")
}, error = identity)
stopifnot(inherits(res, "error"))

## Invalid components; should be either one or two characters
res <- tryCatch({
  alphabet <- make_alphabet("ATT")
}, error = identity)
stopifnot(inherits(res, "error"))

## Incompatible components
res <- tryCatch({
  alphabet <- make_alphabet("AT,T")
}, error = identity)
stopifnot(inherits(res, "error"))

truth <- c(get_alphabet("{DNA}"), a = "u", u = "a")
alphabet <- get_alphabet("{DNA},au,ua")
stopifnot(identical(sort(alphabet), sort(truth)))
