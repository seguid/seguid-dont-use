library(seguid)

stopifnot(
  seguid("ACGTACGTACGT") == "If6HIvcnRSQDVNiAoefAzySc6i4",
  seguid("acgtACGTacgt") == "If6HIvcnRSQDVNiAoefAzySc6i4"
)
