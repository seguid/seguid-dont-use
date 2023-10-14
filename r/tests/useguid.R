library(seguid)

stopifnot(
  useguid("aaa") == "YG7G6b2Kj_KtFOX63j8mRHHoIlE",
  useguid("aAA") == "YG7G6b2Kj_KtFOX63j8mRHHoIlE"
)

