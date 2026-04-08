buildStoredStates <- function() {
  Rcpp::sourceCpp(here::here("data-raw", "Stirling2FlintARB.cpp"))
  logs <- lapply(seq(1e3, 5e4, 1e3), logstirling2_raw_flintARB)
  usethis::use_data(logs, internal = TRUE, overwrite = TRUE)
}
