# job::job({
  # path <-
  #   "C:/Users/jon.blood/OneDrive - US Army/H Drive/Miscellaneous/R/Stirling2"
#   Rcpp::sourceCpp(file.path(path, "Stirling2FlintARB.cpp"))
#     logs <- lapply(seq(1e3, 5e4, 1e3), \(n) {cat(n, Sys.time(), "\n");
#       logstirling2_raw_flintARB(n)})
#     logs <- unlist(lapply(logs, \(x)
#                           `dim<-`(`dim<-`(x, c(16, length(x)/16))[1:10,], NULL)))
#     fst::write.fst(data.table::data.table(s = logs),
#                    file.path(path, "logS3.fst"), 100)
#     job::export(NULL)
# }, import = NULL)
#  
# job::job({
#   library(future.apply)
#   # plan(multisession, workers = parallelly::availableCores() - 1)
#   plan(sequential, split = TRUE)
#   path <- "C:/Users/jon.blood/OneDrive - US Army/H Drive/Miscellaneous/R/Stirling2"
# 
#   # 1. Define a cache directory for the compiled DLL
#   cpp_file <- file.path(path, "Stirling2FlintARB.cpp")
#   cache_dir <- file.path(path, "rcpp_cache")
#   dir.create(cache_dir, showWarnings = FALSE)
# 
#   # 2. Compile in the main session FIRST. This builds the DLL and saves it to the cache.
#   Rcpp::sourceCpp(cpp_file, cacheDir = cache_dir)
# 
#   # 3. Run the parallel loop
#   logs <- future.apply::future_lapply(seq(1e3, 5e4, 1e3), function(n) {
#     # The worker looks in the cache, skips compilation, and loads the DLL into its memory
#     Rcpp::sourceCpp(cpp_file, cacheDir = cache_dir)
# 
#     # Now the pointer is valid in this specific background process!
#     logstirling2_raw_flintARB(n)
#   }, future.seed = NULL)
# 
#   # 4. Process and save
#   logs <- unlist(lapply(logs, \(x)
#                         `dim<-`(`dim<-`(x, c(16, length(x)/16))[1:10,], NULL)))
# 
#   fst::write.fst(data.table::data.table(s = logs),
#                  file.path(path, "logS2.fst"), 100)
#   job::export(NULL)
# }, import = NULL)

buildStoredStates <- function() {
  Rcpp::sourceCpp(here::here("Stirling2FlintARB.cpp"))
  logs <- future.apply::future_lapply(seq(1e3, 5e4, 1e3),
                                      logstirling2_raw_flintARB)
  logs <- unlist(lapply(logs, \(x)
                        `dim<-`(`dim<-`(x, c(16, length(x)/16))[1:10,], NULL)))
  saveRDS(logs, here::here("logS.rds"))
}

readStoredStates <- function(f) {
  x <- readRDS(f)
  x <- rbind(`dim<-`(x, c(10, length(x)/10)),
             `dim<-`(vector("raw", 0.6*length(x)), c(6, length(x)/10)))
  dim(x) <- NULL
  unname(split(x, rep.int(1:50, seq(15984, 799984, 16e3))))
}

#' Stirling Numbers of the Second Kind
#'
#' @param n Set size (length-1, row).
#' @param k Subset size (length-1, column).
#'
#' @export
stirling2direct <- function(n, k) {
  # replacement for buggy gmp::Stirling2(method = "direct")
  gmp::as.bigz(sum(gmp::chooseZ(-2:-(k + 1), (k - 1):0)*
                     gmp::as.bigz(1:k)^n)/gmp::factorialZ(k))
}

#' log-Stirling Numbers of the Second Kind
#'
#' @param n Vector of set sizes (rows).
#' @param k Vector of subset sizes (columns). If NULL, returns full rows.
#' @param as.matrix Logical; if TRUE, returns a matrix with `length(n)` rows and
#'   `length(k)` columns (or `max(n)` columns if `k` is NULL). If FALSE, -Inf
#'   values (k > n) are omitted.
#' @param ones Logical; if FALSE, excludes `k = 1` and `k = n` from the output
#'   (set to TRUE if `as.matrix == TRUE`, `!is.null(k)`, or `any(n < 3)`).
#' @param parallel Logical; if TRUE, uses future.apply for parallel processing
#'   across blocks of `n`.
#'
#' @export
logStirling2 <- function(n, k = NULL, as.matrix = TRUE, ones = TRUE,
                      parallel = FALSE) {
  # check inputs
  if (length(i <- which(n < 1))) {
    n <- n[-i]
    stopifnot("All n are less than 1." = length(n) > 0)
    warning("n < 1 dropped.")
  }
  
  if (any(n%%1 != 0)) {
    warning("n coerced to natural numbers.")
    n <- floor(n)
  }
  
  if (!is.null(k)) {
    if (length(i <- which(k < 1))) {
      k <- k[-i]
      stopifnot("All k are less than 1." = length(k) > 0)
      warning("k < 1 dropped.")
    }
    
    if (any(k%%1 != 0)) {
      warning("k coerced to natural numbers.")
      k <- floor(k)
    }
    
    if (length(n) == 1 && length(k) == 1) { # use `stirling2direct`
      S <- log(stirling2direct(n, k))
      return(if (as.matrix) matrix(S, 1, 1, 0, list(n, k)) else S)
    }
  }
  
  # Group n by cache blocks (0-999, 1000-1999, etc.)
  nu <- sort(unique(n))
  if (nu[length(nu)] < 3) {
    if (is.null(k)) k <- 1:nu[length(nu)]
    S <- lapply(n, \(i) numeric(i)[k])
    if (as.matrix) {
      S <- `dimnames<-`(do.call(rbind, S), list(n, k))
    } else {
      S <- unlist(S)
      S <- S[!is.na(S)]
    }
    return(S)
  }
  
  if (length(i <- which(nu < 3))) {
    nu0 <- nu[i]
    nu <- nu[-i]
    ones <- TRUE
  } else nu0 <- vector(mode(nu), 0)
  
  blocks <- split(nu, nu%/%1e3)
  states <- c(list(NULL), get("logS_states", envir = pkg_env))
  bu <- pmin(as.integer(names(blocks)) + 1, length(states))
  sb <- states[bu]
  ni <- match(n, c(nu0, nu))
  
  # Set up execution wrapper
  fapply <- if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_lapply
  } else base::lapply
  
  # Process blocks
  S <- unlist(fapply(1:length(blocks), function(i) {
    if (length(blocks[[i]]) == 1) { # single row
      logStirling2Row_C(blocks[[i]], state = sb[[i]])
    } else if (blocks[[i]][1] == (bu[i] - 1)*1e3 + (bu[i] == 1)*3 &&
               all(diff(blocks[[i]]) == 1)) { # all rows
      logStirling2All_C(blocks[[i]][length(blocks[[i]])], state = sb[[i]])
    } else logStirling2Mult_C(blocks[[i]], state = sb[[i]])
  }), FALSE, FALSE)
  
  S <- unname(split(S, rep.int(nu, nu - 2)))
  
  if (is.null(k)) {
    if (!as.matrix && !ones) return(unlist(S[ni])) else k <- 1:nu[length(nu)]
  }
  
  S <- c(lapply(nu0, \(i) numeric(i)[k]), lapply(S, \(s) c(0, s, 0)[k]))[ni]
  
  if (as.matrix) {
    `dimnames<-`(do.call(rbind, S), list(c(nu0, nu)[ni], k))
  } else unlist(lapply(S, \(s) s[!is.na(s)]), FALSE, FALSE)
}

#' log-Stirling Numbers of the Second Kind (Temme approximation)
#'
#' @param n Vector of set sizes (rows).
#' @param k Vector of subset sizes (columns). If NULL, returns full rows.
#' @param as.matrix Logical; if TRUE, returns a matrix with `length(n)` rows and
#'   `length(k)` columns (or `max(n)` columns if `k` is NULL). If FALSE, -Inf
#'   values (k > n) are omitted.
#' @param ones Logical; if FALSE, excludes `k = 1` and `k = n` from the output
#'   (set to TRUE if `as.matrix == TRUE`, `!is.null(k)`, or `any(n < 3)`).
#'
#' @export
logStirling2Temme <- function(n, k = NULL, as.matrix = TRUE, ones = TRUE) {
  # Returns log(S(n, k)) using Temme's approximation for
  # k = 4:(n - 2) and exact formulas for k = c(2, 3, n - 2, n - 1).
  # https://core.ac.uk/download/pdf/301651745.pdf
  # https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind#Asymptotic_approximation
  
  # check inputs
  if (length(i <- which(n < 1))) {
    n <- n[-i]
    stopifnot("All n are less than 1." = length(n) > 0)
    warning("n < 1 dropped.")
  }
  
  if (any(n%%1 != 0)) {
    warning("n coerced to natural numbers.")
    n <- floor(n)
  }
  
  rn <- range(n)
  
  if (rn[2] < 3) {
    if (is.null(k)) k <- 1:rn[2]
    S <- lapply(n, \(i) numeric(i)[k])
    
    if (as.matrix) {
      S <- `dimnames<-`(do.call(rbind, S), list(n, k))
    } else {
      S <- unlist(S)
      S <- S[!is.na(S)]
    }
    
    return(S)
  }
  
  if (!is.null(k)) {
    if (length(i <- which(k < 1))) {
      k <- k[-i]
      stopifnot("All k are less than 1." = length(k) > 0)
      warning("k < 1 dropped.")
    }
    
    if (any(k%%1 != 0)) {
      warning("k coerced to natural numbers.")
      k <- floor(k)
    }
    
    ones <- TRUE
    rk <- range(k)
    stopifnot(rk[2] < rn[2] || rk[1] > 0)
  } else {
    if (as.matrix) ones <- TRUE
    rk <- c(2 - ones, rn[2] - 1 + ones)
    k <- rk[1]:rk[2]
  }
  
  if (as.matrix) {
    k0 <- k
    n0 <- n
  }
  
  nk <- length(k)
  nn <- length(n)
  
  if (as.matrix) {
    n <- rep.int(n, nk)
    k <- rep(k, each = nn)
  } else {
    n <- rep(n, each = nk)
    k <- rep.int(k, nn)
    if (ones) i <- k <= n else i <- k < n
    n <- n[i]
    k <- k[i]
  }
  
  s <- numeric(length(n))
  n1 <- n - 1
  i <- k != 1 & k != n
  
  if (as.matrix && rk[2] > rn[1]) {
    j <- which(k > n)
    i[j] <- FALSE
    s[j] <- NA_real_
  }
  # exact formula for edges
  j <- which(k == 2 & i)
  i[j] <- FALSE
  s[j] <- n1[j]*log(2) + log1p(-2^-n1[j])
  j <- which(k == 3 & i)
  i[j] <- FALSE
  s[j] <- n1[j]*log(3) + log1p(3^-n1[j] - 2*(2/3)^n1[j]) - log(2)
  j <- which(n - k == 2 & i)
  i[j] <- FALSE
  s[j] <- log1p(3*n[j] - 6) + lchoose(n[j], 3) - log(4)
  j <- which(n - k == 1 & i)
  i[j] <- FALSE
  s[j] <- log(n1[j]) + log1p(n1[j]) - log(2)
  
  if (any(i)) {
    n <- n[i]
    k <- k[i]
    v <- n/k
    G <- lamW::lambertW0(v/exp(v))
    G1 <- 1 - G
    s[i] <- ((logv1 <- log(v - 1)) - log(v*G1))/2 + lchoose(n, k) +
      (n - k)*(logv1 - log(v - G)) + n*log(k) + k*(G1 - log(n))
  }
  
  if (as.matrix) {
    dim(s) <- c(nn, nk)
    dimnames(s) <- list(n0, k0)
  }
  
  s
}

pkg_env <- new.env(parent = emptyenv())
Rcpp::sourceCpp("C:/Users/jon.blood/OneDrive - US Army/H Drive/Miscellaneous/R/Stirling2/Stirling2.cpp")
local_rds_path <- "C:/Users/jon.blood/OneDrive - US Army/H Drive/Miscellaneous/R/Stirling2/logS.rds"
assign("logS_states", readStoredStates(local_rds_path), envir = pkg_env)
n <- sample(5e3, 20, TRUE)

library(testthat)

test_that("Input sanitation and warnings work", {
  # Decimals should be dropped with a warning
  expect_warning(res <- logStirling2(n = c(5, 5.5)), "n coerced to natural numbers.")
  expect_equal(rownames(res), c("5", "5"))
  
  # Zero and negatives should be dropped
  expect_warning(logStirling2(n = c(10, 0, -5)), "n < 1 dropped.")
  
  # All invalid inputs should trigger the stopifnot
  expect_error(suppressWarnings(logStirling2(n = c(-1, 0.5))), "All n are less than 1.")
  expect_error(suppressWarnings(logStirling2(n = 5, k = 0)), "All k are less than 1.")
})

test_that("Single value dispatch (stirling2direct bypass) is correct", {
  # n = 10, k = 5
  res_mat <- logStirling2(10, 5, as.matrix = TRUE)
  res_vec <- logStirling2(10, 5, as.matrix = FALSE)
  
  expect_true(is.matrix(res_mat))
  expect_false(is.matrix(res_vec))
  expect_equal(as.numeric(res_mat), res_vec)
})

test_that("Small n (n < 3) bypass logic is correct", {
  # Tests the nu0 logic and dimension matching
  res <- logStirling2(n = 1:5, k = 1:2, as.matrix = TRUE)
  expect_equal(dim(res), c(5, 2))
  
  # Log(S(1,1)) and Log(S(2,1)) are 0
  expect_equal(unname(res[1:2, 1]), c(0, 0)) 
})

test_that("Backend C++ Routing Dispatches Correctly", {
  # We assume standard correctness, but we want to ensure the R logic 
  # doesn't crash when hitting the different branches.
  
  # 1. Row_C (Single row in a cache block)
  expect_no_error(logStirling2(n = 1050, k = 5))
  
  # 2. All_C (Dense sequence crossing a block boundary)
  # 1000:1005 should trigger All_C for block 1000
  expect_no_error(res_all <- logStirling2(n = 998:1002, as.matrix = TRUE))
  expect_equal(nrow(res_all), 5)
  
  # 3. Mult_C (Sparse sequence in a block)
  expect_no_error(res_mult <- logStirling2(n = c(1010, 1050, 1090), as.matrix = TRUE))
  expect_equal(nrow(res_mult), 3)
})

test_that("Formatting flags shape the output correctly", {
  # k = NULL tests
  res_null_mat <- logStirling2(n = c(4, 5), k = NULL, as.matrix = TRUE)
  expect_equal(dim(res_null_mat), c(2, 5)) # max(n) columns
  
  # ones = FALSE test
  res_no_ones <- logStirling2(n = 5, k = NULL, as.matrix = FALSE, ones = FALSE)
  # S(5,1) and S(5,5) should be removed. Length should be 3.
  expect_equal(length(res_no_ones), 3)
})

test_that("Parallel execution matches sequential execution", {
  # Requires future setup
  future::plan(future::multisession, workers = 2)
  
  n_test <- c(50, 1050, 2050) # Crosses multiple cache blocks
  
  res_seq <- logStirling2(n_test, parallel = FALSE)
  res_par <- logStirling2(n_test, parallel = TRUE)
  
  expect_identical(res_seq, res_par)
  
  future::plan(future::sequential) # Reset
})
