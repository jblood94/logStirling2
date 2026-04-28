#' @useDynLib logStirling2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# Internal environment to hold cached states across function calls
.state_env <- new.env(parent = emptyenv())

#' Logarithms of Stirling Numbers of the Second Kind
#'
#' Calculates the natural logarithm of Stirling numbers of the second kind,
#' \eqn{S(n, k)}, which represent the number of ways to partition a set of
#' \eqn{n} elements into \eqn{k} non-empty subsets.
#'
#' @param n Integer vector of set sizes. Coerced to natural numbers (floor).
#' @param k Integer vector of subset sizes. Coerced to natural numbers (floor).
#'   If \code{NULL}, returns all available \eqn{k} for each \eqn{n}.
#' @param as.matrix Logical; if \code{TRUE}, returns a matrix where rows
#'   correspond to \code{n} and columns to \code{k}. If \code{FALSE}, returns a
#'   flat vector.
#' @param ones Logical; if \code{FALSE}, excludes the trivial cases where \eqn{k
#'   = 1} and \eqn{k = n} (where \eqn{S(n, k) = 1}). This is automatically set
#'   to \code{TRUE} if \code{as.matrix} is \code{TRUE}, \code{k} is explicitly
#'   provided, or if \code{any(n < 3)} is \code{TRUE}.
#' @param twoterms Logical; if \code{TRUE}, uses Temme's two-term approximation.
#'   If \code{FALSE}, uses the one-term approximation.
#'
#' @details The function dispatches to one of three C++ routines (\code{Row_C},
#'   \code{All_C}, or \code{Mult_C}) depending on the sparsity of the input
#'   vector \code{n}.
#'
#'   For systems supporting 16-byte \code{long double} precision, if \eqn{n \ge
#'   1000}, the function automatically searches for pre-computed state blocks.
#'   If found in the package namespace or the user's data directory
#'   (\code{tools::R_user_dir}), these blocks are used to dramatically
#'   accelerate calculations. If missing, the full table is computed on-the-fly.
#'   If unsupported (e.g., Apple Silicon/ARM64), the full table is computed
#'   using standard double precision.
#'
#'   \code{logStirling2Temme} provides a high-speed asymptotic approximation
#'   based on Temme's method, which is functionally identical in interface but
#'   trades exactness for performance at very large \eqn{n}.
#'
#' @return A numeric matrix or vector containing \eqn{\ln(S(n, k))}. For \eqn{k
#'   > n}, values are returned as \code{NA_real_}.
#'
#' @references Temme, N. M. (1993). Asymptotic estimates of Stirling numbers.
#'   \emph{Studies in Applied Mathematics}, 89(3), 233-243.
#'
#' @examples
#' # 1. Matrix output for specified n and k
#' logStirling2(n = 5:8, k = 2:5, as.matrix = TRUE)
#'
#' # 2. Vector output with 'ones' filtered
#' # This returns only the "non-trivial" values (1 < k < n)
#' logStirling2(n = 8:10, k = NULL, as.matrix = FALSE, ones = FALSE)
#'
#' # 3. Full row with large n
#' s <- logStirling2(n = 1e3, as.matrix = FALSE)
#' length(s)
#' s[10:13]
#'
#' # 4. Temme's asymptotic approximation — fast even for very large n
#' s <- logStirling2Temme(n = 1e5, as.matrix = FALSE)
#' s[1000:1003]
#'
#' @name logStirling2
NULL

#' @rdname logStirling2
#' @export
logStirling2 <- function(n, k = NULL, as.matrix = TRUE, ones = TRUE) {
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

  states <- if (.Machine$sizeof.longdouble == 16) {
    ns <- asNamespace("logStirling2")

    if (!is.null(.state_env$logS_states)) {
      .state_env$logS_states
    } else if (exists("logS_states", envir = ns, inherits = FALSE)) {
      get("logS_states", envir = ns)
    } else {
      cache_file <- file.path(tools::R_user_dir("logStirling2", "data"),
                              "logStirling2_cache_v01.rds")
      if (file.exists(cache_file)) {
        .state_env$logS_states <- readRDS(cache_file)
      } else {
        warning(paste("Cached data not found. The full Stirling table will be",
                      "calculated. Run `logStirling2::get_state_data()`",
                      "to download cached data from the Git repository."),
                call. = FALSE)
        NULL
      }
    }
  } else {
    warning(paste("16-byte extended precision (long double) is not supported",
                  "by this system. The full Stirling table will be calculated",
                  "with double precision."),
            call. = FALSE)
    NULL
  }

  if (is.null(states)) {
    blocks <- list(nu)
    bu <- list(3L)
  } else {
    state_lens <- c(3, lengths(states)/16 + 1)
    blocks <- split(nu, findInterval(nu, state_lens))
    bu <- as.numeric(names(blocks))
    states <- states[bu - 1]
    if (bu[1] == 1) states <- c(list(NULL), states)
  }

  # Process blocks
  S <- unlist(lapply(1:length(blocks), function(i) {
    if (length(blocks[[i]]) == 1) { # single row
      logStirling2Row_C(blocks[[i]], state = states[[i]])
    } else if (blocks[[i]][1] == state_lens[bu[i]] &&
               all(diff(blocks[[i]]) == 1)) { # all rows
      logStirling2All_C(blocks[[i]][length(blocks[[i]])], state = states[[i]])
    } else logStirling2Mult_C(blocks[[i]], state = states[[i]])
  }), FALSE, FALSE)

  S <- unname(split(S, rep.int(nu, nu - 2)))
  ni <- match(n, c(nu0, nu))

  if (is.null(k)) {
    if (!as.matrix && !ones) return(unlist(S[ni])) else k <- 1:nu[length(nu)]
  }

  S <- c(lapply(nu0, \(i) numeric(i)[k]), lapply(S, \(s) c(0, s, 0)[k]))[ni]

  if (as.matrix) {
    `dimnames<-`(do.call(rbind, S), list(c(nu0, nu)[ni], k))
  } else unlist(lapply(S, \(s) s[!is.na(s)]), FALSE, FALSE)
}

#' @rdname logStirling2
#' @export
logStirling2Temme <- function(n, k = NULL, as.matrix = TRUE, ones = TRUE,
                              twoterms = TRUE) {
  # Returns log(S(n, k)) using Temme's approximation for
  # k = 4:(n - 3) and exact formulas for k = c(2, 3, n - 2, n - 1).
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
    s[i] <- if (twoterms) {
      logStirling2Temme2(n[i], k[i])
    } else {
      # 1-term asymptotic expansion
      n <- n[i]
      k <- k[i]
      v <- n/k
      G <- -lambertW0(-v*exp(-v))
      G1 <- 1 - G
      0.5*((logv1 <- log(v - 1)) - log(v*G1)) + lchoose(n, k) +
        (n - k)*(logv1 - log(v - G)) + n*log(k) + k*(G1 - log(n))
    }
  }

  if (as.matrix) {
    dim(s) <- c(nn, nk)
    dimnames(s) <- list(n0, k0)
  }

  s
}

logStirling2Temme2 <- function(n, k) {
  # Returns log(S(n, k)) using Temme's 2-term asymptotic approximation
  x <- n/k
  t <- x - 1
  i <- which(x < 37.43)

  if (length(i)) {
    v <- x[i]
    y <- lambertW0(-v*exp(-v)) + v
    j <- which(abs(y/expm1(-y) + v) > 5e-15)

    if (length(j)) {
      z <- y[j]
      z <- z - (z + v[j]*expm1(-z))/(1 - z/expm1(z))
      y[j] <- z
    }

    x[i] <- y
  }

  ex <- exp(-x)
  ey <- -1/expm1(-x)
  x0 <- k*t
  x4 <- ey^2
  x8 <- 1/x
  x1 <- x8^2
  x2 <- k*ex*x4
  x3 <- 1/(n*x1 - x2)
  x5 <- k/t
  x6 <- x5*x3
  x7 <- sqrt(x6)
  x9 <- t*x8
  x11 <- x6*x7
  x12 <- 1/t^2
  x13 <- n*x8*x1
  x14 <- k*(ex + 1)*ex*ey*x4
  x15 <- (x11*(x13 - 0.5*x14) - k*x12)/3
  x16 <- x15*x3
  x19 <- t*x3
  lchoose(n, k) + k*(x - log(ey)) - n*log(x) - x0 + (n - k)*log(x0) +
    log(
      x7*x9 - x9*(
        t*x1*x11 + 2*x16/x7 - 3*x16*x9 - x6*x8 - 3*x19*(
          (k*x3)^2*x12*(n/(4*x^4) - k*ex*x4^2*(ex*(ex + 4) + 1)/24) -
            k/(4*t^3) + 0.5*x15^2*x19*(n*x1 - x2)/k +
            x15*x5*x3^2*(0.5*x14 - x13)/x7
        )/x7
      )/k
    )
}

lambertW0 <- function(x) {
  # Principal branch of Lambert's W function for `x > 1`
  w <- (w <- log1p(x))*(1 - log1p(w)/(w + 2))

  for (. in 1:2) {
    ew <- exp(w)
    f <- w*ew - x
    w1 <- w + 1
    w <- w - f/(ew*w1 - (w1 + 1)*f/(2*w1))
  }

  w
}
