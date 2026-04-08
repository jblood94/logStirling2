#' @useDynLib logStirling2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

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
#'
#' @details The function dispatches to one of three C++ routines (\code{Row_C},
#'   \code{All_C}, or \code{Mult_C}) depending on the sparsity of the input
#'   vector \code{n}. If \code{length(n) == 1 && length(k) == 1} is \code{TRUE},
#'   the function calls \code{stirling2direct}. If \eqn{n \ge 1000}, the
#'   function automatically loads and utilizes the pre-computed long-double
#'   state blocks stored in \code{sysdata.rda} to accelerate calculations and
#'   maintain precision.
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
#' # 1. Single value calculation (uses stirling2direct)
#' logStirling2(10, 5)
#'
#' # 2. Matrix output for a range of n and specific k
#' # Useful for comparing growth across set sizes
#' logStirling2(n = 5:10, k = 1:5, as.matrix = TRUE)
#'
#' # 3. Vector output with 'ones' filtered
#' # This returns only the "non-trivial" values (1 < k < n)
#' logStirling2(n = 5, k = NULL, as.matrix = FALSE, ones = FALSE)
#'
#' # 4. Using the Temme approximation for large n
#' # Fast even with very large n
#' logStirling2Temme(n = 1e5)
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
  states <- if (nu[length(nu)] > 999) c(list(NULL), logS_states) else list(NULL)
  bu <- pmin(as.integer(names(blocks)) + 1, length(states))
  sb <- states[bu]
  ni <- match(n, c(nu0, nu))

  # Process blocks
  S <- unlist(lapply(1:length(blocks), function(i) {
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

#' @rdname logStirling2
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
