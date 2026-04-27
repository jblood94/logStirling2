#' @importFrom gmp as.bigz chooseZ factorialZ
NULL

#' Stirling Numbers of the Second Kind (Exact)
#'
#' Calculates the exact value of \eqn{S(n,k)} using \code{bigz} integers.
#'
#' @param n Positive integer set size.
#' @param k Integer subset size in \code{1:n}.
#'
#' @details Implements the explicit formula for positive arguments:
#'
#'   \deqn{S(n,k)=\frac{1}{k!}\sum_{j=1}^k(-1)^{k-j}\binom{k}{j}j^n}
#'   \deqn{=\frac{1}{k!}\sum_{j=1}^k\binom{-(j+1)}{k-j}j^n}
#'
#' This is a "direct" calculation similar to \code{gmp::Stirling2(method =
#' "direct")}, but without cancellation errors for "large" n.
#'
#' @return A \code{bigz} object.
#'
#' @seealso \code{\link{logStirling2}} for log-scale calculations accepting
#' vectors for \code{n} and \code{k}.
#'
#' @examples
#' # Basic usage
#' stirling2direct(5, 3)
#'
#' # Comparison with the log version
#' mapply(\(k) log(stirling2direct(200, k)), 10:20)
#' logStirling2(200, 10:20)
#'
#' @export
stirling2direct <- function(n, k) {
  n <- as.bigz(n)
  j <- seq(k - 1, 0, -2)
  abs(as.bigz(sum(chooseZ(k, j)*(j*seq(2, k, 2)^(n - 1) - seq(1, k, 2)^n))/
                factorialZ(k)))
}
