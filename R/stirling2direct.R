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
#' For extreme cases where \eqn{k < 4} or \eqn{n - k < 5}, the function utilizes
#' fast closed-form polynomial expressions. Additionally, for large \eqn{k} where 
#' the direct summation is slow, the function dynamically branches to a high-performance 
#' C++ algorithm utilizing second-order Eulerian numbers.
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
  if (k < 4 || n - k < 5) {
    if (k == n || k == 1) return(as.bigz(1))
    if (k == n - 1) return(chooseZ(n, 2))
    if (k == 2) return(as.bigz(2)^(n - 1) - 1)
    if (k == 3) return((as.bigz(3)^(n - 1) + 1) %/% 2 - as.bigz(2)^(n - 1))
    if (k == n - 2) return(chooseZ(n, 3) + 3 * chooseZ(n, 4))
    if (k == n - 3) return(chooseZ(n, 4) + 10 * chooseZ(n, 5) + 15 * chooseZ(n, 6))
    if (k == n - 4) return(chooseZ(n, 5) + 25 * chooseZ(n, 6) + 105 * chooseZ(n, 7) + 105 * chooseZ(n, 8))
  }
  
  if (k > n - 2.5 * n^0.74) {
    m <- n - k
    e <- Eulerian2All(m)
    j <- 0:(m - 1)
    return(sum(e * chooseZ(n + m - 1 - j, 2 * m)))
  }
  
  n <- as.bigz(n)
  j <- seq(k - 1, 0, -2)
  abs(as.bigz(sum(chooseZ(k, j)*(j*seq(2, k, 2)^(n - 1) - seq(1, k, 2)^n))/
                factorialZ(k)))
}
