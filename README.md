# logStirling2

An R package for computing logarithms of Stirling numbers of the second kind,
S(n, k), accurately and efficiently for arbitrarily large n.

## Motivation

The Stirling number of the second kind S(n, k) counts the number of ways to
partition a set of n distinct objects into k non-empty subsets. These numbers
grow extremely rapidly and overflow 64-bit floating point for moderate n,
making direct computation impractical. `logStirling2` sidesteps this by
working on the log scale via a numerically stable recurrence in C++, enabling
fast and accurate computation for full rows for n past 50,000 — well past the
practical limit of existing packages.

| Approach | Handles large n? | Vectorized? | Precision | Speed |
|---|:---:|:---:|:---:|:---:|
| `copula::Stirling2` | no (overflows n ≥ 220) | no | double (direct formula) | slow (stores values) |
| `gmp::Stirling2` | yes | no | exact (big.z) | very slow |
| `occupancy::logStirling` | yes | yes | double (row recursion) | slow |
| `logStirling2` | yes | yes | long double + checkpoints (row recursion) | fast |
| `logStirling2Temme` | yes | yes | asymptotic | fastest |

## Installation

The package is not yet on CRAN. Install the development version from GitHub
with:

```r
# install.packages("remotes")
remotes::install_github("jblood94/logStirling2")
```

A C++ compiler (Rtools on Windows, Xcode Command Line Tools on macOS) is
required to build the package from source.

## Functions

The package exports three user-facing functions:

**`logStirling2(n, k = NULL, as.matrix = TRUE, ones = TRUE)`** — the main
workhorse. Accepts vectors for both `n` and `k` and returns results as a
matrix (rows = n, columns = k) or a flat vector. When `k = NULL`, all valid
k for each n are returned. Setting `ones = FALSE` drops the trivial edge cases
S(n, 1) = S(n, n) = 1 (i.e., log = 0) from the output. For a single (n, k)
pair, it delegates to `stirling2direct` for an exact result prior to taking the
natural logarithm.

**`logStirling2Temme(n, k = NULL, as.matrix = TRUE, ones = TRUE)`** — same
interface as `logStirling2`, but uses Temme's (1993) asymptotic approximation
for k = 4:(n−2), with exact closed-form expressions at the edges (k = 2, 3,
n−2, n−1). Much faster than the recurrence for very large n, at the cost of
exactness.

**`stirling2direct(n, k)`** — computes the exact value of S(n, k) as a `bigz`
integer via the explicit formula using arbitrary-precision arithmetic from
`gmp`. Called automatically by `logStirling2` for scalar (n, k) inputs.

## Usage

```r
library(logStirling2)

# Single value: log S(10, 5)
logStirling2(10, 5)

# A slice of the Stirling triangle as a matrix
logStirling2(n = 5:8, k = 2:5, as.matrix = TRUE)

# All non-trivial k for a single n, as a vector
logStirling2(n = 10, k = NULL, as.matrix = FALSE, ones = FALSE)

# Full row for large n (no overflow)
logStirling2(n = 5000, k = NULL, as.matrix = FALSE, ones = FALSE)

# Temme approximation — fast even for very large n
logStirling2Temme(n = 1e5)

# Exact big-integer result for a scalar pair
stirling2direct(200, 10)
```

## Algorithm

`logStirling2` applies the standard Stirling recurrence on the log scale.
Writing s(n, k) = log S(n, k), the recurrence S(n+1, k) = k·S(n, k) + S(n, k−1)
becomes:

$$s(n+1,\, k) = s(n,\, k) + \log\!\bigl(k + e^{\,s(n,\, k-1) - s(n,\, k)}\bigr)$$

The C++ implementation evaluates this with `long double` precision using
`log1pl` and `expl` for numerical stability, sweeping k from high to low
within each row so the update is in-place.

For n ≥ 1,000, the function automatically loads pre-computed starting states
stored in `sysdata.rda` — 50 rows of the triangle (every 1,000th row from
n = 1,000 to n = 50,000), computed externally with FLINT and ARB and stored as
raw vectors at `long double` precision. These checkpoints both accelerate
computation and preserve precision for large n by minimizing the number of
recurrence steps needed.

Depending on the structure of the input vector `n`, the C++ layer dispatches
to one of three routines:

- **`logStirling2Row_C`** — a single row of the triangle
- **`logStirling2All_C`** — a contiguous range of rows
- **`logStirling2Mult_C`** — an arbitrary sorted set of rows

`logStirling2Temme` uses the approximation of Temme (1993) via the Lambert W
function (supplied by the `lamW` package).

## Dependencies

- [`Rcpp`](https://CRAN.R-project.org/package=Rcpp) — C++ integration
- [`gmp`](https://CRAN.R-project.org/package=gmp) — arbitrary-precision
  arithmetic for `stirling2direct`
- [`lamW`](https://CRAN.R-project.org/package=lamW) — Lambert W function for
  `logStirling2Temme`

## References

Temme, N. M. (1993). Asymptotic estimates of Stirling numbers. *Studies in
Applied Mathematics*, 89(3), 233–243.

## License

GPL (≥ 3)
