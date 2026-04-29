# Submission comments

This is a new release of `logStirling2`.

## Test environments
* local Windows installation, R-release [cite: 1.2]
* ubuntu-latest (on GitHub Actions), R-release [cite: 1.4, 2.2]
* macos-latest (on GitHub Actions), R-release [cite: 1.4, 2.2]
* macos-arm64 (R-hub), R-release [cite: 1.4, 2.3]
* Fedora Linux, R-devel, clang-ASAN [cite: 1.4, 2.3]
* Ubuntu Linux, R-devel, nold (No Long Double) [cite: 1.4, 2.3]

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. [cite: 1.2]

## Major Changes
* Implemented a tiered data strategy (RAM, local disk, and C++ fallback) to handle high-precision calculations of log Stirling numbers of the second kind without exceeding CRAN package size limits. [cite: 1.2, 2.1]
* Moved large internal data objects (`sysdata.rda`) to `data-raw` to keep the source package lightweight. [cite: 1.2, 3.2]
* Added platform-specific logic to detect `long double` (16-byte) support, providing high-precision results where hardware allows and falling back gracefully to standard double precision on other architectures (e.g., ARM64). [cite: 1.2, 2.3]
* Exported `logStirling2Temme()` as a standalone function for users requiring fast asymptotic approximations for very large n, distinct from the exact solver in `logStirling2()`. [cite: 2.1]
* Removed the `lamW` dependency, replaced with internal C++ implementations. [cite: 1.2, 2.1]
