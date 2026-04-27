#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
RawVector logStirling2State_C(int n, Nullable<RawVector> state = R_NilValue) {
  // Returns as a raw vector elements k = 1 to k = n - 1 of the nth row of the
  // log-Stirling numbers of the second kind. Accepts a length-16 raw vector
  // containing the values of a row to start the calculations.

  // 1. Determine precision-safe vector size
  // long double is usually 12 or 16 bytes depending on OS
  size_t ld_size = sizeof(long double);
  // long double intermediate calculations
  std::vector<long double> s(n - 1, 0.0);
  int n0;
  // 2. Start from state if provided
  if (state.isNotNull()) {
    RawVector rv(state);
    n0 = rv.size()/ld_size;
    std::memcpy(s.data(), rv.begin(), rv.size());
  } else n0 = 1;

  long double* p_s = s.data();

  // 3. Core computation loop
  for (int j = n0; j < n - 1; j++) {
    for (int k = j; k > 0; k--) {
      p_s[k] += log1pl(k + expl(p_s[k - 1] - p_s[k]));
    }
  }

  RawVector state_out(s.size()*ld_size);
  std::memcpy(state_out.begin(), s.data(), s.size()*ld_size);
  return state_out;
}

// [[Rcpp::export]]
NumericVector logStirling2Row_C(int n, Nullable<RawVector> state = R_NilValue) {
  // Returns elements k = 2 to k = n - 1 of the nth row of the log-Stirling
  // numbers of the second kind. Accepts a length-16 raw vector containing the
  // values of a row to start the calculations.
  if (n < 3) return NumericVector(0);
  // long double intermediate calculations
  std::vector<long double> s(n - 1, 0.0);
  int n0;
  // 2. Start from state if provided
  if (state.isNotNull()) {
    RawVector rv(state);
    n0 = rv.size()/16;
    std::memcpy(s.data(), rv.begin(), rv.size());
  } else n0 = 1;

  long double* p_s = s.data();

  // 3. Core computation loop
  for (int j = n0; j < n - 1; j++) {
    for (int k = j; k > 0; k--) {
      p_s[k] += log1pl(k + expl(p_s[k - 1] - p_s[k]));
    }
  }

  return NumericVector(std::next(s.begin()), s.end());
}

// [[Rcpp::export]]
NumericVector logStirling2All_C(int n, Nullable<RawVector> state = R_NilValue) {
  // Returns a slice of rows from the log Stirling 2 triangle. Initialized with
  // the previous row provided as a length-16 raw vector (or starting from n =
  // 3, if omitted).
  if (n < 3) return NumericVector(0);
  // long double intermediate calculations
  std::vector<long double> s(n - 1, 0.0);
  int n0;
  // Resume from state if provided
  if (state.isNotNull()) {
    RawVector rv(state);
    n0 = rv.size()/16 + 1;
    std::memcpy(s.data(), rv.begin(), rv.size());
  } else n0 = 2;

  // Initialize the results table
  NumericVector sall(n/2.0*(n - 3) - n0/2.0*(n0 - 3) + n0 - 2);
  long double* p_s = s.data();
  double* p_sall = sall.begin();
  long long i = n0 - 3;

  // Pre-fill sall with the initial state row
  // (If state is NULL, n0 = 2, so this loop is naturally skipped)
  for (int k = 1; k <= n0 - 2; k++) p_sall[k - 1] = p_s[k];

  // Core computation loop
  for (int j = n0 - 1; j < n - 1; j++) {
    for (int k = j; k > 0; k--) {
      p_s[k] += log1pl(k + expl(p_s[k - 1] - p_s[k]));
      p_sall[i + k] = p_s[k];
    }
    i += j;
  }

  return sall;
}

// [[Rcpp::export]]
NumericVector logStirling2Mult_C(IntegerVector n,
                                 Nullable<RawVector> state = R_NilValue) {
  // Returns a slice of rows from the log Stirling 2 triangle. Initialized with
  // a stored row provided as a length-16 raw vector (or starting from n = 3, if
  // omitted). `n` is assumed to be sorted ascending.

  int n_len = n.size();
  // long double intermediate calculations
  std::vector<long double> s(n[n_len - 1] - 1, 0.0);
  int n0;
  // 2. Resume from state if provided
  if (state.isNotNull()) {
    RawVector rv(state);
    n0 = rv.size()/16 + 1;
    std::memcpy(s.data(), rv.begin(), rv.size());
  } else n0 = 2;

  // Initialize the results table
  NumericVector sall(sum(n) - 2*n_len);
  long double* p_s = s.data();
  double* p_sall = sall.begin();
  long long i = 0;

  // 3. Core computation loop
  for (int idx = 0; idx < n_len; idx++) {
    int nextn = n[idx];
    for (int j = n0 - 1; j < nextn - 1; j++) {
      for (int k = j; k > 0; k--) {
        p_s[k] += log1pl(k + expl(p_s[k - 1] - p_s[k]));
      }
    }

    if (nextn > n0) n0 = nextn;
    // i++;
    for (int k = 1; k < n0 - 1; k++) p_sall[i++] = p_s[k];
  }


  return sall;
}

// [[Rcpp::export]]
RawVector create_ld_state_C(CharacterVector x) {
  size_t n = x.size();
  size_t ld_size = sizeof(long double);

  // Create a buffer for long doubles
  std::vector<long double> buffer(n);

  for(size_t i = 0; i < n; ++i) {
    // Convert R string to C++ long double
    buffer[i] = std::strtold(x[i], nullptr);
  }

  // Serialize to RawVector
  RawVector rv(n * ld_size);
  std::memcpy(rv.begin(), buffer.data(), n * ld_size);

  return rv;
}
