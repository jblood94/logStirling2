#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector logStirling2Row(int n) {
  // Returns elements k = 1 to k = n - 1 of the nth row of the log-Stirling
  // numbers of the second kind.
  
  std::vector<long double> s(n - 1, 0.0);
  long double* p_s = s.data();

  for (int j = 1; j < n - 1; ++j) {
    for (int k = j; k > 0; --k) {
      p_s[k] += log1pl(k + expl(p_s[k - 1] - p_s[k]));
    }
  }
  
  Rcpp::NumericVector out(s.begin(), s.end());
  return out;
}

// [[Rcpp::export]]
NumericVector logStirling2All(int n) {
  // Returns the first `n` rows of the log-Stirling numbers of the second kind.
  // Elements are store row-wise from k = 1 to k = n - 1.
  
  long long n_ll = n;
  NumericVector sall(n_ll*(n_ll - 3)/2 + 1);
  std::vector<long double> s(n - 1, 0.0);
  double* p_sall = sall.begin();
  long double* p_s = s.data();
  long long i = -1;
  
  for (int j = 1; j < n - 1; ++j) {
    for (int k = j; k > 0; --k) {
      long double sk = p_s[k];
      long double sk1 = p_s[k - 1];
      long double updated = sk + log1pl(k + expl(sk1 - sk));
      p_s[k] = updated;
      p_sall[i + k] = updated;
    }
    i += j;
  }
  
  return sall;
}

// [[Rcpp::export]]
List logStirling2Row_Step(int n, Nullable<RawVector> state = R_NilValue) {
  // Returns non-zero elements of row n of the table of log-Stirling
  // numbers of the second kind (k = 2 to k = n - 1). Also returns a state
  // vector (raw) containing the values with long double precision.
  
  // 1. Determine precision-safe vector size
  // long double is usually 12 or 16 bytes depending on OS
  size_t ld_size = sizeof(long double);
  // long double intermediate calculations
  std::vector<long double> s(n - 1, 0.0);
  int n0;
  // 2. Resume from state if provided
  if (state.isNotNull()) {
    RawVector rv(state);
    n0 = rv.size()/ld_size + 1;
    std::memcpy(s.data(), rv.begin(), rv.size());
  } else n0 = 2;
  
  // Initialize the results table
  NumericVector sall(n/2.0*(n - 3) - n0/2.0*(n0 - 3));
  long double* p_s = s.data();
  long long i = -1;
  
  // 3. Core computation loop
  for (int j = n0 - 1; j < n - 1; j++) {
    for (int k = j; k > 0; k--) {
      p_s[k] += log1pl(k + expl(p_s[k - 1] - p_s[k]));
    }
    i += j;
  }
  
  // 4. Serialize the final long double state into a RawVector
  RawVector state_out(s.size()*ld_size);
  std::memcpy(state_out.begin(), s.data(), s.size()*ld_size);
  
  return List::create(
    Named("row") = s,
    Named("state") = state_out,
    Named("last_n") = n
  );
}

// [[Rcpp::export]]
RawVector create_ld_state(CharacterVector x) {
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