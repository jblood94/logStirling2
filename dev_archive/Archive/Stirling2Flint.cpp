#include <Rcpp.h>
#include <cmath>

extern "C" {
  #include <flint/flint.h>
  #include <flint/fmpz.h>
  #include <flint/fmpz_poly.h>
}

using namespace Rcpp;

// Notice we are returning a CharacterVector now
// [[Rcpp::export]]
CharacterVector stirling2_row_flint_fft(int n)
{
  fmpz_poly_t A, B, C;
  fmpz_poly_init(A);
  fmpz_poly_init(B);
  fmpz_poly_init(C);
  
  fmpz_t term_A, term_B, n_over_fact;
  fmpz_init(term_A);
  fmpz_init(term_B);
  fmpz_init(n_over_fact);
  
  fmpz_one(n_over_fact); 
  
  for (int i = n; i >= 0; i--) {
    fmpz_set_ui(term_A, i);
    fmpz_pow_ui(term_A, term_A, n);
    fmpz_mul(term_A, term_A, n_over_fact);
    fmpz_poly_set_coeff_fmpz(A, i, term_A);
    fmpz_set(term_B, n_over_fact);
    if (i % 2) fmpz_neg(term_B, term_B);
    fmpz_poly_set_coeff_fmpz(B, i, term_B);
    if (i > 0) fmpz_mul_ui(n_over_fact, n_over_fact, i);
  }
  
  fmpz_poly_mul(C, A, B);
  fmpz_pow_ui(n_over_fact, n_over_fact, 2);
  
  // Create an Rcpp CharacterVector to hold the raw strings
  CharacterVector result(n + 1);
  fmpz_t val;
  fmpz_init(val);
  
  for (int k = 0; k <= n; k++) {
    if (k < C->length) fmpz_set(val, C->coeffs + k); else fmpz_zero(val);
    fmpz_fdiv_q(val, val, n_over_fact);
    char *str = fmpz_get_str(NULL, 10, val);
    result[k] = str; // Implicitly converts to Rcpp::String
    flint_free(str);
  }
  
  // 1. Clear your specific variables
  fmpz_poly_clear(A);
  fmpz_poly_clear(B);
  fmpz_poly_clear(C);
  
  fmpz_clear(term_A);
  fmpz_clear(term_B);
  fmpz_clear(n_over_fact);
  fmpz_clear(val);
  
  // 2. THE FIX: Tell FLINT to dump its global/thread caches
  flint_cleanup();
  
  return result;
}

// [[Rcpp::export]]
CharacterVector stirling2_row_flint_fft2(int n)
{
  fmpz_poly_t A, B, C;
  fmpz_poly_init(A);
  fmpz_poly_init(B);
  fmpz_poly_init(C);
  
  fmpz_t term_A, term_B, n_over_fact;
  fmpz_init(term_A);
  fmpz_init(term_B);
  fmpz_init(n_over_fact);
  
  fmpz_one(n_over_fact); 
  
  for (int i = n; i >= 0; i--) {
    fmpz_set_ui(term_A, i);
    fmpz_pow_ui(term_A, term_A, n);
    fmpz_mul(term_A, term_A, n_over_fact);
    fmpz_poly_set_coeff_fmpz(A, i, term_A);
    fmpz_set(term_B, n_over_fact);
    if (i % 2) fmpz_neg(term_B, term_B);
    fmpz_poly_set_coeff_fmpz(B, i, term_B);
    if (i > 0) fmpz_mul_ui(n_over_fact, n_over_fact, i);
  }
  
  fmpz_poly_mul(C, A, B);
  fmpz_pow_ui(n_over_fact, n_over_fact, 2);
  
  // Create an Rcpp CharacterVector to hold the raw strings
  CharacterVector result(n + 1);
  fmpz_t val;
  fmpz_init(val);
  
  for (int k = 0; k <= n; k++) {
    if (k < C->length) fmpz_set(val, C->coeffs + k); else fmpz_zero(val);
    fmpz_fdiv_q(val, val, n_over_fact);
    char *str = fmpz_get_str(NULL, 16, val);
    std::string s_val = str;
    result[k] = "0x" + s_val;
    flint_free(str);
  }
  
  // 1. Clear your specific variables
  fmpz_poly_clear(A);
  fmpz_poly_clear(B);
  fmpz_poly_clear(C);
  
  fmpz_clear(term_A);
  fmpz_clear(term_B);
  fmpz_clear(n_over_fact);
  fmpz_clear(val);
  
  // 2. THE FIX: Tell FLINT to dump its global/thread caches
  flint_cleanup();
  
  return result;
}

// [[Rcpp::export]]
RawVector stirling2_log_raw_flint(int n) {
  if (n < 2) return RawVector(0); // Boundary check
  
  fmpz_poly_t A, B, C;
  fmpz_poly_init(A);
  fmpz_poly_init(B);
  fmpz_poly_init(C);
  
  fmpz_t term, nf;
  fmpz_init(term);
  fmpz_init(nf);
  fmpz_one(nf);
  
  // Calculate n! to use for the scaling factor and the loop
  for (long i = 1; i <= n; i++) fmpz_mul_ui(nf, nf, i);
  
  // Scaling factor for the logs: 2*log(n!)
  slong e_nf;
  double m_nf = fmpz_get_d_2exp(&e_nf, nf);
  long double log_nf = std::logl((long double)m_nf) +
    (long double)e_nf*std::logl(2.0L);
  long double log_divisor = 2.0L*log_nf;
  
  // Building Polys: This still needs to run 0 to n for the FFT/Product to be correct
  fmpz_t current_f;
  fmpz_init_set_ui(current_f, 1); // This will track n!/i!
  
  for (int i = n; i >= 0; i--) {
    // a_i = i^n * (n! / i!)
    fmpz_set_ui(term, i);
    fmpz_pow_ui(term, term, n);
    fmpz_mul(term, term, current_f);
    fmpz_poly_set_coeff_fmpz(A, i, term);
    
    // b_i = (-1)^i * (n! / i!)
    fmpz_set(term, current_f);
    if (i % 2 != 0) fmpz_neg(term, term);
    fmpz_poly_set_coeff_fmpz(B, i, term);
    
    if (i > 0) fmpz_mul_ui(current_f, current_f, i);
  }
  
  fmpz_poly_mul(C, A, B);
  
  // EXTRACTION: Only k = 1 to n-1
  // Result size is (n-1) long doubles
  RawVector result((n - 1) * sizeof(long double));
  long double* res_ptr = reinterpret_cast<long double*>(RAW(result));
  fmpz_t val;
  fmpz_init(val);
  
  for (int k = 1; k < n; k++) {
    fmpz_poly_get_coeff_fmpz(val, C, k);
    if (fmpz_sgn(val) <= 0) {
      res_ptr[k - 1] = -INFINITY;
    } else {
      slong e_v;
      double m_v = fmpz_get_d_2exp(&e_v, val);
      // Index is k-1 because we skipped k=0
      res_ptr[k - 1] = (std::logl((long double)m_v) +
        (long double)e_v*std::logl(2.0L)) - log_divisor;
    }
  }
  
  // Cleanup
  fmpz_poly_clear(A); fmpz_poly_clear(B); fmpz_poly_clear(C);
  fmpz_clear(term); fmpz_clear(nf); fmpz_clear(val); fmpz_clear(current_f);
  flint_cleanup();
  return result;
}