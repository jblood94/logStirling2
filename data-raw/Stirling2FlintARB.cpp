#include <Rcpp.h>
#include <cmath>
#include <mpfr.h>      
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/arb_poly.h>
#include <flint/arf.h>

using namespace Rcpp;

// [[Rcpp::export]]
RawVector logstirling2_raw_flintARB(int n, slong prec = 256) {
  if (n < 2) return RawVector(0);
  
  // Peak intermediate term ~ e^n requires ~1.44n bits to represent;
  // factor of 2 gives comfortable headroom above that floor.
  slong wp = std::max(prec, (slong)(2 * n + 64));
  
  arb_poly_t A, B, C;
  arb_poly_init(A); arb_poly_init(B); arb_poly_init(C);
  arb_poly_fit_length(A, n + 1);
  arb_poly_fit_length(B, n + 1);
  
  arb_t val, fac, j_pow_n, inv_fac;
  arb_init(val); arb_init(fac);
  arb_init(j_pow_n); arb_init(inv_fac);
  
  arb_zero(val);
  arb_poly_set_coeff_arb(A, 0, val);
  arb_one(val);
  arb_poly_set_coeff_arb(B, 0, val);
  
  arb_one(fac);
  for (int j = 1; j <= n; j++) {
    arb_mul_ui(fac, fac, j, wp);          // wp throughout
    arb_inv(inv_fac, fac, wp);
    
    arb_set_ui(j_pow_n, j);
    arb_pow_ui(j_pow_n, j_pow_n, n, wp);
    arb_mul(val, j_pow_n, inv_fac, wp);
    arb_poly_set_coeff_arb(A, j, val);
    
    if (j % 2 != 0) arb_neg(inv_fac, inv_fac);
    arb_poly_set_coeff_arb(B, j, inv_fac);
  }
  
  arb_poly_mullow(C, A, B, n, wp);        // wp here too
  
  RawVector result((n - 1) * 16);
  long double* res_ptr = reinterpret_cast<long double*>(RAW(result));
  
  mpfr_t m;
  mpfr_init2(m, 80);
  
  for (int k = 1; k < n; k++) {
    arb_poly_get_coeff_arb(val, C, k);
    if (arb_is_positive(val)) {
      arb_log(val, val, 80);
      arf_get_mpfr(m, arb_midref(val), MPFR_RNDN);
      res_ptr[k - 1] = mpfr_get_ld(m, MPFR_RNDN);
    } else {
      res_ptr[k - 1] = -INFINITY;
    }
  }
  
  mpfr_clear(m);
  arb_poly_clear(A); arb_poly_clear(B); arb_poly_clear(C);
  arb_clear(val); arb_clear(fac); arb_clear(j_pow_n); arb_clear(inv_fac);
  
  return result;
}