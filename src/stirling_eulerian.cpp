#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdint>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// Blazing fast custom hex encoder for exactly 8 characters (32 bits)
inline void write_hex_8(char* dest, uint32_t val) {
  const char* hex_digits = "0123456789abcdef";
  for (int i = 7; i >= 0; --i) {
    dest[i] = hex_digits[val & 0xF];
    val >>= 4;
  }
}

//' Second-Order Eulerian Numbers
//'
//' Calculates the exact values of the second-order Eulerian numbers 
//' \eqn{\langle\langle n, k \rangle\rangle} for a fixed row \eqn{n} and all 
//' \eqn{k \in \{0, \dots, n-1\}}.
//'
//' @param n A positive integer specifying the row index.
//' @return A \code{bigz} vector of length \code{n} containing the sequence 
//' \eqn{\langle\langle n, 0 \rangle\rangle, \langle\langle n, 1 \rangle\rangle, \dots, \langle\langle n, n-1 \rangle\rangle}.
//' @export
// [[Rcpp::export]]
RObject Eulerian2All(const int n) {
  if (n <= 0) stop("n must be a positive integer");
  
  // 1. Calculate max limbs using the exact log-gamma of the row sum (2n-1)!!
  double ln2 = std::log(2.0);
  double ln_sum = std::lgamma(2.0 * n + 1.0) - ((double)n * ln2) - std::lgamma((double)n + 1.0);
  int d = (int)(ln_sum / (ln2 * 32.0)) + 4; // Convert to base-2 limbs and add safety buffer
  
  // 2. Flat Contiguous Memory Allocation (No vector-of-vector indirection)
  std::vector<uint32_t> e_flat((size_t)n * d, 0);
  uint32_t* e = e_flat.data(); // Grab raw pointer to data block
  
  std::vector<uint32_t> i(n, 1);
  e[0] = 1; // e[0 * d + 0] = 1
  
  // 3. Main Loop via Raw Pointer Arithmetic
  for (int row_idx = 1; row_idx < n; ++row_idx) {
    for (int k = row_idx; k > 0; --k) {
      uint32_t carry = 0;
      uint32_t ii = std::max(i[k], i[k - 1]);
      
      // Pre-calculate exact memory addresses for row pointers
      uint32_t* ek = e + (k * d);
      uint32_t* ekminus1 = e + ((k - 1) * d);
      
      // Hoist 64-bit loop invariants out of the inner loop
      uint64_t c1 = (uint64_t)(k + 1);
      uint64_t c2 = (uint64_t)(2 * row_idx - k + 1);
      
      // Contiguous cache-friendly inner loop
      for (uint32_t j = 0; j < ii; ++j) {
        uint64_t prod = (uint64_t)ek[j] * c1 + (uint64_t)ekminus1[j] * c2 + carry;
        ek[j] = (uint32_t)prod; 
        carry = (uint32_t)(prod >> 32); 
      }
      
      if (carry != 0) {
        ek[ii] = carry;
        i[k] = ii + 1;
      } else {
        i[k] = ii;
      }
    }
  }
  
  // 4. Naked High-Speed Hex Serialization
  CharacterVector str_res(n);
  const char* hex_digits = "0123456789abcdef";
  
  for (int row = 0; row < n; ++row) {
    uint32_t num_limbs = i[row];
    uint32_t* ek = e + (row * d);
    
    if (num_limbs == 1 && ek[0] == 0) {
      str_res[row] = "0x0";
      continue;
    }
    
    // Calculate exact string size needed to allocate memory instantly
    uint32_t top_limb = ek[num_limbs - 1];
    int top_limb_chars = 0;
    uint32_t temp = top_limb;
    while (temp > 0) { top_limb_chars++; temp >>= 4; }
    if (top_limb_chars == 0) top_limb_chars = 1;
    
    int total_chars = 2 + top_limb_chars + 8 * (num_limbs - 1);
    std::string s(total_chars, '0');
    s[0] = '0'; s[1] = 'x';
    
    // Write top limb manually (avoids leading zeros)
    temp = top_limb;
    for (int idx = 1 + top_limb_chars; idx >= 2; --idx) {
      s[idx] = hex_digits[temp & 0xF];
      temp >>= 4;
    }
    
    // Stream all remaining lower limbs using our ultra-fast 8-char block writer
    char* dest = &s[2 + top_limb_chars];
    for (int j = (int)num_limbs - 2; j >= 0; --j) {
      write_hex_8(dest, ek[j]);
      dest += 8;
    }
    
    str_res[row] = s;
  }
  
  Environment gmp_env = Environment::namespace_env("gmp");
  Function as_bigz = gmp_env["as.bigz"];
  
  return as_bigz(str_res);
}

//' Stirling Numbers of the Second Kind
//'
//' Calculates the exact values of the Stirling numbers of the second kind 
//' \eqn{S(n, k)} for a fixed \eqn{n} and all \eqn{k \in \{1, \dots, n\}}.
//'
//' @param n A positive integer specifying the set size.
//' @return A \code{bigz} vector of length \code{n} containing the Stirling 
//' numbers of the second kind \eqn{S(n, 1), S(n, 2), \dots, S(n, n)}.
//' @export
// [[Rcpp::export]]
RObject Stirling2All(const int n) {
  if (n <= 0) stop("n must be a positive integer");
  
  // 1. Calculate max limbs using Berend-Tassa explicit bound for Bell numbers
  double bound_base = 0.792 * (double)n / std::log((double)n + 1.0);
  int d = (int)((double)n * std::log2(bound_base) / 32.0) + 4;
  
  // 2. Flat Contiguous Memory Allocation
  std::vector<uint32_t> e_flat((size_t)n * d, 0);
  uint32_t* e = e_flat.data();
  
  std::vector<uint32_t> i(n, 1);
  e[0] = 1; // S(1, 1) = 1
  
  // 3. Main Loop via Raw Pointer Arithmetic
  for (int row_idx = 1; row_idx < n; ++row_idx) {
    for (int k = row_idx; k > 0; --k) {
      uint32_t carry = 0;
      uint32_t ii = std::max(i[k], i[k - 1]);
      
      uint32_t* ek = e + (k * d);
      uint32_t* ekminus1 = e + ((k - 1) * d);
      
      // Hoist loop invariant multiplier c1 = K = k + 1
      uint64_t c1 = (uint64_t)(k + 1);
      
      // Contiguous cache-friendly inner loop
      for (uint32_t j = 0; j < ii; ++j) {
        uint64_t prod = (uint64_t)ek[j] * c1 + (uint64_t)ekminus1[j] + carry;
        ek[j] = (uint32_t)prod; 
        carry = (uint32_t)(prod >> 32); 
      }
      
      if (carry != 0) {
        ek[ii] = carry;
        i[k] = ii + 1;
      } else {
        i[k] = ii;
      }
    }
  }
  
  // 4. Naked High-Speed Hex Serialization
  CharacterVector str_res(n);
  const char* hex_digits = "0123456789abcdef";
  
  for (int row = 0; row < n; ++row) {
    uint32_t num_limbs = i[row];
    uint32_t* ek = e + (row * d);
    
    if (num_limbs == 1 && ek[0] == 0) {
      str_res[row] = "0x0";
      continue;
    }
    
    // Calculate exact string size needed
    uint32_t top_limb = ek[num_limbs - 1];
    int top_limb_chars = 0;
    uint32_t temp = top_limb;
    while (temp > 0) { top_limb_chars++; temp >>= 4; }
    if (top_limb_chars == 0) top_limb_chars = 1;
    
    int total_chars = 2 + top_limb_chars + 8 * (num_limbs - 1);
    std::string s(total_chars, '0');
    s[0] = '0'; s[1] = 'x';
    
    // Write top limb manually
    temp = top_limb;
    for (int idx = 1 + top_limb_chars; idx >= 2; --idx) {
      s[idx] = hex_digits[temp & 0xF];
      temp >>= 4;
    }
    
    // Stream remaining lower limbs using write_hex_8
    char* dest = &s[2 + top_limb_chars];
    for (int j = (int)num_limbs - 2; j >= 0; --j) {
      write_hex_8(dest, ek[j]);
      dest += 8;
    }
    
    str_res[row] = s;
  }
  
  Environment gmp_env = Environment::namespace_env("gmp");
  Function as_bigz = gmp_env["as.bigz"];
  
  return as_bigz(str_res);
}
