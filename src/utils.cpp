#include "utils.h"
#include <cmath>
using namespace Rcpp;

// Logarithm of the binomial coefficient C(n, k) via lgamma.
// lgamma(x) returns log(Gamma(x)) = log((x-1)!), so:
// log(C(n, k)) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
double cpp_lchoose(int n, int k) {
  if (k < 0 || k > n) return R_NegInf;
  if (k == 0 || k == n) return 0.0;
  return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}

// This function calculates the binomial coefficient C(n, k),
// which represents the number of ways to choose k objects from a set
// of n objects without regard to order. The function uses an optimized
// approach to calculate the binomial coefficient for faster execution,
// leveraging the symmetry property C(n, k) = C(n, n-k) to reduce the
// number of operations.
double cpp_choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  if (k == 0 || k == n) return 1;

  // Optimized calculation of the binomial coefficient
  k = std::min(k, n - k); // Using symmetry C(n, k) = C(n, n-k)
  double result = 1.0;
  for (int i = 1; i <= k; ++i) {
    result *= (n - i + 1) / static_cast<double>(i);
  }
  return result;
}

static std::vector<int> column_as_int(DataFrame& df, const char* name) {
  IntegerVector col = as<IntegerVector>(df[name]);
  return as<std::vector<int>>(col);
}

BaseCounts read_base_counts(DataFrame& df) {
  BaseCounts out;
  out.A = column_as_int(df, "A");
  out.C = column_as_int(df, "C");
  out.G = column_as_int(df, "G");
  out.T = column_as_int(df, "T");
  CharacterVector reference = as<CharacterVector>(df["reference"]);
  out.reference = as<std::vector<std::string>>(reference);
  return out;
}

double log_sum_exp(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  double m = std::max(a, b);
  return m + std::log1p(std::exp(-std::abs(a - b)));
}
