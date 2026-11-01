#include <Rcpp.h>
using namespace Rcpp;

// Logarithm of the binomial coefficient C(n, k) via lgamma.
// lgamma(x) returns log(Gamma(x)) = log((x-1)!), so:
// log(C(n, k)) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
double cpp_lchoose(int n, int k) {
  if (k < 0 || k > n) return R_NegInf;
  if (k == 0 || k == n) return 0.0;
  return R::lgammafn(n + 1.0) - R::lgammafn(k + 1.0) - R::lgammafn(n - k + 1.0);
}

// This function calculates the binomial coefficient C(n, k),
// which represents the number of ways to choose k objects from a set of n objects without regard to order.
// The function uses an optimized approach to calculate the binomial coefficient for faster execution,
// leveraging the symmetry property C(n, k) = C(n, n-k) to reduce the number of operations.
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
