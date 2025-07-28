#include <Rcpp.h>
using namespace Rcpp;

// This function calculates the binomial coefficient C(n, k),
// which represents the number of ways to choose k objects from a set of n objects without regard to order.
// The function uses an optimized approach to calculate the binomial coefficient for faster execution,
// leveraging the symmetry property C(n, k) = C(n, n-k) to reduce the number of operations.
double cpp_choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  if (k == 0 || k == n) return 1;
  
  // Оптимизированный расчет биномиального коэффициента
  k = std::min(k, n - k); // Используем симметрию C(n, k) = C(n, n-k)
  double result = 1.0;
  for (int i = 1; i <= k; ++i) {
    result *= (n - i + 1) / static_cast<double>(i);
  }
  return result;
}