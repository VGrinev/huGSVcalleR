#include "snp_fisher_cpp_parallel.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <limits>

// Two-sided Fisher's exact test
// The sample row is the reference count against every other base
// The control row has the same depth and only the reference base
// The denominator is the ordinary hypergeometric C(n, column)
static double fisher_twosided(int ref_count, int alt_count) {
  if (alt_count <= 0) {
    return 1.0;
  }
  int total = ref_count + alt_count;
  int r1 = total;
  int r2 = total;
  int k1 = ref_count + total;
  int n = r1 + r2;
  int lo = std::max(0, r1 + k1 - n);
  int hi = std::min(r1, k1);
  auto log_prob = [&](int x) {
    return cpp_lchoose(r1, x) + cpp_lchoose(r2, k1 - x) - cpp_lchoose(n, k1);
  };
  double log_obs = log_prob(ref_count);
  double log_sum = R_NegInf;
  for (int x = lo; x <= hi; ++x) {
    double log_x = log_prob(x);
    if (log_x <= log_obs + 1e-7) {
      log_sum = log_sum_exp(log_sum, log_x);
    }
  }
  if (!std::isfinite(log_sum)) {
    return 1.0;
  }
  double p = std::exp(log_sum);
  if (!std::isfinite(p) || p > 1.0) {
    return 1.0;
  }
  // exp underflows to 0 long before the probability is actually zero
  if (p == 0.0) {
    return std::numeric_limits<double>::min();
  }
  return p;
}

struct FisherWorker : public Worker {
  const std::vector<int>& A;
  const std::vector<int>& C;
  const std::vector<int>& G;
  const std::vector<int>& T;
  const std::vector<std::string>& reference;
  RVector<double> fisher;

  FisherWorker(const std::vector<int>& A, const std::vector<int>& C,
               const std::vector<int>& G, const std::vector<int>& T,
               const std::vector<std::string>& reference,
               NumericVector& fisher)
    : A(A), C(C), G(G), T(T), reference(reference), fisher(fisher) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      int total = A[j] + C[j] + G[j] + T[j];
      int ref_count = (reference[j] == "A") * A[j] +
        (reference[j] == "C") * C[j] +
        (reference[j] == "G") * G[j] +
        (reference[j] == "T") * T[j];
      fisher[j] = fisher_twosided(ref_count, total - ref_count);
    }
  }
};

//' Parallel SNP Fisher's Exact Test Calculation
//' @name snp_fisher_cpp_parallel
//' @description
//' This function performs a parallel computation of Fisher's Exact Test
//' for SNP (single nucleotide polymorphism) data. It calculates the
//' Fisher's exact test statistic based on nucleotide counts (A, C, G, T)
//' at each SNP position. The function utilizes parallel processing to
//' efficiently handle large datasets by distributing the computation
//' across multiple CPU cores using RcppParallel.
//'
//' @param df_d A data frame containing the following columns:
//'   - `reference`: A character vector with the reference alleles
//'     ("A", "C", "G", "T").
//'   - `A`, `C`, `G`, `T`: Integer vectors representing the counts of
//'     each nucleotide (A, C, G, T) at each SNP position.
//'
//' @return A numeric vector containing the p-values from Fisher's Exact
//'   Test for each SNP position. A p-value of 1 indicates no variation,
//'   while smaller values suggest greater variation and statistical
//'   significance.
//'
//' @details
//' The sample row is the reference count versus all other bases. The
//' control row has the same depth and contains only the reference base
//' The p-value is the two-sided hypergeometric probability, with the
//' ordinary denominator C(n, column total). It does not depend on a
//' random draw
//'
//' @author Dzianis D. Sarnatski, Mikalai M. Yatskou
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector snp_fisher_cpp_parallel(DataFrame& df_d)
{
  int n = df_d.nrows();
  NumericVector fisher(n);
  BaseCounts counts = read_base_counts(df_d);
  FisherWorker worker(counts.A, counts.C, counts.G, counts.T,
                      counts.reference, fisher);
  parallelFor(0, n, worker);
  return fisher;
}
