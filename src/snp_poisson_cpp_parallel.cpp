#include "snp_poisson_cpp_parallel.h"
#include "utils.h"

 // The Worker framework for parallel computing of the Poisson test
 struct SNPPoissonWorker : public Worker {
  const std::vector<std::string>& reference;
  const std::vector<int>& A;
  const std::vector<int>& C;
  const std::vector<int>& G;
  const std::vector<int>& T;

  RVector<double> results;

  SNPPoissonWorker(const std::vector<std::string>& reference,
              const std::vector<int>& A, const std::vector<int>& C,
              const std::vector<int>& G, const std::vector<int>& T,
              NumericVector& results)
    : reference(reference), A(A), C(C), G(G), T(T),
      results(results) {}

  void operator()(size_t begin, size_t end) {
    double p = 1.0 / 1000.0;
    for (size_t j = begin; j < end; ++j) {
      int N = A[j] + C[j] + G[j] +T[j];
      int ref_count = A[j] * (reference[j] == "A") +
        C[j] * (reference[j] == "C") +
        G[j] * (reference[j] == "G") +
        T[j] * (reference[j] == "T");

      int K = N - ref_count;

      double log_p = std::log(p);
      double log_q = std::log(1.0 - p);
      double log_sum = R_NegInf;
      for (int k = K; k <= N; ++k) {
        double log_temp = cpp_lchoose(N, k) + k * log_p +
          (N - k) * log_q;
        log_sum = log_sum_exp(log_sum, log_temp);
      }

      results[j] = std::exp(log_sum);
    }
  }
};

//' Parallel Poisson Test for SNP Detection
//' @name snp_poisson_cpp_parallel
//' @description
//' This function calculates the Poisson-based probability for detecting
//' SNPs at each position in genomic data. The probability is computed
//' for observing the non-reference allele count under a Poisson
//' distribution. Parallel processing is employed to accelerate
//' computation across large datasets.
//'
//' @param df_d A data frame containing the following columns:
//'   - `reference`: A character vector with the reference alleles
//'     ("A", "C", "G", "T").
//'   - `A`, `C`, `G`, `T`: Integer vectors representing the counts of
//'     each nucleotide at each position.
//'
//' @return A numeric vector containing the Poisson-based probability for
//'         each position. This probability represents the likelihood of
//'         observing the non-reference allele count under the Poisson
//'         model with a fixed mutation rate of 1/1000
//'
//' @details
//' The function performs the following steps:
//' 1. For each position, the total number of nucleotide counts
//'    (`A[i]`, `C[i]`, `G[i]`, `T[i]`) is computed.
//' 2. The count of the reference allele is determined based on the
//'    `reference` column.
//' 3. The non-reference allele count is calculated as the total count
//'    minus the reference allele count.
//' 4. Using the Poisson distribution, the probability of observing at
//'    least the non-reference allele count is computed.
//' 5. The calculations are performed in parallel for efficiency on large
//'    datasets.
//'
//' @note The mutation rate used in the Poisson model is fixed at 1/1000
//' The reported probability is the binomial upper tail at that rate
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector snp_poisson_cpp_parallel(DataFrame& df_d)
{
  int n = df_d.nrows();

  NumericVector results(n);  // Vector of results

  BaseCounts counts = read_base_counts(df_d);
  SNPPoissonWorker worker(counts.reference, counts.A, counts.C,
                          counts.G, counts.T, results);

  // Parallel processing
  parallelFor(0, n, worker);

  return results;  // Returning the probability vector
}


