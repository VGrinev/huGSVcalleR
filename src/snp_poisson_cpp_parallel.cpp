#include "snp_poisson_cpp_parallel.h"
#include "utils.h"

 // The Worker framework for parallel computing of the Poisson test
 struct SNPPoissonWorker : public Worker {
  const CharacterVector reference;
  const RVector<int> A;
  const RVector<int> C;
  const RVector<int> G;
  const RVector<int> T;

  RVector<double> results;

  SNPPoissonWorker(const CharacterVector& reference,
              const IntegerVector& A, const IntegerVector& C,
              const IntegerVector& G, const IntegerVector& T,
              const NumericVector& results)
    : reference(reference), A(A), C(C), G(G), T(T),
      results(results) {}

  void operator()(size_t begin, size_t end) {
    double p = 9.0 / 20000.0; // TODO: Clarify
    for (size_t j = begin; j < end; ++j) {
      int N = A[j] + C[j] + G[j] +T[j];  // Сумма A, C, G, T
      std::string ref = as<std::string>(reference[j]);
      int ref_count = A[j] * (reference[j] == "A") +
        C[j] * (reference[j] == "C") +
        G[j] * (reference[j] == "G") +
        T[j] * (reference[j] == "T");

      int K = N - ref_count;

      double Pr = 0.0;
      for (int k = K; k <= N; ++k) {
        Pr += cpp_choose(N, k) * pow(p, k) * pow(1 - p, N - k);
      }

      results[j] = Pr;  // Saving the result for the current SNP
    }
  }
};

//' Parallel Poisson Test for SNP Detection
//' @name snp_poisson_cpp_parallel
//' @description
//' This function calculates the Poisson-based probability for detecting SNPs at each position in genomic data.
//' The probability is computed for observing the non-reference allele count under a Poisson distribution.
//' Parallel processing is employed to accelerate computation across large datasets.
//'
//' @param df_d A data frame containing the following columns:
//'   - `reference`: A character vector with the reference alleles ("A", "C", "G", "T").
//'   - `A`, `C`, `G`, `T`: Integer vectors representing the counts of each nucleotide at each position.
//'
//' @return A numeric vector containing the Poisson-based probability for each position.
//'         This probability represents the likelihood of observing the non-reference allele count
//'         under the Poisson model with a fixed mutation rate (default: 9/20000).
//'
//' @details
//' The function performs the following steps:
//' 1. For each position, the total number of nucleotide counts (`A[i]`, `C[i]`, `G[i]`, `T[i]`) is computed.
//' 2. The count of the reference allele is determined based on the `reference` column.
//' 3. The non-reference allele count is calculated as the total count minus the reference allele count.
//' 4. Using the Poisson distribution, the probability of observing at least the non-reference allele count is computed.
//' 5. The calculations are performed in parallel for efficiency on large datasets.
//'
//' @note The mutation rate used in the Poisson model is fixed at 9/20000 by default.
//'
//' @examples
//' df <- data.frame(reference = c("A", "C", "G", "T"),
//'                  A = c(10, 0, 5, 0),
//'                  C = c(2, 10, 0, 0),
//'                  G = c(0, 0, 10, 5),
//'                  T = c(0, 0, 0, 10))
//'
//' result <- snp_poisson_cpp_parallel(df)
//' print(result)
//'
// [[Rcpp::export]]
NumericVector snp_poisson_cpp_parallel(DataFrame& df_d) {
  int n = df_d.nrows();

  NumericVector results(n);  // Vector of results

  SNPPoissonWorker worker(df_d["reference"], df_d["A"], df_d["C"],
                          df_d["G"], df_d["T"], results);

  // Parallel processing
  parallelFor(0, n, worker);

  return results;  // Returning the probability vector
}


