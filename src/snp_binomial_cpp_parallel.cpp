#include "snp_binomial_cpp_parallel.h"
#include "utils.h"

struct SNPBinomialWorker : public Worker
{
  const std::vector<int>& A;
  const std::vector<int>& C;
  const std::vector<int>& G;
  const std::vector<int>& T;
  const std::vector<std::string>& reference;
  const double p;
  RVector<double> binom;

  SNPBinomialWorker(const std::vector<int>& A, const std::vector<int>& C,
                    const std::vector<int>& G, const std::vector<int>& T,
                    const std::vector<std::string>& reference,
                    const double& p, NumericVector& binom)
    : A(A), C(C), G(G), T(T), reference(reference), p(p), binom(binom) {}

  void operator()(size_t begin, size_t end)
  {
    for (size_t j = begin; j < end; ++j)
    {

      int total = A[j] + C[j] + G[j] + T[j];
      int ref_count = A[j] * (reference[j] == "A") +
        C[j] * (reference[j] == "C") +
        G[j] * (reference[j] == "G") +
        T[j] * (reference[j] == "T");

      int cm = total - ref_count;
      double log_p = std::log(p);
      double log_q = std::log(1.0 - p);
      double log_sum = R_NegInf;

      for (int k = 0; k <= total; ++k)
      {
        if (cm <= k)
        {
          double log_temp = cpp_lchoose(total, k) + k * log_p +
            (total - k) * log_q;
          log_sum = log_sum_exp(log_sum, log_temp);
        }
      }
      binom[j] = std::exp(log_sum);
    }
  }
};

//' Parallel Calculation of Binomial Probabilities for SNPs
//' @name snp_binomial_cpp_parallel
//' @description
//' This function calculates the cumulative binomial probability for each
//' SNP in the provided dataset. It takes into account the reference
//' allele at each position and computes the binomial probability for
//' observing a certain number of non-reference alleles, given a
//' probability of mutation `p`. The function utilizes parallel
//' computation to speed up the process for large datasets.
//'
//' @param df_d A data frame with the following columns:
//'   - `A`, `C`, `G`, `T`: Integer vectors representing the counts of
//'     each nucleotide at each position.
//'   - `reference`: A character vector containing the reference alleles
//'     ("A", "C", "G", "T") at each position.
//'
//' @return A numeric vector of binomial probabilities, where each
//'   element corresponds to a SNP and represents the cumulative binomial
//'   probability of observing the reference allele count given a
//'   mutation probability `p`.
//'
//' @details
//' The function works by first determining the total number of alleles
//' (`A[i]`, `C[i]`, `G[i]`, `T[i]`) for each SNP, and calculating how
//' many of those are the reference allele. The binomial probability is
//' computed for different possible outcomes, considering the mutation
//' probability `p` and the number of non-reference alleles. The
//' calculation is performed in parallel for faster execution on large
//' datasets.
//'
//' The binomial probability for each SNP is calculated as the sum of
//' binomial probabilities for each possible count of non-reference
//' alleles (`k`), using the binomial coefficient formula and the
//' probability mass function for the binomial distribution.
//'
//' @author Dzianis D. Sarnatski, Mikalai M. Yatskou
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector snp_binomial_cpp_parallel(DataFrame& df_d)
{
  int n = df_d.nrows();

  double p = 1e-5;
  NumericVector binom(n);

  BaseCounts counts = read_base_counts(df_d);
  SNPBinomialWorker worker(counts.A, counts.C, counts.G, counts.T,
                           counts.reference, p, binom);
  parallelFor(0, n, worker);

  return binom;
}
