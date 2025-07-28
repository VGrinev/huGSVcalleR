#include "snp_binomial_cpp_parallel.h"
#include "utils.h"

struct SNPBinomialWorker : public Worker
{
  const RVector<int> A;
  const RVector<int> C;
  const RVector<int> G;
  const RVector<int> T;
  const CharacterVector reference;
  const double p;
  RVector<double> binom;

  SNPBinomialWorker(const IntegerVector& A, const IntegerVector& C, const IntegerVector& G,
                    const IntegerVector& T, const CharacterVector& reference,
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
      binom[j] = 0;

      for (int k = 0; k <= total; ++k)
      {
        if (total >= k && cm <= k)
        {
          double temp = cpp_choose(total, k) * pow(p, k) * pow(1 - p, total - k);
          binom[j] += temp;
        }
      }
    }
  }
};

//' Parallel Calculation of Binomial Probabilities for SNPs
//' @name snp_binomial_cpp_parallel
//' @description
//' This function calculates the cumulative binomial probability for each SNP in the provided dataset.
//' It takes into account the reference allele at each position and computes the binomial probability
//' for observing a certain number of non-reference alleles, given a probability of mutation `p`.
//' The function utilizes parallel computation to speed up the process for large datasets.
//'
//' @param df_d A data frame with the following columns:
//'   - `A`, `C`, `G`, `T`: Integer vectors representing the counts of each nucleotide at each position.
//'   - `reference`: A character vector containing the reference alleles ("A", "C", "G", "T") at each position.
//'
//' @return A numeric vector of binomial probabilities, where each element corresponds to a SNP
//'   and represents the cumulative binomial probability of observing the reference allele count
//'   given a mutation probability `p`.
//'
//' @details
//' The function works by first determining the total number of alleles (`A[i]`, `C[i]`, `G[i]`, `T[i]`)
//' for each SNP, and calculating how many of those are the reference allele. The binomial probability is
//' computed for different possible outcomes, considering the mutation probability `p` and the number of
//' non-reference alleles. The calculation is performed in parallel for faster execution on large datasets.
//'
//' The binomial probability for each SNP is calculated as the sum of binomial probabilities for each
//' possible count of non-reference alleles (`k`), using the binomial coefficient formula and the probability
//' mass function for the binomial distribution.
//'
//' @author Dzianis D. Sarnatski, Mikalai M. Yatskou
//' @examples
//' df <- data.frame(reference = c("A", "C", "G", "T"),
//'                  A = c(10, 0, 5, 0),
//'                  C = c(2, 10, 0, 0),
//'                  G = c(0, 0, 10, 5),
//'                  T = c(0, 0, 0, 10))
//'
//' result <- snp_binomial_cpp_parallel(df)
//' print(result)
// [[Rcpp::export]]
NumericVector snp_binomial_cpp_parallel(DataFrame& df_d)
{
  int n = df_d.nrows();

  double p = 1e-5;
  NumericVector binom(n);

  // Creating a Worker for parallel calculations
  SNPBinomialWorker worker(df_d["A"], df_d["C"], df_d["G"], df_d["T"],
                           df_d["reference"], p, binom);
  parallelFor(0, n, worker);

  return binom;
}
