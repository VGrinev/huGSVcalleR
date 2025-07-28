#include "snp_entropy_cpp_parallel.h"

struct SNPEntropyWorker : public Worker
{
  const RVector<int> A;
  const RVector<int> C;
  const RVector<int> G;
  const RVector<int> T;
  const bool p_value;
  const CharacterVector reference;
  RVector<double> results;

  const CharacterVector nucleotides;

  SNPEntropyWorker(const CharacterVector& reference, const IntegerVector& A,
                   const IntegerVector& C, const IntegerVector& G,
                   const IntegerVector& T, const bool& p_value_arg,
                   NumericVector& results)
    : A(A), C(C), G(G), T(T), p_value(p_value_arg), reference(reference), results(results),
      nucleotides(CharacterVector::create("A", "C", "G", "T")) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (std::size_t i = begin; i < end; ++i)
    {
     std::array<double, 4> SampleCite = {static_cast<double>(A[i]), static_cast<double>(C[i]),
                                         static_cast<double>(G[i]), static_cast<double>(T[i])};

     double NumOfCovers = SampleCite[0] + SampleCite[1] + SampleCite[2] + SampleCite[3];

      if (NumOfCovers == 0)
      {
        results[i] = 1.0;
        continue;
      }

      // Find the index of the maximum value
      int maxIndex = 0;
      double max = SampleCite[0];
      for (int j = 1; j < 4; ++j)
      {
       if (SampleCite[j] > max)
       {
         max = SampleCite[j];
         maxIndex = j;
        }
      }

      std::string refCite = Rcpp::as<std::string>(reference[i]);

      // Correct SampleCite if the reference does not match the maximum
      if (refCite != Rcpp::as<std::string>(nucleotides[maxIndex]))
      {
         double m = SampleCite[maxIndex];
         SampleCite[maxIndex] = m / 2;
         SampleCite[refCite == "A" ? 0 :
           refCite == "C" ? 1 :
           refCite == "G" ? 2 : 3] += m / 2;
      }

       // Normalize SampleCite and compute entropy
       double entropy = 0.0;
       for (int k = 0; k < 4; ++k)
       {
         if (SampleCite[k] > 0)
         {
           SampleCite[k] /= NumOfCovers;
           double log_value = log10(SampleCite[k]);
           if (std::isfinite(log_value))
           {
             entropy -= SampleCite[k] * log_value;
           }
         }
       }
       if (p_value) {
         // Compute chi-square statistics and p-value
         double chiStats = 2 * NumOfCovers * entropy;
         double p = R::pchisq(chiStats, 4, 0, false);
         results[i] = p;
       } else {
         results[i] = entropy;
       }
    }
  }
};

//' Parallel SNP Entropy and P-value Calculation
//' @name snp_entropy_cpp_parallel
//' @description
//' This function calculates the entropy and/or p-value for each SNP (single nucleotide polymorphism) position
//' in the provided genomic data. Entropy quantifies the variability in nucleotide distributions, while the
//' p-value assesses the statistical significance of the entropy using a chi-square test. The function supports
//' optional calculation of either entropy or p-values and utilizes parallel processing for efficiency on large datasets.
//'
//' @param df_d A data frame containing the following columns:
//'   - reference: A character vector specifying the reference allele for each position (values: "A", "C", "G", "T").
//'   - A, C, G, T: Integer vectors representing the observed counts of each nucleotide at each SNP position.
//' @param p_value A boolean flag (default: TRUE) indicating whether to calculate p-values. If FALSE, the function
//'   returns entropy values instead.
//'
//' @return
//' A numeric vector containing either:
//'   - Entropy values (if p_value = FALSE).
//'   - P-values of the chi-square tests for the entropy values (if p_value = TRUE).
//'   A p-value of 1 indicates no significant variation, while smaller values suggest greater variation and statistical significance.
//'
//' @details
//' The function performs the following steps:
//' 1. Retrieves the nucleotide counts (A, C, G, T) for each SNP position.
//' 2. Computes the total coverage (sum of counts for all nucleotides).
//' 3. If the total coverage is zero, assigns a p-value or entropy value of 1.
//' 4. Adjusts the nucleotide distribution if the reference allele does not match the most frequent nucleotide:
//'    - Splits the count of the most frequent nucleotide between it and the reference allele.
//' 5. Normalizes nucleotide counts to probabilities.
//' 6. Computes entropy using the formula:
//'    Entropy = -sum(p_i * log10(p_i))
//'    where p_i is the normalized count of nucleotide i.
//' 7. If p_value = TRUE, computes the chi-square statistic using the entropy and total coverage, then derives the p-value:
//'    chi^2 = 2 * total coverage * Entropy
//'    The p-value is calculated using the chi-square distribution with 4 degrees of freedom.
//' 8. Executes the calculations in parallel across multiple cores for faster processing of large datasets.
//'
//' @examples
//' # Example input data
//' df <- data.frame(reference = c("A", "C", "G", "T"),
//'                  A = c(10, 0, 5, 0),
//'                  C = c(2, 10, 0, 0),
//'                  G = c(0, 0, 10, 5),
//'                  T = c(0, 0, 0, 10))
//'
//' # Run SNP entropy and p-value calculation
//' result <- snp_entropy_cpp_parallel(df, p_value = TRUE)
//' print(result)
//'
//' # Run SNP entropy calculation only
//' entropy <- snp_entropy_cpp_parallel(df, p_value = FALSE)
//' print(entropy)
//'
// [[Rcpp::export]]
NumericVector snp_entropy_cpp_parallel(DataFrame& df_d, bool p_value = false)
{
  int n = df_d.nrows();
  NumericVector results(n);

  SNPEntropyWorker worker(df_d["reference"], df_d["A"], df_d["C"],
                           df_d["G"], df_d["T"], p_value, results);
  parallelFor(0, n, worker);

  return results;
}
