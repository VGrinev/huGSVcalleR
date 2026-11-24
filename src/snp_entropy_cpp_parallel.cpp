#include "snp_entropy_cpp_parallel.h"
#include "utils.h"
#include <cmath>

struct SNPEntropyWorker : public Worker
{
  const std::vector<int>& A;
  const std::vector<int>& C;
  const std::vector<int>& G;
  const std::vector<int>& T;
  const bool p_value;
  const std::vector<std::string>& reference;
  RVector<double> results;
  std::vector<double>& chi;

  SNPEntropyWorker(const std::vector<std::string>& reference,
                   const std::vector<int>& A, const std::vector<int>& C,
                   const std::vector<int>& G, const std::vector<int>& T,
                   const bool& p_value_arg, NumericVector& results,
                   std::vector<double>& chi)
    : A(A), C(C), G(G), T(T), p_value(p_value_arg), reference(reference),
      results(results), chi(chi) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (std::size_t i = begin; i < end; ++i)
    {
     std::array<double, 4> SampleCite = {static_cast<double>(A[i]),
                                         static_cast<double>(C[i]),
                                         static_cast<double>(G[i]),
                                         static_cast<double>(T[i])};

     double NumOfCovers = SampleCite[0] +
       SampleCite[1] +
       SampleCite[2] +
       SampleCite[3];

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

      static const char* nucleotides[] = {"A", "C", "G", "T"};
      std::string refCite = reference[i];

      // Correct SampleCite if the reference does not match the maximum
      if (refCite != nucleotides[maxIndex])
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
           double log_value = std::log(SampleCite[k]);
           if (std::isfinite(log_value))
           {
             entropy -= SampleCite[k] * log_value;
           }
         }
       }
       if (p_value) {
         chi[i] = 2 * NumOfCovers * entropy;
       } else {
         results[i] = entropy;
       }
    }
  }
};

//' Parallel SNP Entropy and P-value Calculation
//' @name snp_entropy_cpp_parallel
//' @description
//' This function calculates the entropy and/or p-value for each SNP
//' (single nucleotide polymorphism) position in the provided genomic
//' data. Entropy quantifies the variability in nucleotide distributions,
//' while the p-value assesses the statistical significance of the
//' entropy using a chi-square test. The function supports optional
//' calculation of either entropy or p-values and utilizes parallel
//' processing for efficiency on large datasets.
//'
//' @param df_d A data frame containing the following columns:
//'   - reference: A character vector specifying the reference allele for
//'     each position (values: "A", "C", "G", "T").
//'   - A, C, G, T: Integer vectors representing the observed counts of
//'     each nucleotide at each SNP position.
//' @param p_value A boolean flag (default: TRUE) indicating whether to
//'   calculate p-values. If FALSE, the function returns entropy values
//'   instead.
//'
//' @return
//' A numeric vector containing either:
//'   - Entropy values (if p_value = FALSE).
//'   - P-values of the chi-square tests for the entropy values (if
//'     p_value = TRUE).
//'   A p-value of 1 indicates no significant variation, while smaller
//'   values suggest greater variation and statistical significance.
//'
//' @details
//' The function performs the following steps:
//' 1. Retrieves the nucleotide counts (A, C, G, T) for each SNP position.
//' 2. Computes the total coverage (sum of counts for all nucleotides).
//' 3. If the total coverage is zero, assigns a p-value or entropy value
//'    of 1.
//' 4. Adjusts the nucleotide distribution if the reference allele does
//'    not match the most frequent nucleotide:
//'    - Splits the count of the most frequent nucleotide between it and
//'      the reference allele.
//' 5. Normalizes nucleotide counts to probabilities.
//' 6. Computes entropy using the formula:
//'    Entropy = -sum(p_i * log(p_i))
//'    where p_i is the normalized count of nucleotide i.
//' 7. If p_value = TRUE, computes the chi-square statistic using the
//'    entropy and total coverage, then derives the p-value:
//'    chi^2 = 2 * total coverage * Entropy
//'    The p-value is calculated using the chi-square distribution with
//'    3 degrees of freedom.
//' 8. Executes the calculations in parallel across multiple cores for
//'    faster processing of large datasets.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector snp_entropy_cpp_parallel(DataFrame& df_d, bool p_value = false)
{
  int n = df_d.nrows();
  NumericVector results(n);
  std::vector<double> chi(n, NA_REAL);
  BaseCounts counts = read_base_counts(df_d);

  SNPEntropyWorker worker(counts.reference, counts.A, counts.C,
                           counts.G, counts.T, p_value, results, chi);
  parallelFor(0, n, worker);

  if (p_value) {
    for (int i = 0; i < n; ++i) {
      if (R_IsNA(chi[i])) {
        continue;
      }
      results[i] = R::pchisq(chi[i], 3, 0, false);
    }
  }

  return results;
}
