#include "snp_counting_cpp_parallel.h"

 struct SNPCounting : public Worker
 {
   const CharacterVector reference;
   const RVector<int> A;
   const RVector<int> C;
   const RVector<int> G;
   const RVector<int> T;

   RVector<double> count;

   SNPCounting(const CharacterVector& reference,
               const IntegerVector& A, const IntegerVector& C,
               const IntegerVector& G, const IntegerVector& T,
               const NumericVector& count)
     : reference(reference), A(A), C(C), G(G), T(T),
       count(count) {}

   void operator()(size_t begin, size_t end)
   {
     for (size_t i = begin; i < end; ++i)
     {
       int i_A = (reference[i] == "A");
       int i_C = (reference[i] == "C");
       int i_G = (reference[i] == "G");
       int i_T = (reference[i] == "T");

       int total = A[i] + C[i] + G[i] + T[i];
       int ref_count = A[i] * i_A + C[i] * i_C + G[i] * i_G + T[i] * i_T;

       if (total != 0) {  // To avoid zero division
         count[i] = static_cast<double>(ref_count) / total;
       } else {
         count[i] = NA_REAL;  // If total == 0, return NA
       }
     }
   }
 };

//' Parallel SNP Counting in a DataFrame
//' @name snp_counting_cpp_parallel
//' @description
//' This function calculates the reference allele count ratio for each position in the provided genomic data.
//' It uses parallel processing for faster execution by distributing the computation across multiple cores.
//' The ratio is computed as the reference allele count divided by the total allele count at each position.
//' The function uses the Rcpp and RcppParallel libraries to efficiently handle the calculations in parallel.
//'
//' @param df_d A data frame containing the following columns:
//'   - `reference`: A character vector with the reference alleles ("A", "C", "G", "T").
//'   - `A`, `C`, `G`, `T`: Integer vectors representing the counts of each nucleotide at each position.
//'
//' @return A numeric vector containing the ratio of the reference allele count to the total allele count for each position.
//'   If the total count is zero, the result will be `NA`.
//'
//' @details
//' The function performs the following steps:
//' 1. For each position, the total number of nucleotide counts (`A[i]`, `C[i]`, `G[i]`, `T[i]`) is computed.
//' 2. The reference allele count is determined based on the `reference` column, comparing it to the nucleotide counts.
//' 3. The ratio of the reference allele count to the total count is calculated for each position.
//' 4. The calculation is performed in parallel across multiple cores to speed up the process for large datasets.
//'
//' @author Dzianis D. Sarnatski, Mikalai M. Yatskou
//' @examples
//' df <- data.frame(reference = c("A", "C", "G", "T"),
//'                  A = c(10, 0, 5, 0),
//'                  C = c(2, 10, 0, 0),
//'                  G = c(0, 0, 10, 5),
//'                  T = c(0, 0, 0, 10))
//'
//' result <- snp_counting_cpp_parallel(df)
//' print(result)
//'
// [[Rcpp::export]]
NumericVector snp_counting_cpp_parallel(DataFrame& df_d)
{
  int n = df_d.nrows();

  NumericVector count(n);

  SNPCounting snpCounting(df_d["reference"], df_d["A"], df_d["C"],
                          df_d["G"], df_d["T"], count);

  parallelFor(0, n, snpCounting);

  return count;
}
