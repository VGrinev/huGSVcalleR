#include "snp_fisher_cpp_parallel.h"
#include "utils.h"

struct FisherWorker : public Worker {
  const RVector<int> A;
  const RVector<int> C;
  const RVector<int> G;
  const RVector<int> T;
  const CharacterVector reference;

  RVector<int> cn, ct, cr, cm, ctr, ctr_max;

  FisherWorker(const IntegerVector& A, const IntegerVector& C, const IntegerVector& G,
               const IntegerVector& T, const CharacterVector& reference,
               IntegerVector& cn, IntegerVector& ct,
               IntegerVector& cr, IntegerVector& cm, IntegerVector& ctr, IntegerVector& ctr_max)
    : A(A), C(C), G(G), T(T), reference(reference),
      cn(cn), ct(ct), cr(cr), cm(cm), ctr(ctr), ctr_max(ctr_max) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    static thread_local std::mt19937 gen(std::random_device{}());

    for (std::size_t j = begin; j < end; ++j)
    {
      int total = A[j] + C[j] + G[j] + T[j];
      int ref_count = (reference[j] == "A") * A[j] +
        (reference[j] == "C") * C[j] +
        (reference[j] == "G") * G[j] +
        (reference[j] == "T") * T[j];

      std::uniform_real_distribution<> dis1(0.0, total % 2 + 1);
      cn[j] = static_cast<int>(std::min(total / 2 + dis1(gen), 1.25 * ref_count + 1));
      ct[j] = total - cn[j];
      cr[j] = ref_count;
      cm[j] = total - ref_count;

      std::uniform_real_distribution<> dis2(0.8 * cm[j], cm[j] + 1.0);
      ctr[j] = ct[j] - std::min(static_cast<int>(dis2(gen)), ct[j]);
      ctr_max[j] = ctr[j];
    }
  }
};

struct FisherCalcWorker : public Worker
{
  const RVector<int> cn;
  const RVector<int> ct;
  const RVector<int> cr;
  const RVector<int> cm;
  const RVector<int> ctr_max;

  RVector<int> ctr;

  RVector<double> fisher;

  FisherCalcWorker(const IntegerVector& cn, const IntegerVector& ct,
                   const IntegerVector& cr, const IntegerVector& cm,
                   const IntegerVector& ctr_max, const IntegerVector& ctr,
                   NumericVector& fisher)
    : cn(cn), ct(ct), cr(cr), cm(cm), ctr_max(ctr_max), ctr(ctr), fisher(fisher) {}

  void operator()(std::size_t begin, std::size_t end)
  {
    for (std::size_t j = begin; j < end; ++j)
    {
      fisher[j] = 0; // Инициализируем значение fisher
      for (int i = 0; i <= ctr_max[j]; ++i)
      {
        if (ctr[j] >= i)
        {
          double temp = cpp_choose(ct[j], i) * cpp_choose(cn[j], cr[j] - i) / cpp_choose(ct[j] + cm[j], cm[j]);
          fisher[j] += temp;
        }
      }
    }
  }
};

//' Parallel SNP Fisher's Exact Test Calculation
//' @name snp_fisher_cpp_parallel
//' @description
//' This function performs a parallel computation of Fisher's Exact Test for SNP (single nucleotide polymorphism) data.
//' It calculates the Fisher's exact test statistic based on nucleotide counts (A, C, G, T) at each SNP position.
//' The function utilizes parallel processing to efficiently handle large datasets by distributing the computation
//' across multiple CPU cores using RcppParallel.
//'
//' @param df_d A data frame containing the following columns:
//'   - `reference`: A character vector with the reference alleles ("A", "C", "G", "T").
//'   - `A`, `C`, `G`, `T`: Integer vectors representing the counts of each nucleotide (A, C, G, T) at each SNP position.
//'
//' @return A numeric vector containing the p-values from Fisher's Exact Test for each SNP position.
//'   A p-value of 1 indicates no variation, while smaller values suggest greater variation and statistical significance.
//'
//' @details
//' The function performs the following steps:
//' 1. For each SNP position, the nucleotide counts for A, C, G, T are extracted from the input data frame.
//' 2. Intermediate quantities such as `cn`, `ct`, `cr`, and others are computed using parallel processing.
//' 3. Fisher's exact test statistics are calculated for each SNP position using binomial coefficients.
//' 4. The function is optimized for parallel execution, distributing the workload across multiple cores.
//'
//' @author Dzianis D. Sarnatski, Mikalai M. Yatskou
//' @examples
//' df <- data.frame(reference = c("A", "C", "G", "T"),
//'                  A = c(10, 0, 5, 0),
//'                  C = c(2, 10, 0, 0),
//'                  G = c(0, 0, 10, 5),
//'                  T = c(0, 0, 0, 10))
//'
//' result <- snp_fisher_cpp_parallel(df)
//' print(result)
//'
// [[Rcpp::export]]
NumericVector snp_fisher_cpp_parallel(DataFrame& df_d)
{

  int n = df_d.nrows();
  NumericVector fisher(n);
  IntegerVector cn(n), ct(n), cr(n), cm(n), ctr(n), ctr_max(n);

  // 1. Параллельные вычисления для cn, ct, cr, cm, ctr
  FisherWorker worker(df_d["A"], df_d["C"], df_d["G"], df_d["T"],
                      df_d["reference"], cn, ct, cr, cm, ctr, ctr_max);
  parallelFor(0, n, worker);

  // 2. Параллельные вычисления для fisher
  FisherCalcWorker fisherCalcWorker(cn, ct, cr, cm, ctr, ctr_max, fisher);
  parallelFor(0, n, fisherCalcWorker);

  return fisher;
}
