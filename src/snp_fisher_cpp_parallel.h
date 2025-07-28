#ifndef SNP_FISHER_H
#define SNP_FISHER_Hs

#include <Rcpp.h>
#include <RcppParallel.h>
#include <random>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]

// Function to compute fisher for SNPs
NumericVector snp_fisher_cpp_parallel(DataFrame& df_d);

#endif // SNP_FISHER_H
