#ifndef SNP_ENTROPY_H
#define SNP_ENTROPY_H

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]

// Function to compute entropy for SNPs
NumericVector snp_entropy_cpp_parallel(DataFrame& df_d, bool p_value);

#endif // SNP_ENTROPY_H