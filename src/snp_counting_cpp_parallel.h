#ifndef SNP_COUNTING_H
#define SNP_COUNTING_H

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]

// Function to compute counting criterions for SNPs
NumericVector snp_counting_cpp_parallel(DataFrame& df_d);

#endif // SNP_COUNTING_H