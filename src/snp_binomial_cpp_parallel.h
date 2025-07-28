#ifndef SNP_BINOMIAL_H
#define SNP_BINOMIAL_H

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]

// Function to compute cumulative binomial probabilities for SNPs
NumericVector snp_binomial_cpp_parallel(DataFrame& df_d);

#endif // SNP_BINOMIAL_H
