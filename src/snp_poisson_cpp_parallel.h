#ifndef SNP_POISSON_H
#define SNP_POISSON_H

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]

// Function to compute Poisson criterions for SNPs
NumericVector snp_poisson_cpp_parallel(DataFrame& df_d);

#endif // SNP_POISSON_H