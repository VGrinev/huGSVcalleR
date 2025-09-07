#ifndef SNP_FISHER_H
#define SNP_FISHER_Hs

#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <string>
#include "snp_entropy_cpp_parallel.h"
#include "snp_binomial_cpp_parallel.h"
#include "snp_fisher_cpp_parallel.h"
#include "snp_poisson_cpp_parallel.h"

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]

// Function to extract feastures for SNPs
DataFrame vectorise_data_cpp_parallel(DataFrame& df_d, CharacterVector features);

#endif // SNP_FISHER_H
