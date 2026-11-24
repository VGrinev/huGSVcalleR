#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;

double cpp_lchoose(int n, int k);
double cpp_choose(int n, int k);
double log_sum_exp(double a, double b);

#endif // UTILS_H
