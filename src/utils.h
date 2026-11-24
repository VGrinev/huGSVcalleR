#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <RcppParallel.h>
#include <string>
#include <vector>

using namespace Rcpp;

double cpp_lchoose(int n, int k);
double cpp_choose(int n, int k);
double log_sum_exp(double a, double b);

struct BaseCounts {
  std::vector<int> A;
  std::vector<int> C;
  std::vector<int> G;
  std::vector<int> T;
  std::vector<std::string> reference;
};

BaseCounts read_base_counts(DataFrame& df);

#endif // UTILS_H
