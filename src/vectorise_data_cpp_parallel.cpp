#include "vectorise_data_cpp_parallel.h"

bool check_match_with_non_ref_quantity(const std::vector<int>& counts, int ref_index) {
  int ref_reads = counts[ref_index];
  int non_ref_reads = std::accumulate(counts.begin(), counts.end(), 0) - ref_reads;
  return (ref_reads == non_ref_reads) ||
    (ref_reads + 1 == non_ref_reads) ||
    (ref_reads - 1 == non_ref_reads);
}

bool check_match_with_max_coverage(const std::vector<int>& counts, int ref_index) {
  int ref_reads = counts[ref_index];
  int max_reads = *std::max_element(counts.begin(), counts.end());
  return ref_reads == max_reads;
}

double read_position_mean(const std::vector<int>& counts) {
  std::vector<int> positions = {1, 2, 3, 4};
  int total_reads = std::accumulate(counts.begin(), counts.end(), 0);
  if (total_reads == 0) return 0.0;

  double weighted_sum = 0.0;
  for (size_t i = 0; i < counts.size(); ++i) {
    weighted_sum += positions[i] * counts[i];
  }

  return weighted_sum / total_reads;
}

double read_position_variance(const std::vector<int>& counts) {
  std::vector<int> positions = {1, 2, 3, 4};
  std::vector<int> variant_positions;

  for (size_t i = 0; i < counts.size(); ++i) {
    if (counts[i] > 0) {
      variant_positions.push_back(positions[i]);
    }
  }

  if (variant_positions.size() <= 1) return 0.0;

  double mean = std::accumulate(variant_positions.begin(), variant_positions.end(), 0.0) / variant_positions.size();

  double variance = 0.0;
  for (int pos : variant_positions) {
    variance += std::pow(pos - mean, 2);
  }

  return variance / (variant_positions.size() - 1);
}

double log_error_prob(const std::vector<int>& counts, const int& ref_index) {
  const double p1 = 0.05;
  const double p2 = 0.5;

  int ref_reads = counts[ref_index];
  int non_ref_reads = std::accumulate(counts.begin(), counts.end(), 0) - ref_reads;
  int total_reads = std::accumulate(counts.begin(), counts.end(), 0);

  double prob1 = std::exp(std::lgamma(total_reads + 1) - std::lgamma(non_ref_reads + 1) - std::lgamma(total_reads - non_ref_reads + 1) + non_ref_reads * std::log(p1) + (total_reads - non_ref_reads) * std::log(1 - p1));
  double prob2 = std::exp(std::lgamma(total_reads + 1) - std::lgamma(non_ref_reads + 1) - std::lgamma(total_reads - non_ref_reads + 1) + non_ref_reads * std::log(p2) + (total_reads - non_ref_reads) * std::log(1 - p2));

  return std::log(prob1 / prob2);
}

struct VectoriseWorker : public Worker {
  const RVector<int> A;
  const RVector<int> C;
  const RVector<int> G;
  const RVector<int> T;
  const CharacterVector reference;
  RMatrix<double> results;
  std::unordered_set<std::string> selected_features;

  VectoriseWorker(const IntegerVector& A, const IntegerVector& C,
                  const IntegerVector& G, const IntegerVector& T,
                  const CharacterVector& reference, NumericMatrix& results,
                  const std::unordered_set<std::string>& selected_features)
    : A(A), C(C), G(G), T(T),
      reference(reference), results(results),
      selected_features(selected_features) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      std::vector<int> counts(4);
      counts[0] = A[i];
      counts[1] = C[i];
      counts[2] = G[i];
      counts[3] = T[i];

      std::string ref = as<std::string>(reference[i]);
      int ref_index = (ref == "A" ? 0 : ref == "C" ? 1 : ref == "G" ? 2 : 3);

      size_t col = 0;

      if (selected_features.count("x1")) {
        results(i, col++) = counts[ref_index];
      }

      if (selected_features.count("x2") || selected_features.count("x3") || selected_features.count("x4")) {
        std::vector<int> non_ref_counts = counts;
        non_ref_counts[ref_index] = 0;
        std::sort(non_ref_counts.begin(), non_ref_counts.end(), std::greater<int>());
        if (selected_features.count("x2")) results(i, col++) = non_ref_counts[0];
        if (selected_features.count("x3")) results(i, col++) = non_ref_counts[1];
        if (selected_features.count("x4")) results(i, col++) = non_ref_counts[2];
      }

      if (selected_features.count("x5")) {
        col++;
      }

      if (selected_features.count("x6")) {
        col++;
      }

      if (selected_features.count("x7")) {
        col++;
      }

      if (selected_features.count("x8")) {
        col++;
      }

      if (selected_features.count("x9")) {
        col++;
      }

      if (selected_features.count("x10")) {
        results(i, col++) = log_error_prob(counts, ref_index);
      }

      if (selected_features.count("x11")) {
        results(i, col++) = read_position_variance(counts);
      }

      if (selected_features.count("x12")) {
        results(i, col++) = read_position_mean(counts);
      }

      if (selected_features.count("x13")) {
        results(i, col++) = static_cast<int>(check_match_with_max_coverage(counts, ref_index));
      }

      if (selected_features.count("x14")) {
        results(i, col++) = static_cast<int>(check_match_with_non_ref_quantity(counts, ref_index));
      }
    }
  }
};

// [[Rcpp::export]]
DataFrame vectorise_data_cpp_parallel(DataFrame& df_d, CharacterVector features) {
  std::unordered_set<std::string> selected_features;
  for (int i = 0; i < features.size(); ++i) {
    selected_features.insert(as<std::string>(features[i]));
  }
  std::size_t n = df_d.nrows();
  std::size_t num_features = selected_features.size();

  NumericMatrix results(n, num_features);

  size_t col = 0;

  if (selected_features.count("x1")) {
    col++;
  }
  if (selected_features.count("x2")) {
    col++;
  }
  if (selected_features.count("x3")) {
    col++;
  }
  if (selected_features.count("x4")) {
    col++;
  }

  if (selected_features.count("x5")) {
    results(_, col++) = snp_entropy_cpp_parallel(df_d, false);
  }
  if (selected_features.count("x6")) {
    results(_, col++) = snp_entropy_cpp_parallel(df_d, true);
  }
  if (selected_features.count("x7")) {
    results(_, col++) = snp_binomial_cpp_parallel(df_d);
  }
  if (selected_features.count("x8")) {
    results(_, col++) = snp_fisher_cpp_parallel(df_d);
  }
  if (selected_features.count("x9")) {
    results(_, col++) = snp_poisson_cpp_parallel(df_d);
  }

  VectoriseWorker worker(df_d["A"], df_d["C"], df_d["G"], df_d["T"],
                         df_d["reference"], results, selected_features);
  parallelFor(0, n, worker);

  DataFrame output = DataFrame(results);

  CharacterVector col_names(features.size());
  for (R_xlen_t i = 0; i < features.size(); ++i) {
    col_names[i] = features[i];
  }
  output.attr("names") = col_names;

  return output;
}
