#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <utility>
#include <map>
#include <unordered_set>
#include "CppStats.h"
#include "CppDistances.h"
#include "DataStruct.h"
#include "IntersectionalCardinality.h"
#include <RcppThread.h>

/**
 * @brief Computes the Cross Mapping Cardinality (CMC) causal strength score.
 *
 * This function evaluates the directional causal influence from one time series
 * to another using Cross Mapping Cardinality. It performs state-space reconstruction,
 * neighbor searching, and statistical evaluation across a range of library sizes.
 *
 * @param embedding_x State-space reconstructed time series of the potential cause.
 * @param embedding_y State-space reconstructed time series of the potential effect.
 * @param lib_sizes A vector of library sizes to use for subsampling during CMC analysis.
 * @param lib Indices of library points (0-based).
 * @param pred Indices of prediction points (0-based).
 * @param num_neighbors Number of neighbors used in cross mapping.
 * @param n_excluded Number of temporally excluded neighbors (Theiler window).
 * @param dist_metric Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param threads Number of threads for parallel processing.
 * @param parallel_level Level of parallelism to control nested parallel execution.
 * @param progressbar Boolean flag to show or hide a progress bar.
 *
 * @return CMCRes A struct containing:
 *   - cross_mapping: A vector of AUC values for the largest library size.
 *   - causal_strength: A 2D vector with rows [library size, mean AUC] across all lib sizes.
 */
CMCRes CMC(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<size_t>& lib_sizes,
    const std::vector<size_t>& lib,
    const std::vector<size_t>& pred,
    size_t num_neighbors = 4,
    size_t n_excluded = 0,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0,
    bool progressbar = true) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Filter valid prediction points (exclude those with all NaN values)
  std::vector<size_t> valid_pred;
  for (size_t idx : pred) {
    if (idx < 0 || idx >= embedding_x.size()) continue;

    bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                             [](double v) { return std::isnan(v); });
    bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                             [](double v) { return std::isnan(v); });
    if (!x_nan && !y_nan) valid_pred.push_back(idx);
  }

  // Transform to ensure no size exceeds max library number
  size_t max_lib_size = lib.size();
  std::vector<size_t> unique_lib_sizes(lib_sizes.begin(), lib_sizes.end());
  unique_lib_sizes.push_back(max_lib_size);  // Ensure max_lib_size is included

  // // Clamp each size between (num_neighbors + n_excluded) and max_lib_size (need C++17, so commented)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](size_t size) {
  //                  return std::clamp(size, num_neighbors + n_excluded, max_lib_size);
  //                });

  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](size_t size) { return std::min(size, max_lib_size); });

  // Ensure the minimum value in unique_lib_sizes is num_neighbors + n_excluded
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](size_t size) { return std::max(size, num_neighbors + n_excluded ); });

  // Remove duplicates
  std::sort(unique_lib_sizes.begin(), unique_lib_sizes.end());
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
  bool L1norm = (dist_metric == 1);

  // // Precompute neighbors (The earlier implementation based on a serial version)
  // auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, L1norm, true), lib, num_neighbors + n_excluded);
  // auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, L1norm, true), lib, num_neighbors + n_excluded);

  // Precompute neighbors (parallel computation)
  auto nx = CppMatKNNeighbors(embedding_x, lib, num_neighbors + n_excluded, threads_sizet, L1norm);
  auto ny = CppMatKNNeighbors(embedding_y, lib, num_neighbors + n_excluded, threads_sizet, L1norm);

  // Local results for each library
  std::vector<std::vector<IntersectionRes>> local_results(unique_lib_sizes.size());

  if (parallel_level == 0){
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      for(size_t i = 0; i < unique_lib_sizes.size(); ++i){
        local_results[i] = IntersectionalCardinalitySingle(
          nx,ny,unique_lib_sizes[i],lib,valid_pred,
          num_neighbors, n_excluded,
          threads_sizet, parallel_level
        );
        bar++;
      }
    } else {
      for(size_t i = 0; i < unique_lib_sizes.size(); ++i){
        local_results[i] = IntersectionalCardinalitySingle(
          nx,ny,unique_lib_sizes[i],lib,valid_pred,
          num_neighbors, n_excluded,
          threads_sizet, parallel_level
        );
      }
    }
  } else {
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        local_results[i] = IntersectionalCardinalitySingle(
          nx,ny,unique_lib_sizes[i],lib,valid_pred,
          num_neighbors, n_excluded,
          threads_sizet, parallel_level
        );
        bar++;
      }, threads_sizet);
    } else {
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        local_results[i] = IntersectionalCardinalitySingle(
          nx,ny,unique_lib_sizes[i],lib,valid_pred,
          num_neighbors, n_excluded,
          threads_sizet, parallel_level
        );
      }, threads_sizet);
    }
  }

  // Output vector to store pairs of (libsize, AUC)
  std::vector<std::pair<int, double>> x_xmap_y;

  // Merge all local results into one vector
  std::vector<IntersectionRes> H1_vector;
  for (const auto& local_result : local_results) {
    H1_vector.insert(H1_vector.end(), local_result.begin(), local_result.end());
  }

  for (const auto& h1 : H1_vector) {
    // Compute AUC metrics from the intersection values using CppCMCTest
    std::vector<double> auc_result = CppCMCTest(h1.Intersection, ">");

    // If the result is non-empty, store the first AUC value along with libsize
    if (!auc_result.empty()) {
      x_xmap_y.emplace_back(static_cast<int>(h1.libsize), auc_result[0]);
    }
  }

  // Group by the first int(libsize) and compute the mean (auc)
  std::map<int, std::vector<double>> grouped_results;
  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].push_back(result.second);
  }

  std::vector<std::vector<double>> mean_aucs;
  for (const auto& group : grouped_results) {
    double mean_value = CppMean(group.second, true);
    mean_aucs.push_back({static_cast<double>(group.first), mean_value});
  }

  // Find the largest valid libsize
  size_t largest_libsize = unique_lib_sizes.back();

  // Locate the corresponding PartialCorRes from H1_vector
  std::vector<double> result_auc;
  for (const auto& h1 : H1_vector) {
    if (h1.libsize == largest_libsize) {
      // Run full CppCMCTest on this intersection vector
      result_auc = CppCMCTest(h1.Intersection, ">");
      result_auc.insert(result_auc.begin(), static_cast<double>(num_neighbors));
      break;
    }
  }

  CMCRes result;
  result.cross_mapping = result_auc;
  result.causal_strength = mean_aucs;

  return result;
}
