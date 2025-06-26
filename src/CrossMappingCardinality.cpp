#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <utility>
#include <unordered_set>
#include "CppStats.h"
#include "tEDMDataStruct.h"
#include "IntersectionCardinality.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

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
 * @param threads Number of threads for parallel processing.
 * @param parallel_level Level of parallelism to control nested parallel execution.
 * @param progressbar Boolean flag to show or hide a progress bar.
 *
 * @return CMCRes A struct containing:
 *   - cross_mapping: A vector of AUC values for the largest library size.
 *   - causal_strength: A 2D vector with rows [library size, mean AUC] across all lib sizes.
 */
CMCRes CrossMappingCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<size_t>& lib_sizes,
    const std::vector<size_t>& lib,
    const std::vector<size_t>& pred,
    size_t num_neighbors = 4,
    size_t n_excluded = 0,
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

  // Precompute neighbors
  auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, false, true), lib);
  auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, false, true), lib);

  // Local results for each library
  std::vector<std::vector<IntersectionRes>> local_results(unique_lib_sizes.size());

  if (parallel_level == 0){
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      for(size_t i = 0; i < unique_lib_sizes.size(); ++i){
        local_results[i] = IntersectionCardinalitySingle(
          nx,ny,unique_lib_sizes[i],lib,valid_pred,
          num_neighbors, n_excluded,
          threads_sizet, parallel_level
        );
        bar++;
      }
    } else {
      for(size_t i = 0; i < unique_lib_sizes.size(); ++i){
        local_results[i] = IntersectionCardinalitySingle(
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
        local_results[i] = IntersectionCardinalitySingle(
          nx,ny,unique_lib_sizes[i],lib,valid_pred,
          num_neighbors, n_excluded,
          threads_sizet, parallel_level
        );
        bar++;
      }, threads_sizet);
    } else {
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        local_results[i] = IntersectionCardinalitySingle(
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

// // Previous version of CCM implementation, retained for comparison
// #include <vector>
// #include <cmath>
// #include <algorithm>
// #include <numeric>
// #include <limits>
// #include <utility>
// #include <unordered_set>
// #include "CppStats.h"
// #include <RcppThread.h>
//
// // [[Rcpp::depends(RcppThread)]]
//
// /**
//  * Computes the Cross Mapping Cardinality (CMC) causal strength score (adjusted based on Python logic).
//  *
//  * Parameters:
//  *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
//  *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
//  *   lib: Library index vector (1-based in R, converted to 0-based).
//  *   pred: Prediction index vector (1-based in R, converted to 0-based).
//  *   num_neighbors: Vector of numbers of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
//  *   n_excluded: Vector of numbers of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
//  *   threads: Number of parallel threads.
//  *   parallel_level: the level of parallelization
//  *   progressbar: Whether to display a progress bar.
//  *
//  * Returns:
//  *   A vector the results of the DeLong test for the AUC values: [number of neighbors, IC score, p-value, confidence interval upper bound, confidence interval lower bound] one for each entry in num_neighbors.
//  *   The result contains multiple rows, each corresponding to a different number of neighbors.
//  */
// std::vector<std::vector<double>> CrossMappingCardinality(
//     const std::vector<std::vector<double>>& embedding_x,
//     const std::vector<std::vector<double>>& embedding_y,
//     const std::vector<int>& lib,
//     const std::vector<int>& pred,
//     const std::vector<int>& num_neighbors,
//     const std::vector<int>& n_excluded,
//     int threads,
//     int parallel_level = 0,
//     bool progressbar = true) {
//   // Store results for each num_neighbors
//   std::vector<std::vector<double>> results(num_neighbors.size(), std::vector<double>(5,std::numeric_limits<double>::quiet_NaN()));
//
//   // Input validation
//   if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
//     return results;
//   }
//
//   // Filter valid prediction points (exclude those with all NaN values)
//   std::vector<int> valid_pred;
//   for (int idx : pred) {
//     if (idx < 0 || static_cast<size_t>(idx) >= embedding_x.size()) continue;
//
//     bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
//                              [](double v) { return std::isnan(v); });
//     bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
//                              [](double v) { return std::isnan(v); });
//     if (!x_nan && !y_nan) valid_pred.push_back(idx);
//   }
//   if (valid_pred.empty()) return results;
//
//   // Configure threads
//   size_t threads_sizet = static_cast<size_t>(std::abs(threads));
//   threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);
//
//   // Precompute distance matrices (corresponding to _dismats in python package crossmapy)
//   auto dist_x = CppMatDistance(embedding_x, false, true);
//   auto dist_y = CppMatDistance(embedding_y, false, true);
//
//   if (parallel_level == 0){
//     auto CMCSingle = [&](size_t j) {
//       const size_t k = static_cast<size_t>(num_neighbors[j]);
//       const size_t n_excluded_sizet = static_cast<size_t>(n_excluded[j]);
//       const size_t max_r = k + n_excluded_sizet; // Total number of neighbors = actual used + excluded ones
//
//       // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
//       std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, std::numeric_limits<double>::quiet_NaN()));
//
//       // Perform the operations using RcppThread
//       RcppThread::parallelFor(0, valid_pred.size(), [&](size_t i) {
//         const int idx = valid_pred[i];
//
//         // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
//         auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r, lib);
//         if (neighbors_x.size() > n_excluded_sizet) {
//           neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded_sizet);
//         }
//         neighbors_x.resize(k); // Keep only the k actual neighbors
//
//         // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
//         auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r, lib);
//         if (neighbors_y.size() > n_excluded_sizet) {
//           neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded_sizet);
//         }
//         neighbors_y.resize(k); // Keep only the k actual neighbors
//
//         // Precompute y-neighbors set for fast lookup
//         std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());
//
//         // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
//         std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
//         for (size_t nx : neighbors_x) {
//           mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k, lib);
//         }
//
//         // Compute intersection ratio between mapped x-neighbors and original y-neighbors
//         for (size_t ki = 0; ki < k; ++ki) {
//           size_t count = 0;
//           for (size_t nx : neighbors_x) {
//             if (ki < mapped_neighbors[nx].size()) {
//               auto& yn = mapped_neighbors[nx];
//               // Check if any of first ki+1 mapped neighbors exist in y's original neighbors
//               for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
//                 if (y_neighbors_set.count(yn[pos])) {
//                   ++count;
//                   break; // Count each x-neighbor only once if has any intersection
//                 }
//               }
//             }
//           }
//           if (!neighbors_x.empty()) {
//             ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
//           }
//         }
//       }, threads_sizet);
//
//       std::vector<double> H1sequence;
//       for (size_t col = 0; col < k; ++col) {
//         std::vector<double> mean_intersect;
//         for (size_t row = 0; row < ratio_curves.size(); ++row){
//           mean_intersect.push_back(ratio_curves[row][col]);
//         }
//         H1sequence.push_back(CppMean(mean_intersect,true));
//       }
//
//       std::vector<double> dp_res = CppCMCTest(H1sequence,">");
//       dp_res.insert(dp_res.begin(), k);
//       results[j] = dp_res;
//     };
//
//     if (progressbar) {
//       RcppThread::ProgressBar bar(num_neighbors.size(), 1);
//       for (size_t i = 0; i < num_neighbors.size(); ++i){
//         CMCSingle(i);
//         bar++;
//       }
//     } else {
//       for (size_t i = 0; i < num_neighbors.size(); ++i){
//         CMCSingle(i);
//       }
//     }
//   } else {
//     auto CMCSingle = [&](size_t j) {
//       const size_t k = static_cast<size_t>(num_neighbors[j]);
//       const size_t n_excluded_sizet = static_cast<size_t>(n_excluded[j]);
//       const size_t max_r = k + n_excluded_sizet; // Total number of neighbors = actual used + excluded ones
//
//       // Store mapping ratio curves for each prediction point (corresponding to ratios_x2y in python package crossmapy)
//       std::vector<std::vector<double>> ratio_curves(valid_pred.size(), std::vector<double>(k, std::numeric_limits<double>::quiet_NaN()));
//
//       for (size_t i = 0; i < valid_pred.size(); ++i){
//         const int idx = valid_pred[i];
//
//         // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
//         auto neighbors_x = CppDistKNNIndice(dist_x, idx, max_r, lib);
//         if (neighbors_x.size() > n_excluded_sizet) {
//           neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded_sizet);
//         }
//         neighbors_x.resize(k); // Keep only the k actual neighbors
//
//         // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
//         auto neighbors_y = CppDistKNNIndice(dist_y, idx, max_r, lib);
//         if (neighbors_y.size() > n_excluded_sizet) {
//           neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded_sizet);
//         }
//         neighbors_y.resize(k); // Keep only the k actual neighbors
//
//         // Precompute y-neighbors set for fast lookup
//         std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());
//
//         // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
//         std::vector<std::vector<size_t>> mapped_neighbors(embedding_x.size());
//         for (size_t nx : neighbors_x) {
//           mapped_neighbors[nx] = CppDistKNNIndice(dist_y, nx, k, lib);
//         }
//
//         // Compute intersection ratio between mapped x-neighbors and original y-neighbors
//         for (size_t ki = 0; ki < k; ++ki) {
//           size_t count = 0;
//           for (size_t nx : neighbors_x) {
//             if (ki < mapped_neighbors[nx].size()) {
//               auto& yn = mapped_neighbors[nx];
//               // Check if any of first ki+1 mapped neighbors exist in y's original neighbors
//               for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
//                 if (y_neighbors_set.count(yn[pos])) {
//                   ++count;
//                   break; // Count each x-neighbor only once if has any intersection
//                 }
//               }
//             }
//           }
//           if (!neighbors_x.empty()) {
//             ratio_curves[i][ki] = static_cast<double>(count) / neighbors_x.size();
//           }
//         }
//       }
//
//       std::vector<double> H1sequence;
//       for (size_t col = 0; col < k; ++col) {
//         std::vector<double> mean_intersect;
//         for (size_t row = 0; row < ratio_curves.size(); ++row){
//           mean_intersect.push_back(ratio_curves[row][col]);
//         }
//         H1sequence.push_back(CppMean(mean_intersect,true));
//       }
//
//       std::vector<double> dp_res = CppCMCTest(H1sequence,">");
//       dp_res.insert(dp_res.begin(), k);
//       results[j] = dp_res;
//     };
//
//     if (progressbar) {
//       RcppThread::ProgressBar bar(num_neighbors.size(), 1);
//       RcppThread::parallelFor(0, num_neighbors.size(), [&](size_t i) {
//         CMCSingle(i);
//         bar++;
//       }, threads_sizet);
//     } else {
//       RcppThread::parallelFor(0, num_neighbors.size(), CMCSingle, threads_sizet);
//     }
//   }
//
//   return results; // Return the vector of results
// }
