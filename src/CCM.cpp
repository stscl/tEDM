#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include <map>
#include "CppStats.h"
#include "Embed.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/*
 * Perform convergent cross mapping on a single lib and pred.
 *
 * Parameters:
 *   - x_vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - y: time series used as the target (should align with x_vectors).
 *   - lib_size: Size of the library used for cross mapping.
 *   - lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of indices indicating which states to predict from.
 *   - b: Number of neighbors to use for simplex projection.
 *   - simplex: If true, uses simplex projection for prediction; otherwise, uses s-mapping.
 *   - theta: Distance weighting parameter for local neighbors in the manifold (used in s-mapping).
 *   - threads: The number of threads to use for parallel processing.
 *   - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components.
 *
 * Returns:
 *   A vector of pairs, where each pair consists of:
 *   - An integer representing the library size.
 *   - A double representing the Pearson correlation coefficient (rho) between predicted and actual values.
 */
std::vector<std::pair<int, double>> CCMSingle(
    const std::vector<std::vector<double>>& x_vectors,
    const std::vector<double>& y,
    int lib_size,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int b,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    int dist_metric,
    bool dist_average
) {
  int max_lib_size = lib_indices.size();

  // No possible library variation if using all vectors
  if (lib_size == max_lib_size) {
    std::vector<std::pair<int, double>> x_xmap_y;

    // Run cross map and store results
    double rho = std::numeric_limits<double>::quiet_NaN();
    if (simplex) {
      rho = SimplexProjection(x_vectors, y, lib_indices, pred_indices, b, dist_metric, dist_average);
    } else {
      rho = SMap(x_vectors, y, lib_indices, pred_indices, b, theta, dist_metric, dist_average);
    }
    x_xmap_y.emplace_back(lib_size, rho);
    return x_xmap_y;
  } else if (parallel_level == 0){
    // Precompute valid indices for the library
    std::vector<std::vector<int>> valid_lib_indices;
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Loop around to beginning of lib indices
      if (start_lib + lib_size > max_lib_size) {
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        // // no wrapping around
        // int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        // for (int i = 0; i < num_vectors_remaining; ++i) {
        //   local_lib_indices.emplace_back(lib_indices[i]);
        // }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }
      valid_lib_indices.emplace_back(local_lib_indices);
    }

    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<std::pair<int, double>> x_xmap_y(valid_lib_indices.size());

    // Perform the operations using RcppThread
    RcppThread::parallelFor(0, valid_lib_indices.size(), [&](size_t i) {
      // Run cross map and store results
      double rho = std::numeric_limits<double>::quiet_NaN();
      if (simplex) {
        rho = SimplexProjection(x_vectors, y, valid_lib_indices[i], pred_indices, b, dist_metric, dist_average);
      } else {
        rho = SMap(x_vectors, y, valid_lib_indices[i], pred_indices, b, theta, dist_metric, dist_average);
      }

      std::pair<int, double> result(lib_size, rho); // Store the product of row and column library sizes
      x_xmap_y[i] = result;
    }, threads);

    return x_xmap_y;
  } else {
    std::vector<std::pair<int, double>> x_xmap_y;

    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Setup changing library
      if (start_lib + lib_size > max_lib_size) { // Loop around to beginning of lib indices
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        // // no wrapping around
        // int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        // for (int i = 0; i < num_vectors_remaining; ++i) {
        //   local_lib_indices.emplace_back(lib_indices[i]);
        // }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }

      // Run cross map and store results
      double rho = std::numeric_limits<double>::quiet_NaN();
      if (simplex) {
        rho = SimplexProjection(x_vectors, y, local_lib_indices, pred_indices, b, dist_metric, dist_average);
      } else {
        rho = SMap(x_vectors, y, local_lib_indices, pred_indices, b, theta, dist_metric, dist_average);
      }
      x_xmap_y.emplace_back(lib_size, rho);
    }

    return x_xmap_y;
  }
}

/**
 * Performs convergent cross mapping on a time series data.
 *
 * Parameters:
 * - x: time series used as the predict variable (**cross mapping from**).
 * - y: time series used as the target variable (**cross mapping to**).
 * - lib_sizes: A vector specifying different library sizes for CCM analysis.
 * - lib: A vector of representing the indices of time points to be the library.
 * - pred: A vector of representing the indices of time points to be predicted.
 * - E: Embedding dimension for attractor reconstruction.
 * - tau: the step of stime lags for prediction.
 * - b: Number of nearest neighbors used for prediction.
 * - simplex: Boolean flag indicating whether to use simplex projection (true) or S-mapping (false) for prediction.
 * - theta: Distance weighting parameter used for weighting neighbors in the S-mapping prediction.
 * - threads: Number of threads to use for parallel computation.
 * - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 *   dist_metric    - Distance metric selector (1: Manhattan, 2: Euclidean).
 *   dist_average   - Whether to average distance by the number of valid vector components.
 * - progressbar: Boolean flag to indicate whether to display a progress bar during computation.
 *
 * Returns:
 *    A 2D vector of results, where each row contains:
 *      - The library size.
 *      - The mean cross-mapping correlation.
 *      - The statistical significance of the correlation.
 *      - The upper bound of the confidence interval.
 *      - The lower bound of the confidence interval.
 */
std::vector<std::vector<double>> CCM(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<int>& lib_sizes,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    int E = 3,
    int tau = 1,
    int b = 4,
    bool simplex = true,
    double theta = 1.0,
    int threads = 8,
    int parallel_level = 0,
    int dist_metric = 2,
    bool dist_average = true,
    bool progressbar = false
) {
  // If b is not provided correctly, default it to E + 1
  if (b <= 0) {
    b = E + 1;
  }

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Generate embeddings
  std::vector<std::vector<double>> x_vectors = Embed(x, E, tau);

  size_t n = pred.size();

  std::vector<int> unique_lib_sizes(lib_sizes.begin(), lib_sizes.end());

  // Transform to ensure no size exceeds max library size
  int max_lib_size = lib.size();

  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size); });

  // Ensure the minimum value in unique_lib_sizes is E + 1 (uncomment this section if required)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](int size) { return std::max(size, E + 1); });

  // Remove duplicates
  std::sort(unique_lib_sizes.begin(), unique_lib_sizes.end());
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Local results for each library
  std::vector<std::vector<std::pair<int, double>>> local_results(unique_lib_sizes.size());

  if (parallel_level == 0){
    // Iterate over each library size
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        local_results[i] = CCMSingle(
          x_vectors,
          y,
          unique_lib_sizes[i],
          lib,
          pred,
          b,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          dist_metric,
          dist_average);
        bar++;
      }
    } else {
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        local_results[i] = CCMSingle(
          x_vectors,
          y,
          unique_lib_sizes[i],
          lib,
          pred,
          b,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          dist_metric,
          dist_average);
      }
    }
  } else {
    // Perform the operations using RcppThread
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        local_results[i] = CCMSingle(
          x_vectors,
          y,
          lib_size,
          lib,
          pred,
          b,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          dist_metric,
          dist_average);
        bar++;
      }, threads_sizet);
    } else {
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        local_results[i] = CCMSingle(
          x_vectors,
          y,
          lib_size,
          lib,
          pred,
          b,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          dist_metric,
          dist_average);
      }, threads_sizet);
    }
  }

  // Initialize the result container
  std::vector<std::pair<int, double>> x_xmap_y;

  // Merge all local results into the final result
  for (const auto& local_result : local_results) {
    x_xmap_y.insert(x_xmap_y.end(), local_result.begin(), local_result.end());
  }

  // Group by the first int(lib_size) and compute the mean (rho)
  std::map<int, std::vector<double>> grouped_results;
  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].push_back(result.second);
  }

  // // Previous implementation calculated significance and confidence intervals using the mean of rho vector only.
  // // This approach is now deprecated and kept here for comparison purposes.
  // std::vector<std::vector<double>> final_results;
  // for (const auto& group : grouped_results) {
  //   double mean_value = CppMean(group.second, true);
  //   final_results.push_back({static_cast<double>(group.first), mean_value});
  // }
  //
  // // Calculate significance and confidence interval for each result
  // for (size_t i = 0; i < final_results.size(); ++i) {
  //   double rho = final_results[i][1];
  //   double significance = CppCorSignificance(rho, n);
  //   std::vector<double> confidence_interval = CppCorConfidence(rho, n);
  //
  //   final_results[i].push_back(significance);
  //   final_results[i].push_back(confidence_interval[0]);
  //   final_results[i].push_back(confidence_interval[1]);
  // }

  // Refactor correlation analysis to compute significance and confidence intervals directly from grouped correlation vectors
  std::vector<std::vector<double>> final_results;
  for (const auto& group : grouped_results) {
    // Calculate the mean correlation coefficient from the group
    double mean_value = CppMean(group.second, true);

    // Compute significance (p-value) using the vector of correlations directly
    double significance = CppMeanCorSignificance(group.second, n);

    // Compute confidence interval using the vector of correlations directly
    std::vector<double> confidence_interval = CppMeanCorConfidence(group.second, n);

    // Store results: group ID, mean correlation, p-value, lower CI, upper CI
    final_results.push_back({
      static_cast<double>(group.first),
      mean_value,
      significance,
      confidence_interval[0],
      confidence_interval[1]
    });
  }

  return final_results;
}
