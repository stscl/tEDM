#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include "CppStats.h"
#include "CppDistances.h"
#include "DataStruct.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/**
 * Computes intersection-based mapping ratio sequences between two neighbor graphs
 * for use in Cross Mapping Cardinality (CMC) or similar causal inference frameworks.
 *
 * Parameters:
 *   neighborsX     - Precomputed sorted neighbor indices for embedding X
 *   neighborsY     - Precomputed sorted neighbor indices for embedding Y
 *   lib_size       - Size of the moving library used in mapping
 *   lib_indices    - Global indices from which to draw the sliding libraries
 *   pred_indices   - Indices at which to perform prediction (evaluation points)
 *   num_neighbors  - Number of neighbors used for mapping (after exclusion)
 *   n_excluded     - Number of nearest neighbors to exclude from the front
 *   threads        - Number of parallel threads for computation
 *   parallel_level - Whether to use multithreaded (0) or serial (1) mode
 *
 * Returns:
 *   A vector of IntersectionRes structures, each containing the average intersection
 *   ratio sequence (IC curve) for a different starting point of the moving library.
 *   If lib_size == lib_indices.size(), returns a single result using full library.
 *
 * Notes:
 *   - Neighbor lists must use std::numeric_limits<size_t>::max() to indicate invalid entries.
 *   - This function assumes that the neighbor vectors are sorted by ascending distance.
 *   - Use in combination with AUC computation to assess causal strength.
 */
std::vector<IntersectionRes> IntersectionCardinalitySingle(
    const std::vector<std::vector<size_t>>& neighborsX,
    const std::vector<std::vector<size_t>>& neighborsY,
    size_t lib_size,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    size_t num_neighbors,
    size_t n_excluded,
    size_t threads,
    int parallel_level
){
  size_t max_lib_size = lib_indices.size();
  const size_t max_r = num_neighbors + n_excluded; // Total number of neighbors = actual used + excluded ones

  auto ICSingle = [&](const std::vector<size_t>& lib) {
    // Store mapping ratio curves for each prediction point
    std::vector<std::vector<double>> ratio_curves(pred_indices.size(), std::vector<double>(num_neighbors, std::numeric_limits<double>::quiet_NaN()));

    // Precompute library set for fast lookup
    std::unordered_set<size_t> lib_set(lib.begin(), lib.end());

    if (parallel_level == 0){
      // Perform the operations using RcppThread
      RcppThread::parallelFor(0, pred_indices.size(), [&](size_t i) {
        const size_t idx = pred_indices[i];

        if ((neighborsX[idx][0] != std::numeric_limits<size_t>::max()) &&
            (neighborsY[idx][0] != std::numeric_limits<size_t>::max())){
          if ((neighborsX[idx].size() >= max_r) && (neighborsY[idx].size() >= max_r)){
            // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
            std::vector<size_t> neighbors_x;
            for (size_t iidx = 0; iidx < neighborsX[idx].size(); ++iidx){
              if(lib_set.count(neighborsX[idx][iidx]) > 0){
                neighbors_x.push_back(neighborsX[idx][iidx]);
                if (neighbors_x.size() >= max_r) break;
              }
            }
            if (neighbors_x.size() > n_excluded) {
              neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded);
            }

            // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
            std::vector<size_t> neighbors_y;
            for (size_t iidx = 0; iidx < neighborsY[idx].size(); ++iidx){
              if(lib_set.count(neighborsY[idx][iidx])){
                neighbors_y.push_back(neighborsY[idx][iidx]);
                if (neighbors_y.size() >= max_r) break;
              }
            }
            if (neighbors_y.size() > n_excluded) {
              neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded);
            }

            // Precompute y-neighbors set for fast lookup
            std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());

            // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
            std::unordered_map<size_t, std::vector<size_t>> mapped_neighbors;

            for (size_t nx : neighbors_x) {
              if (neighborsY[nx][0] != std::numeric_limits<size_t>::max()) {
                for (size_t iidx = 0; iidx < neighborsY[nx].size(); ++iidx) {
                  if (lib_set.count(neighborsY[nx][iidx])) {
                    mapped_neighbors[nx].push_back(neighborsY[nx][iidx]);
                    if (mapped_neighbors[nx].size() >= num_neighbors) break;
                  }
                }
              }
            }

            // Compute intersection ratio between mapped x-neighbors and original y-neighbors
            for (size_t ki = 0; ki < num_neighbors; ++ki) {
              size_t count = 0;

              for (size_t nx : neighbors_x) {
                auto it = mapped_neighbors.find(nx);
                if (it != mapped_neighbors.end() && ki < it->second.size()) {
                  const auto& yn = it->second;

                  // Check if any of the first (ki+1) mapped neighbors exist in y's original neighbors
                  for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
                    if (y_neighbors_set.count(yn[pos])) {
                      ++count;
                      break;  // Count each x-neighbor only once if it intersects
                    }
                  }
                }
              }

              if (neighbors_x.size() > 0) {
                ratio_curves[i][ki] = static_cast<double>(count) / static_cast<double>(neighbors_x.size());
              }
            }
          }
        }
      }, threads);
    } else {
      // Perform the operations one by one
      for (size_t i = 0; i < pred_indices.size(); ++i){
        const size_t idx = pred_indices[i];

        if ((neighborsX[idx][0] != std::numeric_limits<size_t>::max()) &&
            (neighborsY[idx][0] != std::numeric_limits<size_t>::max())){
          if ((neighborsX[idx].size() >= max_r) && (neighborsY[idx].size() >= max_r)){
            // Get the k-nearest neighbors of x (excluding the first n_excluded ones)
            std::vector<size_t> neighbors_x;
            for (size_t iidx = 0; iidx < neighborsX[idx].size(); ++iidx){
              if(lib_set.count(neighborsX[idx][iidx]) > 0){
                neighbors_x.push_back(neighborsX[idx][iidx]);
                if (neighbors_x.size() >= max_r) break;
              }
            }
            if (neighbors_x.size() > n_excluded) {
              neighbors_x.erase(neighbors_x.begin(), neighbors_x.begin() + n_excluded);
            }

            // Get the k-nearest neighbors of y (excluding the first n_excluded ones)
            std::vector<size_t> neighbors_y;
            for (size_t iidx = 0; iidx < neighborsY[idx].size(); ++iidx){
              if(lib_set.count(neighborsY[idx][iidx])){
                neighbors_y.push_back(neighborsY[idx][iidx]);
                if (neighbors_y.size() >= max_r) break;
              }
            }
            if (neighbors_y.size() > n_excluded) {
              neighbors_y.erase(neighbors_y.begin(), neighbors_y.begin() + n_excluded);
            }

            // Precompute y-neighbors set for fast lookup
            std::unordered_set<size_t> y_neighbors_set(neighbors_y.begin(), neighbors_y.end());

            // Retrieve y's neighbor indices by mapping x-neighbors through x->y mapping
            std::unordered_map<size_t, std::vector<size_t>> mapped_neighbors;

            for (size_t nx : neighbors_x) {
              if (neighborsY[nx][0] != std::numeric_limits<size_t>::max()) {
                for (size_t iidx = 0; iidx < neighborsY[nx].size(); ++iidx) {
                  if (lib_set.count(neighborsY[nx][iidx])) {
                    mapped_neighbors[nx].push_back(neighborsY[nx][iidx]);
                    if (mapped_neighbors[nx].size() >= num_neighbors) break;
                  }
                }
              }
            }

            // Compute intersection ratio between mapped x-neighbors and original y-neighbors
            for (size_t ki = 0; ki < num_neighbors; ++ki) {
              size_t count = 0;

              for (size_t nx : neighbors_x) {
                auto it = mapped_neighbors.find(nx);
                if (it != mapped_neighbors.end() && ki < it->second.size()) {
                  const auto& yn = it->second;

                  // Check if any of the first (ki+1) mapped neighbors exist in y's original neighbors
                  for (size_t pos = 0; pos <= ki && pos < yn.size(); ++pos) {
                    if (y_neighbors_set.count(yn[pos])) {
                      ++count;
                      break;  // Count each x-neighbor only once if it intersects
                    }
                  }
                }
              }

              if (neighbors_x.size() > 0) {
                ratio_curves[i][ki] = static_cast<double>(count) / static_cast<double>(neighbors_x.size());
              }
            }
          }
        }
      }
    }

    std::vector<double> H1sequence;
    for (size_t col = 0; col < num_neighbors; ++col) {
      std::vector<double> mean_intersect;
      for (size_t row = 0; row < ratio_curves.size(); ++row){
        mean_intersect.push_back(ratio_curves[row][col]);
      }
      H1sequence.push_back(CppMean(mean_intersect,true));
    }

    return H1sequence;
  };

  // No possible library variation if using all vectors
  if (lib_size == max_lib_size) {
    std::vector<IntersectionRes> x_xmap_y;
    x_xmap_y.emplace_back(lib_size, ICSingle(lib_indices));
    return x_xmap_y;
  } else {
    if (parallel_level == 0){
      parallel_level += 1; // change parallelism level
      // Precompute valid indices for the library
      std::vector<std::vector<size_t>> valid_lib_indices;
      for (size_t start_lib = 0; start_lib < max_lib_size; ++start_lib) {
        std::vector<size_t> local_lib_indices;
        // Loop around to beginning of lib indices
        if (start_lib + lib_size > max_lib_size) {
          for (size_t i = start_lib; i < max_lib_size; ++i) {
            local_lib_indices.emplace_back(lib_indices[i]);
          }
          int num_vectors_remaining = static_cast<int>(lib_size) - (static_cast<int>(max_lib_size) - static_cast<int>(start_lib));
          for (int i = 0; i < num_vectors_remaining; ++i) {
            local_lib_indices.emplace_back(lib_indices[i]);
          }
        } else {
          for (size_t i = start_lib; i < start_lib + lib_size; ++i) {
            local_lib_indices.emplace_back(lib_indices[i]);
          }
        }
        valid_lib_indices.emplace_back(local_lib_indices);
      }

      // Preallocate the result vector to avoid out-of-bounds access
      std::vector<IntersectionRes> x_xmap_y(valid_lib_indices.size());

      // Perform the operations using RcppThread
      RcppThread::parallelFor(0, valid_lib_indices.size(), [&](size_t i) {
        IntersectionRes result(lib_size, ICSingle(valid_lib_indices[i]));
        x_xmap_y[i] = result;
      }, threads);
      return x_xmap_y;
    } else {
      // Serial version
      parallel_level += 1; // change parallelism level
      // Preallocate the result vector to avoid out-of-bounds access
      std::vector<IntersectionRes> x_xmap_y;

      for (size_t start_lib = 0; start_lib < max_lib_size; ++start_lib) {
        std::vector<size_t> local_lib_indices;
        // Setup changing library
        if (start_lib + lib_size > max_lib_size) { // Loop around to beginning of lib indices
          for (size_t i = start_lib; i < max_lib_size; ++i) {
            local_lib_indices.emplace_back(lib_indices[i]);
          }
          int num_vectors_remaining = static_cast<int>(lib_size) - (static_cast<int>(max_lib_size) - static_cast<int>(start_lib));
          for (int i = 0; i < num_vectors_remaining; ++i) {
            local_lib_indices.emplace_back(lib_indices[i]);
          }
        } else {
          for (size_t i = start_lib; i < start_lib + lib_size; ++i) {
            local_lib_indices.emplace_back(lib_indices[i]);
          }
        }

        // Run cross map and store results
        x_xmap_y.emplace_back(lib_size, ICSingle(local_lib_indices));
      }

      return x_xmap_y;
    }
  }
}

/**
 * Computes the Intersection Cardinality (IC) curve for causal inference via cross mapping.
 *
 * This function evaluates the extent to which neighbors of the effect variable Y
 * are preserved when mapped through the neighbors of cause variable X.
 * Specifically, for each number of neighbors from 1 to `num_neighbors`, it computes
 * the intersection count between the k nearest neighbors of Y and the k nearest neighbors of X,
 * for each prediction point.
 *
 * The output is an Intersection Cardinality (IC) curve, which can be further processed
 * (e.g., calculating AUC, statistical significance) outside this function.
 *
 * @param embedding_x     State-space embedding of the potential cause variable (NxE matrix).
 * @param embedding_y     State-space embedding of the potential effect variable (NxE matrix).
 * @param lib             Vector of library indices (shouble be 0-based in C++).
 * @param pred            Vector of prediction indices (shouble be 0-based in C++).
 * @param num_neighbors   Maximum number of neighbors to consider in intersection (e.g., from 1 to k).
 * @param n_excluded      Number of nearest neighbors to exclude (e.g., due to temporal proximity).
 * @param dist_metric     Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param threads         Number of threads used for parallel computation.
 * @param parallel_level  Parallel mode flag: 0 = parallel, 1 = serial.
 *
 * @return A vector of size `num_neighbors`:
 *         - Each element represents the average number of overlapping neighbors
 *           between X and Y across prediction points, for each neighbor count k = 1, 2, ..., num_neighbors.
 *
 *         If inputs are invalid or no valid prediction points exist, the returned vector
 *         is filled with `NaN` values.
 *
 * @note
 *   - This function returns only the raw intersection values. To compute AUC or p-values,
 *     use additional post-processing such as DeLong’s test.
 */
std::vector<double> IntersectionCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<size_t>& lib,
    const std::vector<size_t>& pred,
    size_t num_neighbors,
    size_t n_excluded,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0) {
  std::vector<double> result (num_neighbors, std::numeric_limits<double>::quiet_NaN());
  // Input validation
  if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
    return result;
  }

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
  if (valid_pred.empty()) return result;

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
  bool L1norm = (dist_metric == 1);

  // // Precompute neighbors (The earlier implementation based on a serial version)
  // auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, L1norm, true), lib, num_neighbors + n_excluded);
  // auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, L1norm, true), lib, num_neighbors + n_excluded);

  // Precompute neighbors (parallel computation)
  auto nx = CppMatKNNeighbors(embedding_x, lib, num_neighbors + n_excluded, threads_sizet, L1norm);
  auto ny = CppMatKNNeighbors(embedding_y, lib, num_neighbors + n_excluded, threads_sizet, L1norm);

  // run cross mapping
  std::vector<IntersectionRes> res = IntersectionCardinalitySingle(
    nx,ny,lib.size(),lib,pred,num_neighbors,n_excluded,threads_sizet,parallel_level
  );

  if (res.empty()) return result;

  return res[0].Intersection;
}

/**
 * Computes the Intersection Cardinality (IC) AUC-based causal strength score.
 *
 * This function evaluates the extent to which neighbors of the effect variable Y
 * are preserved when mapped through the neighbors of cause variable X, by calculating
 * the intersection ratio curve and evaluating its AUC (area under curve).
 * Statistical significance (p-value) and confidence interval are computed via DeLong's test.
 *
 * Parameters:
 *   embedding_x    - State-space reconstruction (embedding) of the potential cause variable.
 *   embedding_y    - State-space reconstruction (embedding) of the potential effect variable.
 *   lib            - Library index vector (shouble be 0-based in C++).
 *   pred           - Prediction index vector (shouble be 0-based in C++).
 *   num_neighbors  - Number of neighbors used for cross mapping (after exclusion).
 *   n_excluded     - Number of nearest neighbors to exclude (e.g. temporal).
 *   dist_metric    - Distance metric selector (1: Manhattan, 2: Euclidean).
 *   threads        - Number of threads used in parallel computation.
 *   parallel_level - Whether to use multithreaded (0) or serial (1) mode
 *
 * Returns:
 *   A vector of 4 values:
 *     [0] - AUC (Intersection Cardinality score, bounded [0, 1])
 *     [1] - p-value from DeLong test (testing whether AUC > 0.5)
 *     [2] - Confidence interval lower bound
 *     [3] - Confidence interval upper bound
 */
std::vector<double> IntersectionCardinalityScores(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<size_t>& lib,
    const std::vector<size_t>& pred,
    size_t num_neighbors,
    size_t n_excluded,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0) {
  // Input validation
  if (embedding_x.size() != embedding_y.size() || embedding_x.empty()) {
    return {0, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  }

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
  if (valid_pred.empty()) return {0, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
  bool L1norm = (dist_metric == 1);

  // // Precompute neighbors (The earlier implementation based on a serial version)
  // auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, L1norm, true), lib, num_neighbors + n_excluded);
  // auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, L1norm, true), lib, num_neighbors + n_excluded);

  // Precompute neighbors (parallel computation)
  auto nx = CppMatKNNeighbors(embedding_x, lib, num_neighbors + n_excluded, threads_sizet, L1norm);
  auto ny = CppMatKNNeighbors(embedding_y, lib, num_neighbors + n_excluded, threads_sizet, L1norm);

  // run cross mapping
  std::vector<IntersectionRes> res = IntersectionCardinalitySingle(
    nx,ny,lib.size(),lib,pred,num_neighbors,n_excluded,threads_sizet,parallel_level
  );

  if (res.empty()) return {0, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

  return CppCMCTest(res[0].Intersection,">");
}
