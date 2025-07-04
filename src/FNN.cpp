#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "CppStats.h"
#include "CppDistances.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/*
 * Compute the False Nearest Neighbors (FNN) ratio for time serires data.
 *
 * This function determines whether nearest neighbors identified in a lower-dimensional
 * embedded space (E1) remain close in a higher-dimensional space (E2).
 * If not, the neighbor is considered a "false" neighbor, indicating the need for
 * a higher embedding dimension to accurately capture spatial proximity.
 *
 * The FNN test is computed in two modes:
 * - parallel_level = 0: each prediction is processed in parallel using RcppThreads.
 * - parallel_level = 1: all pairwise distances are precomputed once in advance
 *   (better for repeated queries or small prediction sets).
 *
 * Parameters:
 * - embedding: A matrix (vector of vectors) representing the spatial embedding,
 *              where each row corresponds to a spatial unit's attributes.
 *              Must contain at least E2 columns.
 * - lib: Library index vector (1-based in R, converted to 0-based).
 * - pred: Prediction index vector (1-based in R, converted to 0-based).
 * - E1: The base embedding dimension used to identify the nearest neighbor (E1 < E2).
 * - E2: The full embedding dimension used to test false neighbors (usually E1 + 1).
 * - threads: Number of threads used when parallel_level = 0.
 * - parallel_level: 0 for per-pred parallelism (default), 1 for precomputed full distance matrix.
 * - Rtol: Relative threshold (default 10.0).
 * - Atol: Absolute threshold (default 2.0).
 * - L1norm: Whether to use Manhattan (L1) distance instead of Euclidean (L2).
 *
 * Returns:
 * - A double value indicating the proportion of false nearest neighbors (0–1).
 *   If no valid pairs are found, returns NaN.
 */
double CppSingleFNN(const std::vector<std::vector<double>>& embedding,
                    const std::vector<int>& lib,
                    const std::vector<int>& pred,
                    size_t E1,
                    size_t E2,
                    size_t threads,
                    int parallel_level = 0,
                    double Rtol = 10.0,
                    double Atol = 2.0,
                    bool L1norm = false) {
  if (embedding.empty() || embedding[0].size() < E2) {
    return std::numeric_limits<double>::quiet_NaN();  // Invalid dimensions
  }

  size_t N = embedding.size();

  if (parallel_level >= 1){
    std::vector<std::vector<double>> distmat(N, std::vector<double>(N, std::numeric_limits<double>::quiet_NaN()));

    std::vector<int> merged = lib;
    merged.insert(merged.end(), pred.begin(), pred.end());
    std::sort(merged.begin(), merged.end());
    merged.erase(std::unique(merged.begin(), merged.end()), merged.end());

    // Compute distance between every pair of merged
    for (size_t i = 0; i < merged.size(); ++i) {
      for (size_t j = i+1; j < merged.size(); ++j) {
        std::vector<double> xi_E1(embedding[merged[i]].begin(), embedding[merged[i]].begin() + E1);
        std::vector<double> xj_E1(embedding[merged[j]].begin(), embedding[merged[j]].begin() + E1);
        double distv = CppDistance(xi_E1, xj_E1, L1norm, true);
        distmat[i][j] = distv;  // Correctly assign distance to upper triangle
        distmat[j][i] = distv;  // Mirror the value to the lower triangle
        // distmat[i][j] = distmat[j][i] = CppDistance(xi_E1, xj_E1, L1norm, true);
      }
    }

    int false_count = 0;
    int total = 0;

    // // Brute-force linear search leads to slow performance
    // for (size_t i = 0; i < pred.size(); ++i) {
    //   if (checkOneDimVectorNotNanNum(embedding[pred[i]]) == 0) {
    //     continue;  // Skip rows with all NaNs
    //   }
    //
    //   // Extract E1-dimensional embedding for unit pred[i]
    //   std::vector<double> xi_E1(embedding[pred[i]].begin(), embedding[pred[i]].begin() + E1);
    //
    //   double min_dist = std::numeric_limits<double>::max();
    //   size_t nn_idx = N;  // invalid index placeholder
    //
    //   // Find nearest neighbor of i in E1-dimensional space
    //   for (size_t j = 0; j < lib.size(); ++j) {
    //     if (pred[i] == lib[j] || checkOneDimVectorNotNanNum(embedding[lib[j]]) == 0) continue;
    //
    //     std::vector<double> xj_E1(embedding[lib[j]].begin(), embedding[lib[j]].begin() + E1);
    //     double dist = CppDistance(xi_E1, xj_E1, L1norm, true);  // true: skip NaNs
    //
    //     if (dist < min_dist) {
    //       min_dist = dist;
    //       nn_idx = lib[j];
    //     }
    //   }
    //
    //   if (nn_idx == N || min_dist == 0.0) continue;  // skip degenerate cases
    //
    //   // Compare E2-th coordinate difference (new dimension)
    //   double diff = std::abs(embedding[pred[i]][E2 - 1] - embedding[nn_idx][E2 - 1]);
    //   double ratio = diff / min_dist;
    //
    //   if (ratio > Rtol || diff > Atol) {
    //     ++false_count;
    //   }
    //   ++total;
    // }

    for (size_t i = 0; i < pred.size(); ++i) {
      if (checkOneDimVectorNotNanNum(embedding[pred[i]]) == 0) {
        continue;  // Skip rows with all NaNs
      }

      double min_dist = std::numeric_limits<double>::max();
      size_t nn_idx = N;  // invalid index placeholder

      // Find nearest neighbor of i in E1-dimensional space
      for (size_t j = 0; j < lib.size(); ++j) {
        if (pred[i] == lib[j] || checkOneDimVectorNotNanNum(embedding[lib[j]]) == 0) continue;

        double dist = distmat[pred[i]][lib[j]];

        if (dist < min_dist) {
          min_dist = dist;
          nn_idx = lib[j];
        }
      }

      if (nn_idx == N || min_dist == 0.0) continue;  // skip degenerate cases

      // Compare E2-th coordinate difference (new dimension)
      double diff = std::abs(embedding[pred[i]][E2 - 1] - embedding[nn_idx][E2 - 1]);
      double ratio = diff / min_dist;

      if (ratio > Rtol || diff > Atol) {
        ++false_count;
      }
      ++total;
    }

    return total > 0 ? static_cast<double>(false_count) / total
    : std::numeric_limits<double>::quiet_NaN();
  } else {
    // Parallel version: allocate one slot for each pred[i], thread-safe without locks
    std::vector<int> false_flags(pred.size(), -1); // -1 means skip or invalid, 0 means not a false neighbor, 1 means false neighbor

    RcppThread::parallelFor(0, pred.size(), [&](size_t i) {
      int pidx = pred[i];
      if (checkOneDimVectorNotNanNum(embedding[pidx]) == 0) return;

      double min_dist = std::numeric_limits<double>::max();
      int nn_idx = -1;

      for (size_t j = 0; j < lib.size(); ++j) {
        int lidx = lib[j];
        if (pidx == lidx || checkOneDimVectorNotNanNum(embedding[lidx]) == 0) continue;

        // Compute distance using only the first E1 dimensions
        std::vector<double> xi(embedding[pidx].begin(), embedding[pidx].begin() + E1);
        std::vector<double> xj(embedding[lidx].begin(), embedding[lidx].begin() + E1);
        double dist = CppDistance(xi, xj, L1norm, true);

        if (dist < min_dist) {
          min_dist = dist;
          nn_idx = lidx;
        }
      }

      // Skip if no neighbor found or minimum distance is zero
      if (nn_idx == -1 || min_dist == 0.0) return;

      // Compare the E2-th dimension to check for false neighbors
      double diff = std::abs(embedding[pidx][E2 - 1] - embedding[nn_idx][E2 - 1]);
      double ratio = diff / min_dist;

      // Determine if this is a false neighbor
      if (ratio > Rtol || diff > Atol) {
        false_flags[i] = 1;
      } else {
        false_flags[i] = 0;
      }
    }, threads); // use specified number of threads

    // After parallel section, aggregate results
    int false_count = 0, total = 0;
    for (int flag : false_flags) {
      if (flag >= 0) {
        total++;
        if (flag == 1) false_count++;
      }
    }

    if (total > 0) {
      return static_cast<double>(false_count) / total;
    } else {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }
}

/*
 * Compute False Nearest Neighbor (FNN) ratios across multiple embedding dimensions
 * for time series data.
 *
 * For a given embedding matrix (with each row representing a spatial unit and
 * each column an embedding dimension), this function evaluates the proportion
 * of false nearest neighbors (FNN) as the embedding dimension increases.
 *
 * It iteratively calls `CppSingleFNN` for each embedding dimension pair (E1, E2),
 * where E1 ranges from 1 to D - 1 (D = number of columns), and E2 = E1 + 1.
 * The FNN ratio measures how often a nearest neighbor in dimension E1 becomes
 * distant in dimension E2, suggesting that E1 is insufficient for reconstructing
 * the system.
 *
 * If `parallel_level == 0`, the function executes in serial;
 * otherwise, it uses multithreading to compute FNN ratios for each (E1, E2) pair
 * in parallel.
 *
 * Parameters:
 * - embedding: A vector of vectors where each row is a spatial unit’s embedding.
 *              Must have at least 2 columns (dimensions).
 * - lib: A vector of integer indices indicating the library set (0-based).
 * - pred: A vector of integer indices indicating the prediction set (0-based).
 * - Rtol: A vector of relative distance thresholds (one per E1).
 * - Atol: A vector of absolute distance thresholds (one per E1).
 * - L1norm: If true, use L1 (Manhattan) distance; otherwise, use L2 (Euclidean).
 * - threads: Number of threads to use for parallel computation.
 * - parallel_level: 0 for serial loop over E1, >0 for parallel loop over E1.
 *
 * Returns:
 * - A vector of FNN ratios corresponding to each E1 from 1 to D - 1.
 *   If not computable for a given E1, NaN is returned at that position.
 */
std::vector<double> CppFNN(const std::vector<std::vector<double>>& embedding,
                           const std::vector<int>& lib,
                           const std::vector<int>& pred,
                           const std::vector<double>& Rtol,
                           const std::vector<double>& Atol,
                           bool L1norm = false,
                           int threads = 8,
                           int parallel_level = 0) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  size_t max_E2 = embedding[0].size();
  std::vector<double> results(max_E2 - 1, std::numeric_limits<double>::quiet_NaN());

  if (embedding.empty() || embedding[0].size() < 2) {
    return results;  // Not enough dimensions to compute FNN
  }

  if (parallel_level == 0){
    // Loop through E1 = 1 to max_E2 - 1
    for (size_t E1 = 1; E1 < max_E2; ++E1) {
      size_t E2 = E1 + 1;
      double fnn_ratio = CppSingleFNN(embedding, lib, pred, E1, E2, threads_sizet,
                                      parallel_level, Rtol[E1 - 1], Atol[E1 - 1], L1norm);
      results[E1 - 1] = fnn_ratio;
    }
  } else {
    // Parallel computation
    RcppThread::parallelFor(1, max_E2, [&](size_t E1) {
      size_t E2 = E1 + 1;
      double fnn_ratio = CppSingleFNN(embedding, lib, pred, E1, E2, threads_sizet,
                                      parallel_level, Rtol[E1 - 1], Atol[E1 - 1], L1norm);
      results[E1 - 1] = fnn_ratio;
    }, threads_sizet);
  }

  return results;
}
