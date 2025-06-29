#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "CppStats.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/*
 * Compute the False Nearest Neighbors (FNN) ratio for spatial cross-sectional data.
 *
 * This function determines whether nearest neighbors identified in a lower-dimensional
 * embedded space (E1) remain close in a higher-dimensional space (E2).
 * If not, the neighbor is considered a "false" neighbor, indicating the need for
 * a higher embedding dimension to accurately capture spatial proximity.
 *
 * Parameters:
 * - embedding: A matrix (vector of vectors) representing the spatial embedding,
 *              where each row corresponds to a spatial unit's attributes.
 *              Must contain at least E2 columns.
 * - lib: Library index vector (1-based in R, converted to 0-based).
 * - pred: Prediction index vector (1-based in R, converted to 0-based).
 * - E1: The base embedding dimension used to identify the nearest neighbor (E1 < E2).
 * - E2: The full embedding dimension used to test false neighbors (usually E1 + 1).
 * - Rtol: Relative threshold (default 10.0). If the change in the added dimension is
 *         large relative to the E1 distance, the neighbor is considered false.
 * - Atol: Absolute threshold (default 2.0). If the added dimension changes too much
 *         absolutely, the neighbor is also considered false.
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
                    double Rtol = 10.0,
                    double Atol = 2.0,
                    bool L1norm = false) {
  if (embedding.empty() || embedding[0].size() < E2) {
    return std::numeric_limits<double>::quiet_NaN();  // Invalid dimensions
  }

  size_t N = embedding.size();

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
}

/*
 * Compute False Nearest Neighbor (FNN) ratios for a range of embedding dimensions
 * on spatial cross-sectional data.
 *
 * This function iteratively calls CppSingleFNN using pairs of embedding dimensions:
 * from 1 to D-1 (where D is the number of columns in the embedding matrix), such that
 * each E1 = d and E2 = d + 1.
 *
 * It is used to identify the optimal embedding dimension by observing how the proportion
 * of false nearest neighbors decreases with higher dimensions.
 *
 * Parameters:
 * - embedding: A matrix (vector of vectors) where each row is a spatial unit's
 *              multidimensional embedding. Should have at least 2 columns.
 * - lib: Library index vector (1-based in R, converted to 0-based).
 * - pred: Prediction index vector (1-based in R, converted to 0-based).
 * - Rtol: Vectors of relative distance threshold.
 * - Atol: Vectors of absolute distance threshold.
 * - L1norm: Whether to use L1 (Manhattan) distance instead of L2 (Euclidean).
 * - threads: Number of parallel threads.
 *
 * Returns:
 * - A vector of doubles, where each entry corresponds to the FNN ratio for a given
 *   E1 (from 1 to D - 1), with E2 = E1 + 1. If the result is invalid, NaN is returned.
 *   Returns an NaN vector if the embedding has fewer than 2 columns.
 */
std::vector<double> CppFNN(const std::vector<std::vector<double>>& embedding,
                           const std::vector<int>& lib,
                           const std::vector<int>& pred,
                           const std::vector<double>& Rtol,
                           const std::vector<double>& Atol,
                           bool L1norm, int threads) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  size_t max_E2 = embedding[0].size();
  std::vector<double> results(max_E2 - 1, std::numeric_limits<double>::quiet_NaN());

  if (embedding.empty() || embedding[0].size() < 2) {
    return results;  // Not enough dimensions to compute FNN
  }

  // // Loop through E1 = 1 to max_E2 - 1
  // for (size_t E1 = 1; E1 < max_E2; ++E1) {
  //   size_t E2 = E1 + 1;
  //   double fnn_ratio = CppSingleFNN(embedding, lib, pred, E1, E2,
  //                                   Rtol[E1 - 1], Atol[E1 - 1], L1norm);
  //   results[E1 - 1] = fnn_ratio;
  // }

  // Parallel computation
  RcppThread::parallelFor(1, max_E2, [&](size_t E1) {
    size_t E2 = E1 + 1;
    double fnn_ratio = CppSingleFNN(embedding, lib, pred, E1, E2,
                                    Rtol[E1 - 1], Atol[E1 - 1], L1norm);
    results[E1 - 1] = fnn_ratio;
  }, threads_sizet);

  return results;
}
