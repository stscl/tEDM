#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"

/*
 * Computes predictions using the simplex projection method based on state-space reconstruction.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - target: Time series used as the target (should align with vectors).
 *   - lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of indices indicating which states to predict from.
 *   - num_neighbors: Number of neighbors to use for simplex projection.
 *
 * Returns:
 *   A vector<double> containing predicted target values (target_pred).
 */
std::vector<double> SimplexProjectionPrediction(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors
) {
  size_t N = target.size();
  std::vector<double> pred(N, std::numeric_limits<double>::quiet_NaN());

  // no neighbor to use, return all nan
  if (num_neighbors <= 0){
    return pred;
  }

  // Make predictions
  for (size_t pi = 0; pi < pred_indices.size(); ++pi) {
    int p = pred_indices[pi];

    // Skip if target at prediction index is NaN
    if (std::isnan(target[p])) {
      continue;
    }

    // Create lib set excluding the current prediction index p
    std::vector<int> libs;
    for (int idx : lib_indices) {
      if (idx != p) {
        libs.push_back(idx);
      }
    }

    // Handle the case where libs is empty
    if (libs.empty()) {
      pred[p] = std::numeric_limits<double>::quiet_NaN();
      continue;
    }

    // Adjust num_neighbors if it exceeds libs.size()
    size_t k = std::min(static_cast<size_t>(num_neighbors), libs.size());

    // Compute distances
    std::vector<double> distances;
    for (int i : libs) {
      double sum_sq = 0.0;
      double sum_na = 0.0;
      for (size_t j = 0; j < vectors[p].size(); ++j) {
        if (!std::isnan(vectors[i][j]) && !std::isnan(vectors[p][j])) {
          sum_sq += std::pow(vectors[i][j] - vectors[p][j], 2);
          sum_na += 1.0;
        }
      }

      if (sum_na > 0) {
        distances.push_back(std::sqrt(sum_sq / sum_na));
      } else {
        distances.push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }

    // Find nearest neighbors
    std::vector<size_t> neighbors(libs.size());
    std::iota(neighbors.begin(), neighbors.end(), 0);
    std::partial_sort(neighbors.begin(), neighbors.begin() + k, neighbors.end(),
                      [&](size_t a, size_t b) {
                        return (distances[a] < distances[b]) ||
                          (distances[a] == distances[b] && a < b);
                      });

    double min_distance = distances[neighbors[0]];

    // Compute weights
    std::vector<double> weights(k);
    if (min_distance == 0) { // Perfect match
      std::fill(weights.begin(), weights.end(), 0.000001);
      for (size_t i = 0; i < k; ++i) {
        if (distances[neighbors[i]] == 0) weights[i] = 1.0;
      }
    } else {
      for (size_t i = 0; i < k; ++i) {
        weights[i] = std::exp(-distances[neighbors[i]] / min_distance);
        if (weights[i] < 0.000001) weights[i] = 0.000001;
      }
    }

    double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);

    // Make prediction
    double prediction = 0.0;
    for (size_t i = 0; i < k; ++i) {
      prediction += weights[i] * target[libs[neighbors[i]]];
    }

    pred[p] = prediction / total_weight;
  }

  return pred;
}

/*
 * Computes the Pearson correlation coefficient (rho) using the simplex projection prediction method.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - target: Time series used as the target (should align with vectors).
 *   - lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of indices indicating which states to use for prediction.
 *   - num_neighbors: Number of neighbors to use for simplex projection.
 *
 * Returns:
 *   A double representing the Pearson correlation coefficient (rho) between the predicted and actual target values.
 */
double SimplexProjection(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors
) {
  double rho = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    rho = PearsonCor(target_pred, target, true);
  }
  return rho;
}

/*
 * Computes the simplex projection and evaluates prediction performance.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Time series to be used as the target (should align with vectors).
 *   - lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of indices indicating which states to predict from.
 *   - num_neighbors: Number of neighbors to use for simplex projection.
 *
 * Returns:
 *   A vector<double> containing:
 *     - Pearson correlation coefficient (PearsonCor)
 *     - Mean absolute error (MAE)
 *     - Root mean squared error (RMSE)
 */
std::vector<double> SimplexBehavior(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors
) {
  double pearson = std::numeric_limits<double>::quiet_NaN();
  double mae = std::numeric_limits<double>::quiet_NaN();
  double rmse = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    pearson = PearsonCor(target_pred, target, true);
    mae = CppMAE(target_pred, target, true);
    rmse = CppRMSE(target_pred, target, true);
  }

  return {pearson, mae, rmse};
}
