#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"

/*
 * Perform simplex projection prediction using state-space reconstruction.
 *
 * Given reconstructed vectors and a target time series, this function predicts
 * target values at specified prediction indices by weighting nearest neighbors
 * in the library set. Distances involving NaN components are excluded to ensure
 * numerical stability.
 *
 * Parameters:
 *   vectors     - A 2D vector representing reconstructed state-space vectors.
 *                 Each element vectors[i] is a vector/state at time i.
 *   target      - The time series values corresponding to each vector.
 *   lib_indices - Indices specifying which states to use as library (neighbors).
 *   pred_indices- Indices specifying which states to make predictions for.
 *   num_neighbors - Number of nearest neighbors to consider in prediction.
 *
 * Returns:
 *   A vector<double> of predicted target values aligned with input target size.
 *   Entries are NaN if prediction is not possible for that index.
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

  if (num_neighbors <= 0) {
    return pred;  // no neighbors to use, return all NaNs
  }

  for (size_t pi = 0; pi < pred_indices.size(); ++pi) {
    int p = pred_indices[pi];

    // // Skip if target at prediction index is NaN
    // if (std::isnan(target[p])) {
    //   continue;
    // }

    // Compute distances only for valid vector pairs (exclude NaNs)
    std::vector<double> distances;
    std::vector<int> valid_libs;  // keep track of libs corresponding to valid distances
    for (int i : lib_indices) {
      if (i == p) continue; // Skip self-matching

      double sum_sq = 0.0;
      double count = 0.0;
      for (size_t j = 0; j < vectors[p].size(); ++j) {
        if (!std::isnan(vectors[i][j]) && !std::isnan(vectors[p][j])) {
          sum_sq += std::pow(vectors[i][j] - vectors[p][j], 2);
          count += 1.0;
        }
      }
      if (count > 0) {
        distances.push_back(std::sqrt(sum_sq / count));
        valid_libs.push_back(i);
      }
    }

    // If no valid distances found, prediction is NaN
    if (distances.empty()) {
      continue;
    }

    // Adjust number of neighbors to available valid libs
    size_t k = std::min(static_cast<size_t>(num_neighbors), distances.size());

    // Prepare indices for sorting
    std::vector<size_t> neighbors(distances.size());
    std::iota(neighbors.begin(), neighbors.end(), 0);

    // Partial sort to find k nearest neighbors by distance
    std::partial_sort(
      neighbors.begin(), neighbors.begin() + k, neighbors.end(),
      [&](size_t a, size_t b) {
        return (distances[a] < distances[b]) ||
          (distances[a] == distances[b] && a < b);
      });

    double min_distance = distances[neighbors[0]];

    // Compute weights for neighbors
    std::vector<double> weights(k);
    if (min_distance == 0.0) {  // Perfect match found
      std::fill(weights.begin(), weights.end(), 0.000001);
      for (size_t i = 0; i < k; ++i) {
        if (distances[neighbors[i]] == 0.0) {
          weights[i] = 1.0;
        }
      }
    } else {
      for (size_t i = 0; i < k; ++i) {
        weights[i] = std::exp(-distances[neighbors[i]] / min_distance);
        if (weights[i] < 0.000001) {
          weights[i] = 0.000001;
        }
      }
    }

    double total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);

    // Calculate weighted prediction, ignoring any NaN targets
    // (No NaNs here, as NaN values in the corresponding components of lib and pred are excluded in advance.)
    double prediction = 0.0;
    for (size_t i = 0; i < k; ++i) {
      prediction += weights[i] * target[valid_libs[neighbors[i]]];
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
