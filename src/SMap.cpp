#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"

/*
 * Computes the S-Map prediction using a reconstructed state-space representation.
 *
 * Parameters:
 *   - vectors: A 2D vector where each row is a reconstructed state vector.
 *   - target: A vector of taregt values corresponding to each state.
 *   - lib_indices: A vector of indices indicating which states to use as the library (neighbors).
 *   - pred_indices: A vector of indices indicating which states to make predictions for.
 *   - num_neighbors: Number of nearest neighbors to use in the S-Map algorithm.
 *   - theta: Distance weighting parameter for neighbor contributions.
 *
 * Returns: A vector<double> containing predicted target values at the positions in pred_indices.
 */
std::vector<double> SMapPrediction(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors,
    double theta
) {
  size_t N = target.size();
  std::vector<double> pred(N, std::numeric_limits<double>::quiet_NaN());

  if (num_neighbors <= 0 || lib_indices.empty() || pred_indices.empty()) {
    return pred;
  }

  size_t num_neighbors_sizet = static_cast<size_t>(num_neighbors);

  for (int pred_i : pred_indices) {
    std::vector<size_t> libs;
    for (int lib_i : lib_indices) {
      if (lib_i != pred_i) {
        libs.push_back(static_cast<size_t>(lib_i));
      }
    }

    if (libs.empty()) {
      pred[pred_i] = std::numeric_limits<double>::quiet_NaN();
      continue;
    }

    if (num_neighbors_sizet > libs.size()) {
      num_neighbors_sizet = libs.size();
    }

    // Compute distances
    std::vector<double> distances;
    for (size_t i : libs) {
      double sum_sq = 0.0;
      double sum_na = 0.0;
      for (size_t j = 0; j < vectors[pred_i].size(); ++j) {
        if (!std::isnan(vectors[i][j]) && !std::isnan(vectors[pred_i][j])) {
          sum_sq += std::pow(vectors[i][j] - vectors[pred_i][j], 2);
          sum_na += 1.0;
        }
      }
      distances.push_back((sum_na > 0) ? std::sqrt(sum_sq / sum_na) : std::numeric_limits<double>::quiet_NaN());
    }

    // Compute mean distance
    double mean_distance = 0.0;
    for (double dist : distances) {
      mean_distance += dist;
    }
    mean_distance /= distances.size();

    // Compute weights
    std::vector<double> weights(distances.size());
    for (size_t i = 0; i < distances.size(); ++i) {
      weights[i] = std::exp(-theta * distances[i] / mean_distance);
    }

    // Find nearest neighbors
    std::vector<size_t> neighbors(distances.size());
    std::iota(neighbors.begin(), neighbors.end(), 0);
    std::partial_sort(
      neighbors.begin(), neighbors.begin() + num_neighbors_sizet, neighbors.end(),
      [&](size_t a, size_t b) {
        return (distances[a] < distances[b]) ||
          (distances[a] == distances[b] && a < b);
      });

    // Prepare A matrix and b vector for weighted linear system
    std::vector<std::vector<double>> A(num_neighbors_sizet, std::vector<double>(vectors[pred_i].size() + 1, 0.0));
    std::vector<double> b(num_neighbors_sizet, 0.0);
    for (size_t i = 0; i < num_neighbors_sizet; ++i) {
      size_t idx = libs[neighbors[i]];
      for (size_t j = 0; j < vectors[pred_i].size(); ++j) {
        A[i][j] = vectors[idx][j] * weights[neighbors[i]];
      }
      A[i][vectors[pred_i].size()] = weights[neighbors[i]]; // bias term
      b[i] = target[idx] * weights[neighbors[i]];
    }

    // Perform SVD
    std::vector<std::vector<std::vector<double>>> svd_result = CppSVD(A);
    std::vector<std::vector<double>> U = svd_result[0];
    std::vector<double> S = svd_result[1][0];
    std::vector<std::vector<double>> V = svd_result[2];

    // Invert singular values with threshold
    double max_s = *std::max_element(S.begin(), S.end());
    std::vector<double> S_inv(S.size(), 0.0);
    for (size_t i = 0; i < S.size(); ++i) {
      if (S[i] >= max_s * 1e-5) {
        S_inv[i] = 1.0 / S[i];
      }
    }

    // Compute map coefficients
    std::vector<double> map_coeffs(vectors[pred_i].size() + 1, 0.0);
    for (size_t i = 0; i < V.size(); ++i) {
      for (size_t j = 0; j < S_inv.size(); ++j) {
        map_coeffs[i] += V[i][j] * S_inv[j] * U[j][i];
      }
    }

    // Multiply map coefficients by b
    for (size_t i = 0; i < map_coeffs.size(); ++i) {
      map_coeffs[i] *= b[i];
    }

    // Generate prediction
    double prediction = 0.0;
    for (size_t i = 0; i < vectors[pred_i].size(); ++i) {
      prediction += map_coeffs[i] * vectors[pred_i][i];
    }
    prediction += map_coeffs[vectors[pred_i].size()]; // bias term
    pred[pred_i] = prediction;
  }

  return pred;
}

/*
 * Computes the Rho value using the 'S-Maps' prediction method.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Time series data vector to be predicted.
 *   - lib_indices: Vector of integer indices specifying which states to use for finding neighbors.
 *   - pred_indices: Vector of integer indices specifying which states to predict.
 *   - num_neighbors: Number of neighbors to use for S-Map.
 *   - theta: Weighting parameter for distances.
 *
 * Returns: The Pearson correlation coefficient (Rho) between predicted and actual values.
 */
double SMap(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors,
    double theta
) {
  double rho = std::numeric_limits<double>::quiet_NaN();

  // Call SMapPrediction to get the prediction results
  std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    rho = PearsonCor(target_pred, target, true);
  }
  return rho;
}


/*
 * Computes the S-Map prediction and evaluates prediction performance.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Time series data vector to be predicted.
 *   - lib_indices: Vector of integer indices specifying which states to use for finding neighbors.
 *   - pred_indices: Vector of integer indices specifying which states to predict.
 *   - num_neighbors: Number of neighbors to use for S-Map.
 *   - theta: Weighting parameter for distances.
 *
 * Returns: A vector<double> containing {Pearson correlation, MAE, RMSE}.
 */
std::vector<double> SMapBehavior(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors,
    double theta
) {
  // Initialize PearsonCor, MAE, and RMSE
  double pearson = std::numeric_limits<double>::quiet_NaN();
  double mae = std::numeric_limits<double>::quiet_NaN();
  double rmse = std::numeric_limits<double>::quiet_NaN();

  // Call SMapPrediction to get the prediction results
  std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors, theta);

  if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
    // Compute PearsonCor, MAE, and RMSE
    pearson = PearsonCor(target_pred, target, true);
    mae = CppMAE(target_pred, target, true);
    rmse = CppRMSE(target_pred, target, true);
  }

  // Return the three metrics as a vector
  return {pearson, mae, rmse};
}
