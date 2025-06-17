#ifndef SimplexProjection_H
#define SimplexProjection_H

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
 *   - target: time series used as the target (should align with vectors).
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
);

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
);

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
);

#endif // SimplexProjection_H
