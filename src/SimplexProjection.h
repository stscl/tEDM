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
