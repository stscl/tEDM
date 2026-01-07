#ifndef SimplexProjection_H
#define SimplexProjection_H

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "NumericUtils.h"
#include "CppStats.h"

/*
 * Perform simplex projection prediction using state-space reconstruction.
 *
 * Given reconstructed state-space vectors and corresponding target values,
 * this function predicts target values at specified prediction indices by
 * weighting nearest neighbors from a given library set.
 *
 * Distance calculations exclude NaN components to ensure numerical stability.
 * Supports two distance metrics and optional averaging of distance by the
 * number of valid vector components.
 *
 * Supported distance metrics:
 *   dist_metric = 1: L1 (Manhattan) distance
 *   dist_metric = 2: L2 (Euclidean) distance
 *
 * Parameters:
 *   vectors        - 2D vector of reconstructed state-space vectors; vectors[i] corresponds to time point i.
 *   target         - Target values for each time point.
 *   lib_indices    - Indices specifying states used as the library (neighbors).
 *   pred_indices   - Indices specifying states to predict.
 *   num_neighbors  - Number of nearest neighbors considered for prediction. Default is 4.
 *   dist_metric    - Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   dist_average   - Whether to average distance by the number of valid vector components. Default is true.
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
    int num_neighbors = 4,
    int dist_metric = 2,
    bool dist_average = true
);

/*
 * Computes the Pearson correlation coefficient (rho) using the simplex projection prediction method.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - target: Time series used as the target (should align with vectors).
 *   - lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of indices indicating which states to use for prediction.
 *   - num_neighbors: Number of neighbors to use for simplex projection. Default is 4.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *
 * Returns:
 *   A double representing the Pearson correlation coefficient (rho) between the predicted and actual target values.
 */
double SimplexProjection(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    int dist_metric = 2,
    bool dist_average = true
);

/*
 * Computes the simplex projection and evaluates prediction performance.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Time series to be used as the target (should align with vectors).
 *   - lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of indices indicating which states to predict from.
 *   - num_neighbors: Number of neighbors to use for simplex projection. Default is 4.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
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
    int num_neighbors = 4,
    int dist_metric = 2,
    bool dist_average = true
);

#endif // SimplexProjection_H
