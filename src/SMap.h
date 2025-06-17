#ifndef SMap_H
#define SMap_H

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
);

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
);

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
);

#endif // SMap_H
