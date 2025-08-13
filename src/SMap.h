#ifndef SMap_H
#define SMap_H

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include "CppStats.h"

/**
 * @brief Perform S-Mapping prediction using locally weighted linear regression.
 *
 * This function performs prediction based on a reconstructed state-space (time-delay embedding).
 * For each prediction index, it:
 *   - Finds the nearest neighbors from the library indices, excluding the current prediction index.
 *   - Computes distance-based weights using the S-map weighting parameter (theta).
 *   - Constructs a locally weighted linear regression model using the valid neighbors.
 *   - Predicts the target value using the derived local model.
 *
 * @param vectors        A 2D matrix where each row is a reconstructed state vector (embedding).
 * @param target         A vector of scalar values to predict (e.g., time series observations).
 * @param lib_indices    Indices of the vectors used as the library (neighbor candidates).
 * @param pred_indices   Indices of the vectors used for prediction.
 * @param num_neighbors  Number of nearest neighbors to use in local regression. Default is 4.
 * @param theta          Weighting parameter controlling exponential decay of distances. Default is 1.0.
 * @param dist_metric    Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 * @param dist_average   Whether to average distance by the number of valid vector components. Default is true.
 * @return std::vector<double> Predicted values aligned with the input target vector.
 *         Entries at non-prediction indices or with insufficient valid neighbors are NaN.
 */
std::vector<double> SMapPrediction(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    double theta = 1.0,
    int dist_metric = 2,
    bool dist_average = true
);

/*
 * Computes the Rho value using the 'S-Mapping' prediction method.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Time series data vector to be predicted.
 *   - lib_indices: Vector of integer indices specifying which states to use for finding neighbors.
 *   - pred_indices: Vector of integer indices specifying which states to predict.
 *   - num_neighbors: Number of neighbors to use for S-Map. Default is 4.
 *   - theta: Weighting parameter for distances. Default is 1.0.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *
 * Returns: The Pearson correlation coefficient (Rho) between predicted and actual values.
 */
double SMap(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    double theta = 1.0,
    int dist_metric = 2,
    bool dist_average = true
);

/*
 * Computes the S-Mapping prediction and evaluates prediction performance.
 *
 * Parameters:
 *   - vectors: Reconstructed state-space (each row is a separate vector/state).
 *   - target: Time series data vector to be predicted.
 *   - lib_indices: Vector of integer indices specifying which states to use for finding neighbors.
 *   - pred_indices: Vector of integer indices specifying which states to predict.
 *   - num_neighbors: Number of neighbors to use for S-Map. Default is 4.
 *   - theta: Weighting parameter for distances. Default is 1.0.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *
 * Returns: A vector<double> containing {Pearson correlation, MAE, RMSE}.
 */
std::vector<double> SMapBehavior(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int num_neighbors = 4,
    double theta = 1.0,
    int dist_metric = 2,
    bool dist_average = true
);

#endif // SMap_H
