#ifndef CCM_H
#define CCM_H

#include <vector>
#include <cmath>
#include <algorithm> // Include for std::partial_sort
#include <numeric>
#include <utility>
#include <limits>
#include <map>
#include "CppStats.h"
#include "Embed.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

/*
 * Perform convergent cross mapping on a single lib and pred.
 *
 * Parameters:
 *   - x_vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - y: time series used as the target (should align with x_vectors).
 *   - lib_size: Size of the library used for cross mapping.
 *   - lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of indices indicating which states to predict from.
 *   - b: Number of neighbors to use for simplex projection.
 *   - simplex: If true, uses simplex projection for prediction; otherwise, uses s-mapping.
 *   - theta: Distance weighting parameter for local neighbors in the manifold (used in s-mapping).
 *   - threads: The number of threads to use for parallel processing.
 *   - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components.
 *
 * Returns:
 *   A vector of pairs, where each pair consists of:
 *   - An integer representing the library size.
 *   - A double representing the Pearson correlation coefficient (rho) between predicted and actual values.
 */
std::vector<std::pair<int, double>> CCMSingle(
    const std::vector<std::vector<double>>& x_vectors,
    const std::vector<double>& y,
    int lib_size,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int b,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    int dist_metric,
    bool dist_average
);

/**
 * Performs convergent cross mapping on a time series data.
 *
 * Parameters:
 * - x: time series used as the predict variable (**cross mapping from**).
 * - y: time series used as the target variable (**cross mapping to**).
 * - lib_sizes: A vector specifying different library sizes for CCM analysis.
 * - lib: A vector of representing the indices of time points to be the library.
 * - pred: A vector of representing the indices of time points to be predicted.
 * - E: Embedding dimension for attractor reconstruction.
 * - tau: the step of stime lags for prediction.
 * - b: Number of nearest neighbors used for prediction.
 * - simplex: Boolean flag indicating whether to use simplex projection (true) or S-mapping (false) for prediction.
 * - theta: Distance weighting parameter used for weighting neighbors in the S-mapping prediction.
 * - threads: Number of threads to use for parallel computation.
 * - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 *   dist_metric    - Distance metric selector (1: Manhattan, 2: Euclidean).
 *   dist_average   - Whether to average distance by the number of valid vector components.
 * - progressbar: Boolean flag to indicate whether to display a progress bar during computation.
 *
 * Returns:
 *    A 2D vector of results, where each row contains:
 *      - The library size.
 *      - The mean cross-mapping correlation.
 *      - The statistical significance of the correlation.
 *      - The upper bound of the confidence interval.
 *      - The lower bound of the confidence interval.
 */
std::vector<std::vector<double>> CCM(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<int>& lib_sizes,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    int E = 3,
    int tau = 1,
    int b = 4,
    bool simplex = true,
    double theta = 1.0,
    int threads = 8,
    int parallel_level = 0,
    int dist_metric = 2,
    bool dist_average = true,
    bool progressbar = false
);

#endif // CCM_H
