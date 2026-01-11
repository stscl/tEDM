#ifndef PCM_H
#define PCM_H

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
#include "DataStruct.h"

// Note: <RcppThread.h> is intentionally excluded from this header to avoid
//       unnecessary Rcpp dependencies and potential header inclusion order
//       issues (e.g., R.h being included before Rcpp headers). It should only
//       be included in the corresponding .cpp implementation file.
// #include <RcppThread.h>

/**
 * @brief Computes the partial correlation between the target variable and its simplex projection,
 *        incorporating control variables using a time-delay embedding approach.
 *
 * @param vectors: Reconstructed state-space, where each row represents a separate state vector.
 * @param target: Time series to be used as the target, aligned with 'vectors'.
 * @param controls: Time series data of control variables, stored row-wise.
 * @param lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 * @param pred_indices: Vector of indices indicating which states to predict from.
 * @param conEs: Vector specifying the number of dimensions for attractor reconstruction with control variables.
 * @param taus: Vector specifying the time lag step for constructing lagged state-space vectors with control variables.
 * @param num_neighbors: Vector specifying the numbers of neighbors to use for simplex projection.
 * @param cumulate: Flag indicating whether to cumulatively incorporate control variables.
 * @param dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average: Whether to average distance by the number of valid vector components.
 *
 * @return A std::vector<double> containing:
 *         - rho[0]: Pearson correlation between the target and its simplex projection.
 *         - rho[1]: Partial correlation controlling for the influence of the control variables.
 */
std::vector<double> PartialSimplex4TS(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    bool cumulate = false,
    int dist_metric = 2,
    bool dist_average = true
);

/**
 * @brief Computes the partial correlation between a time series and its prediction
 *        using the S-Map method, incorporating control variables.
 *
 * This function performs state-space reconstruction and S-Map prediction while accounting for
 * control variables in a time series. The process can be either cumulative or independent in
 * terms of incorporating control variables.
 *
 * @param vectors: Reconstructed state-space where each row represents a separate vector/state.
 * @param target: Time series used as the prediction target.
 * @param controls: Time series data of control variables, stored row-wise.
 * @param lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 * @param pred_indices: Vector of indices indicating which states to predict from.
 * @param conEs: Vector specifying the number of dimensions for attractor reconstruction with control variables.
 * @param taus: Vector specifying the stime lag step for constructing lagged state-space vectors with control variables.
 * @param num_neighbors: Vector specifying the numbers of neighbors to use for S-Map prediction.
 * @param theta: Weighting parameter for distances in S-Map.
 * @param cumulate: Boolean flag to determine whether to cumulate the partial correlations.
 * @param dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average: Whether to average distance by the number of valid vector components.
 * @return A vector of size 2 containing:
 *         - rho[0]: Pearson correlation between the target and its predicted values.
 *         - rho[1]: Partial correlation between the target and its predicted values, adjusting for control variables.
 */
std::vector<double> PartialSMap4TS(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    double theta = 1.0,
    bool cumulate = false,
    int dist_metric = 2,
    bool dist_average = true
);

/*
 * Perform partial cross mapping on a single library and prediction set.
 *
 * Parameters:
 *   - x_vectors: Reconstructed state-space (each row represents a separate vector/state).
 *   - y: Time series used as the target (should align with x_vectors).
 *   - controls: Time series data of control variables (stored by row).
 *   - lib_size: Size of the library used for cross mapping.
 *   - lib_indices: Vector of indices indicating which states to include when searching for neighbors.
 *   - pred_indices: Vector of indices indicating which states to predict from.
 *   - conEs: Number of dimensions for attractor reconstruction with control variables.
 *   - taus: Time lag step for constructing lagged state-space vectors with control variables.
 *   - b: A vector specifying the numbers of neighbors to use for simplex projection.
 *   - simplex: If true, uses simplex projection for prediction; otherwise, uses s-mapping.
 *   - theta: Distance weighting parameter for local neighbors in the manifold (used in s-mapping).
 *   - threads: The number of threads to use for parallel processing.
 *   - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 *   - cumulate: Whether to accumulate partial correlations.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components.
 *
 * Returns:
 *   A vector of PartialCorRes objects, where each contains:
 *   - An integer representing the library size.
 *   - A double representing the Pearson correlation coefficient (rho).
 *   - A double representing the Partial correlation coefficient (pratial rho).
 */
std::vector<PartialCorRes> PCMSingle(
    const std::vector<std::vector<double>>& x_vectors,  // Reconstructed state-space (each row is a separate vector/state)
    const std::vector<double>& y,                       // Time series to be used as the target (should line up with vectors)
    const std::vector<std::vector<double>>& controls,   // Time series data of control variables (**stored by row**)
    int lib_size,                                       // Size of the library
    const std::vector<int>& lib_indices,                // Indices of possible library states
    const std::vector<int>& pred_indices,               // Vector of indices indicating which states to predict from
    const std::vector<int>& conEs,                      // Number of dimensions for the attractor reconstruction with control variables
    const std::vector<int>& taus,                       // Time lag step for constructing lagged state-space vectors with control variables
    const std::vector<int>& b,                          // Numbers of neighbors to use for simplex projection
    bool simplex,                                       // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                       // Distance weighting parameter for the local neighbours in the manifold
    size_t threads,                                     // Number of threads to use for parallel processing
    int parallel_level,                                 // Level of parallel computing: 0 for `lower`, 1 for `higher`
    bool cumulate,                                      // Whether to cumulate the partial correlations
    int dist_metric,                                    // Distance metric selector (1: Manhattan, 2: Euclidean)
    bool dist_average                                   // Whether to average distance by the number of valid vector components
);

/**
 * Performs partial cross mapping on a time series dataset.
 *
 * Parameters:
 * - x: Time series used as the predictor variable (**cross mapping from**).
 * - y: Time series series used as the target variable (**cross mapping to**).
 * - controls: Time series data of control variables (**stored by row**).
 * - lib_sizes: A vector specifying different library sizes for SCPCM analysis.
 * - lib: A vector of representing the indices of spatial units to be the library.
 * - pred: A vector of representing the indices of spatial units to be predicted.
 * - Es: A vector specifying the embedding dimensions for attractor reconstruction using x and control variables.
 * - taus: A vector specifying the time lag steps for constructing lagged state-space vectors using x and control variables.
 * - b: A vector specifying the numbers of nearest neighbors used for prediction.
 * - simplex: Boolean flag indicating whether to use simplex projection (true) or S-mapping (false) for prediction.
 * - theta: Distance weighting parameter used for weighting neighbors in the S-mapping prediction.
 * - threads: Number of threads to use for parallel computation.
 * - cumulate: Boolean flag indicating whether to cumulate partial correlations.
 * - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 * - dist_average: Whether to average distance by the number of valid vector components.
 * - single_sig: Whether to estimate significance and confidence intervals using a single rho value.
 * - progressbar: Boolean flag indicating whether to display a progress bar during computation.
 *
 * Returns:
 *    A 2D vector of results, where each row contains:
 *      - The library size.
 *      - The mean pearson cross-mapping correlation.
 *      - The statistical significance of the pearson correlation.
 *      - The lower bound of the pearson correlation confidence interval.
 *      - The upper bound of the pearson correlation confidence interval.
 *      - The mean partial cross-mapping partial correlation.
 *      - The statistical significance of the partial correlation.
 *      - The lower bound of the partial correlation confidence interval.
 *      - The upper bound of the partial correlation confidence interval.
 */
std::vector<std::vector<double>> PCM(
    const std::vector<double>& x,                       // Time series to cross map from
    const std::vector<double>& y,                       // Time series to cross map to
    const std::vector<std::vector<double>>& controls,   // Time series data of control variables (**stored by row**)
    const std::vector<int>& lib_sizes,                  // Vector of library sizes to use
    const std::vector<int>& lib,                        // Vector specifying the library indices
    const std::vector<int>& pred,                       // Vector specifying the prediction indices
    const std::vector<int>& Es,                         // Number of dimensions for the attractor reconstruction with the x and control variables
    const std::vector<int>& taus,                       // Time lag step for constructing lagged state-space vectors with the x and control variables
    const std::vector<int>& b,                          // Numbers of nearest neighbors to use for prediction
    bool simplex,                                       // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                       // Distance weighting parameter for the local neighbours in the manifold
    int threads,                                        // Number of threads used from the global pool
    int parallel_level,                                 // Level of parallel computing: 0 for `lower`, 1 for `higher`
    bool cumulate,                                      // Whether to cumulate the partial correlations
    int dist_metric,                                    // Distance metric selector (1: Manhattan, 2: Euclidean)
    bool dist_average,                                  // Whether to average distance by the number of valid vector components
    bool single_sig,                                    // Whether to estimate significance and confidence intervals using a single rho value
    bool progressbar                                    // Whether to print the progress bar
);

#endif // PCM_H
