#ifndef Forecast4TS_H
#define Forecast4TS_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "Embed.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionCardinality.h"
#include <RcppThread.h>

/*
 * Evaluates prediction performance of different combinations of embedding dimensions and number of nearest neighbors
 * for time series data using simplex projection.
 *
 * Parameters:
 *   - source: A vector to be embedded.
 *   - target: A vector to be predicted.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - E: A vector of embedding dimensions to evaluate.
 *   - b: A vector of nearest neighbor values to evaluate.
 *   - tau: The time lag step for constructing lagged state-space vectors. Default is 1.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *   - threads: Number of threads used from the global pool. Default is 8.
 *
 * Returns:
 *   A 2D vector where each row contains [E, b, rho, mae, rmse] for a given combination of E and b.
 */
std::vector<std::vector<double>> Simplex4TS(const std::vector<double>& source,
                                            const std::vector<double>& target,
                                            const std::vector<int>& lib_indices,
                                            const std::vector<int>& pred_indices,
                                            const std::vector<int>& E,
                                            const std::vector<int>& b,
                                            int tau = 1,
                                            int dist_metric = 2,
                                            bool dist_average = true,
                                            int threads = 8);

/*
 * Evaluates prediction performance of different theta parameters for time series data using the s-mapping method.
 *
 * Parameters:
 *   - source: A vector to be embedded.
 *   - target: A vector to be predicted.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - theta: A vector of weighting parameters for distance calculation in SMap.
 *   - E: The embedding dimension to evaluate. Default is 3.
 *   - tau: The time lag step for constructing lagged state-space vectors. Default is 1.
 *   - b: Number of nearest neighbors to use for prediction. Default is 4.
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *   - threads: Number of threads used from the global pool. Default is 8.
 *
 * Returns:
 *   A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value.
 */
std::vector<std::vector<double>> SMap4TS(const std::vector<double>& source,
                                         const std::vector<double>& target,
                                         const std::vector<int>& lib_indices,
                                         const std::vector<int>& pred_indices,
                                         const std::vector<double>& theta,
                                         int E = 3,
                                         int tau = 1,
                                         int b = 4,
                                         int dist_metric = 2,
                                         bool dist_average = true,
                                         int threads = 8);

/**
 * Perform multivariate Simplex projection prediction across multiple time series datasets
 * for a range of embedding dimensions (E) and neighbor numbers (b), using a parallelized
 * implementation. The function embeds all source time series, combines all observations,
 * and evaluates the forecasting skill on the target time series using given training (lib)
 * and testing (pred) indices.
 *
 * @param source         A vector of time series vectors used for embedding.
 * @param target         A vector of time series vectors representing prediction targets.
 * @param lib_indices    Indices of time series used as training libraries (1-based).
 * @param pred_indices   Indices of time series used as prediction sets (1-based).
 * @param E              Vector of embedding dimensions to evaluate.
 * @param b              Vector of number of nearest neighbors to evaluate.
 * @param tau            Time delay between embedding dimensions. Default is 1.
 * @param dist_metric    Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 * @param dist_average   Whether to average distance by the number of valid vector components. Default is true.
 * @param threads        Number of threads for parallel computation. Default is 8.
 *
 * @return A vector of vectors where each sub-vector contains:
 *         [E, b, Pearson correlation, Mean Absolute Error (MAE), Root Mean Square Error (RMSE)]
 *         for one (E, b) combination.
 */
std::vector<std::vector<double>> MultiSimplex4TS(const std::vector<std::vector<double>>& source,
                                                 const std::vector<std::vector<double>>& target,
                                                 const std::vector<int>& lib_indices,
                                                 const std::vector<int>& pred_indices,
                                                 const std::vector<int>& E,
                                                 const std::vector<int>& b,
                                                 int tau = 1,
                                                 int dist_metric = 2,
                                                 bool dist_average = true,
                                                 int threads = 8);

/**
 * Compute Intersection Cardinality AUC over Lattice Embedding Settings.
 *
 * This function computes the causal strength between two lattice-structured time series
 * (`source` and `target`) by evaluating the Intersection Cardinality (IC) curve, and
 * summarizing it using the Area Under the Curve (AUC) metric.
 *
 * For each combination of embedding dimension `E` and neighbor size `b`, the function:
 *  - Generates state-space embeddings based on lattice neighborhood topology.
 *  - Filters out prediction points with missing (NaN) values.
 *  - Computes neighbor structures and evaluates intersection sizes between the mapped
 *    neighbors of `source` and `target`.
 *  - Aggregates the IC curve and estimates the AUC (optionally using significance test).
 *
 * @param source         Time series values of the potential cause variable (flattened lattice vector).
 * @param target         Time series values of the potential effect variable (same shape as `source`).
 * @param lib_indices    Indices used for library (training) data.
 * @param pred_indices   Indices used for prediction (testing) data.
 * @param E              Vector of embedding dimensions to try.
 * @param b              Vector of neighbor sizes to try.
 * @param tau            Embedding delay (usually 1 for lattice).
 * @param exclude        Number of nearest neighbors to exclude (e.g., temporal or spatial proximity).
 * @param dist_metric    Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param threads        Number of threads for parallel computation.
 * @param parallel_level Flag indicating whether to use multi-threading (0: serial, 1: parallel).
 *
 * @return A vector of size `E.size() * b.size()`, each element is a vector:
 *         [embedding_dimension, neighbor_size, auc_value].
 *         If inputs are invalid or no prediction point is valid, the AUC value is NaN.
 *
 * @note
 *   - Only AUC and p value are returned in current version. Use other utilities to derive CI.
 *   - Library and prediction indices should be adjusted for 0-based indexing before calling.
 *   - Lattice embedding assumes neighborhood-based spatial structure.
 */
std::vector<std::vector<double>> IC4TS(const std::vector<double>& source,
                                       const std::vector<double>& target,
                                       const std::vector<size_t>& lib_indices,
                                       const std::vector<size_t>& pred_indices,
                                       const std::vector<int>& E,
                                       const std::vector<int>& b,
                                       int tau = 1,
                                       int exclude = 0,
                                       int dist_metric = 2,
                                       int threads = 8,
                                       int parallel_level = 0);

#endif // Forecast4TS_H
