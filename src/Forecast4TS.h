#ifndef Forecast4TS_H
#define Forecast4TS_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "Embed.h"
#include "SimplexProjection.h"
#include "SMap.h"
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
 *   - tau: The time lag step for constructing lagged state-space vectors.
 *   - threads: Number of threads used from the global pool.
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
                                            int tau,
                                            int threads);

/*
 * Evaluates prediction performance of different theta parameters for time series data using the s-mapping method.
 *
 * Parameters:
 *   - source: A vector to be embedded.
 *   - target: A vector to be predicted.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - theta: A vector of weighting parameters for distance calculation in SMap.
 *   - E: The embedding dimension to evaluate.
 *   - tau: The time lag step for constructing lagged state-space vectors.
 *   - b: Number of nearest neighbors to use for prediction.
 *   - threads: Number of threads used from the global pool.
 *
 * Returns:
 *   A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value.
 */
std::vector<std::vector<double>> SMap4TS(const std::vector<double>& source,
                                         const std::vector<double>& target,
                                         const std::vector<int>& lib_indices,
                                         const std::vector<int>& pred_indices,
                                         const std::vector<double>& theta,
                                         int E,
                                         int tau,
                                         int b,
                                         int threads);

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
 * @param tau            Time delay between embedding dimensions.
 * @param threads        Number of threads for parallel computation.
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
                                                 int tau,
                                                 int threads);

#endif // Forecast4TS_H
