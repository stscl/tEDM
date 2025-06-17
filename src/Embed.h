#ifndef Embed_H
#define Embed_H

#include <vector>
#include <cmath>
#include <limits>

/**
 * @brief Generate time-delay embeddings for a univariate time series.
 *
 * This function reconstructs the state space of a scalar time series
 * using time-delay embedding with dimension E and lag tau.
 *
 * - When tau = 0, embedding uses lags of 0, 1, ..., E-1.
 * - When tau > 0, embedding uses lags of tau, 2*tau, ..., E*tau.
 *
 * All values are pre-initialized to NaN. Elements are filled only when
 * sufficient non-NaN lagged values are available. Columns containing only
 * NaN values are removed before returning.
 *
 * @param vec The input time series as a vector of doubles.
 * @param E Embedding dimension.
 * @param tau Time lag.
 * @return A 2D vector (matrix) with valid embeddings (rows × cols).
 */
std::vector<std::vector<double>> Embed(
    const std::vector<double>& vec,
    int E,
    int tau
);

#endif // Embed_H
