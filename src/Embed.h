#ifndef Embed_H
#define Embed_H

#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

/**
 * @brief Generate time-delay embeddings for a univariate time series.
 *
 * This function reconstructs the state space of a scalar time series
 * using time-delay embedding with dimension E and lag tau.
 *
 * - When tau = 0, embedding uses lags of 0, 1, ..., E-1.  (original behavior)
 * - When tau > 0 and style = 1, embedding uses lags of tau, 2*tau, ..., E*tau.
 * - When tau > 0 and style = 0, embedding uses lags of 0, tau, 2*tau, ..., (E-1)*tau.
 *
 * Example:
 * Input: vec = {1, 2, 3, 4, 5}, E = 3, tau = 0
 * Output:
 * 1    NaN    NaN
 * 2    1      NaN
 * 3    2      1
 * 4    3      2
 * 5    4      3
 *
 * All values are pre-initialized to NaN. (Elements are filled only when
 * sufficient non-NaN lagged values are available. *Previously bound,
 * now abandoned*) Columns containing only NaN values are removed before
 * returning. If no valid embedding columns remain (due to short input
 * and large E/tau), an exception is thrown.
 *
 * @param vec The input time series as a vector of doubles.
 * @param E Embedding dimension.
 * @param tau Time lag.
 * @param style Lag style when tau > 0:
 *        - style = 1: tau, 2*tau, ..., E*tau
 *        - style = 0: 0, tau, 2*tau, ..., (E-1)*tau
 * @return A 2D vector (matrix) with valid embeddings (rows × cols).
 */
std::vector<std::vector<double>> Embed(
    const std::vector<double>& vec,
    int E = 3,
    int tau = 1,
    int style = 0
);

#endif // Embed_H
