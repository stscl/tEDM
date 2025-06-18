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
 * - When tau = 0, embedding uses lags of 0, 1, ..., E-1.
 * - When tau > 0, embedding uses lags of tau, 2*tau, ..., E*tau.
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
 * @return A 2D vector (matrix) with valid embeddings (rows × cols).
 */
std::vector<std::vector<double>> Embed(
    const std::vector<double>& vec,
    int E,
    int tau
) {
  const size_t N = vec.size();
  const double NaN = std::numeric_limits<double>::quiet_NaN();

  // Preallocate embedding matrix: N rows, E columns
  std::vector<std::vector<double>> mat(N, std::vector<double>(E, NaN));

  for (size_t t = 0; t < N; ++t) {
    for (int j = 0; j < E; ++j) {
      int lag = (tau == 0) ? j : (j + 1) * tau;
      int idx = static_cast<int>(t) - lag;
      if (idx >= 0 && idx < static_cast<int>(N)) {
        mat[t][j] = vec[idx];
        // else leave NaN
      }
    }
  }

  // // Compute the actual lag step
  // for (size_t t = 0; t < N; ++t) {
  //   bool valid = true;
  //   for (int j = 0; j < E; ++j) {
  //     int lag = (tau == 0) ? j : (j + 1) * tau;
  //     int idx = static_cast<int>(t) - lag;
  //
  //     if (idx < 0 || idx >= static_cast<int>(N)) {
  //       valid = false;
  //       break;
  //     }
  //
  //     mat[t][j] = vec[idx];
  //   }
  //   if (!valid) {
  //     for (int j = 0; j < E; ++j) {
  //       mat[t][j] = NaN;
  //     }
  //   }
  // }

  // Check which columns contain at least one non-NaN value
  std::vector<bool> keep(E, false);
  for (int j = 0; j < E; ++j) {
    for (size_t i = 0; i < N; ++i) {
      if (!std::isnan(mat[i][j])) {
        keep[j] = true;
        break;
      }
    }
  }

  // If no columns remain, throw exception
  bool any_column = false;
  for (bool b : keep) {
    if (b) {
      any_column = true;
      break;
    }
  }
  if (!any_column) {
    throw std::invalid_argument(
        "Embedding dimension E and lag tau are too large for input length, "
        "no valid embeddings can be generated."
    );
  }

  // Create cleaned matrix with only columns having valid data
  std::vector<std::vector<double>> cleaned;
  for (size_t i = 0; i < N; ++i) {
    std::vector<double> row;
    for (int j = 0; j < E; ++j) {
      if (keep[j]) {
        row.push_back(mat[i][j]);
      }
    }
    cleaned.push_back(row);
  }

  return cleaned;
}
