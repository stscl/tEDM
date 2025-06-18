#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "Embed.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

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
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
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
                                            int threads) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Unique sorted embedding dimensions and neighbor values
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  // Generate unique (E, b) combinations
  std::vector<std::pair<int, int>> unique_Ebcom;
  for (int e : Es)
    for (int bb : bs)
      unique_Ebcom.emplace_back(e, bb);

  std::vector<std::vector<double>> result(unique_Ebcom.size(), std::vector<double>(5));

  RcppThread::parallelFor(0, unique_Ebcom.size(), [&](size_t i) {
    const int Ei = unique_Ebcom[i].first;
    const int bi = unique_Ebcom[i].second;

    auto embeddings = Embed(source, Ei, tau);
    auto metrics = SimplexBehavior(embeddings, target, lib_indices, pred_indices, bi);

    result[i][0] = Ei;
    result[i][1] = bi;
    result[i][2] = metrics[0]; // rho
    result[i][3] = metrics[1]; // MAE
    result[i][4] = metrics[2]; // RMSE
  }, threads_sizet);

  return result;
}

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
 *   - tau: The spatial lag step for constructing lagged state-space vectors.
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
                                         int threads){
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Generate embeddings
  auto embeddings = Embed(source, E, tau);
  std::vector<std::vector<double>> result(theta.size(), std::vector<double>(4));

  RcppThread::parallelFor(0, theta.size(), [&](size_t i) {
    auto metrics = SMapBehavior(embeddings, target, lib_indices, pred_indices, b, theta[i]);

    result[i][0] = theta[i];   // theta
    result[i][1] = metrics[0]; // rho
    result[i][2] = metrics[1]; // MAE
    result[i][3] = metrics[2]; // RMSE
  }, threads_sizet);

  return result;
}
