#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "CppStats.h"
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
                                                 int threads) {
  // Sort and deduplicate lib and pred
  std::vector<int> lib;
  lib.reserve(lib_indices.size());
  lib.insert(lib.end(), lib_indices.begin(), lib_indices.end());
  std::sort(lib.begin(), lib.end());
  lib.erase(std::unique(lib.begin(), lib.end()), lib.end());

  std::vector<int> pred;
  pred.reserve(pred_indices.size());
  pred.insert(pred.end(), pred_indices.begin(), pred_indices.end());
  std::sort(pred.begin(), pred.end());
  pred.erase(std::unique(pred.begin(), pred.end()), pred.end());

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

    // Combine all spatial observations and generate embeddings
    std::vector<std::vector<double>> all_vectors;
    std::vector<double> all_targets;

    for (size_t i = 0; i < source.size(); ++i) {
      auto emb = Embed(source[i], Ei, tau);
      all_vectors.insert(all_vectors.end(), emb.begin(), emb.end());
      all_targets.insert(all_targets.end(), target[i].begin(), target[i].end());
    }

    // Construct lib_indices and pred_indices
    std::vector<int> com_lib;
    std::vector<int> com_pred;

    int ts_len = source[0].size();
    int max_lag = (tau == 0) ? (Ei - 1) : (Ei * tau);
    for (size_t i = 0; i < lib.size(); ++i) {
      for (size_t j = 0; j < source[0].size(); ++j) {
        if (j > max_lag){
          com_lib.push_back(j + lib[i] * ts_len);
        }
      }
    }
    for (size_t i = 0; i < pred.size(); ++i) {
      for (size_t j = 0; j < source[0].size(); ++j) {
        if (j > max_lag){
          com_pred.push_back(j + pred[i] * ts_len);
        }
      }
    }

    double pearson = std::numeric_limits<double>::quiet_NaN();
    double mae = std::numeric_limits<double>::quiet_NaN();
    double rmse = std::numeric_limits<double>::quiet_NaN();

    std::vector<double> target_pred = SimplexProjectionPrediction(all_vectors, all_targets, com_lib, com_pred, bi);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
      pearson = PearsonCor(target_pred, all_targets, true);
      mae = CppMAE(target_pred, all_targets, true);
      rmse = CppRMSE(target_pred, all_targets, true);
    }

    result[i][0] = Ei;
    result[i][1] = bi;
    result[i][2] = pearson; // rho
    result[i][3] = mae; // MAE
    result[i][4] = rmse; // RMSE
  }, threads_sizet);

  return result;
}
