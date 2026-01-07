#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include "CppStats.h"
#include "CppDistances.h"
#include "Embed.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionalCardinality.h"
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
                                            int threads = 8) {
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
    auto metrics = SimplexBehavior(embeddings, target, lib_indices, pred_indices, bi, dist_metric, dist_average);

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
                                         int threads = 8){
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Generate embeddings
  auto embeddings = Embed(source, E, tau);
  std::vector<std::vector<double>> result(theta.size(), std::vector<double>(4));

  RcppThread::parallelFor(0, theta.size(), [&](size_t i) {
    auto metrics = SMapBehavior(embeddings, target, lib_indices, pred_indices, b, theta[i], dist_metric, dist_average);

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
                                                 int threads = 8) {
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

  RcppThread::parallelFor(0, unique_Ebcom.size(), [&](size_t idx) {
    const int Ei = unique_Ebcom[idx].first;
    const int bi = unique_Ebcom[idx].second;

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
    size_t max_lag = static_cast<size_t>((tau == 0) ? (Ei - 1) : (Ei * tau));
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

    std::vector<double> target_pred = SimplexProjectionPrediction(all_vectors, all_targets, com_lib, com_pred, bi, dist_metric, dist_average);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3) {
      pearson = PearsonCor(target_pred, all_targets, true);
      mae = CppMAE(target_pred, all_targets, true);
      rmse = CppRMSE(target_pred, all_targets, true);
    }

    result[idx][0] = Ei;
    result[idx][1] = bi;
    result[idx][2] = pearson; // rho
    result[idx][3] = mae;     // MAE
    result[idx][4] = rmse;    // RMSE
  }, threads_sizet);

  return result;
}

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
                                       int parallel_level = 0) {
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

  std::vector<std::vector<double>> result(unique_Ebcom.size(), std::vector<double>(4));

  size_t max_num_neighbors = 0;
  if (!bs.empty()) {
    max_num_neighbors = static_cast<size_t>(bs.back() + exclude);
  }

  if (parallel_level == 0){
    for (size_t i = 0; i < Es.size(); ++i) {
      // Generate embeddings
      auto embedding_x = Embed(source, Es[i], tau);
      auto embedding_y = Embed(target, Es[i], tau);

      // Filter valid prediction points (exclude those with all NaN values)
      std::vector<size_t> valid_pred;
      for (size_t idx : pred_indices) {
        if (idx < 0 || idx >= embedding_x.size()) continue;

        bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                                 [](double v) { return std::isnan(v); });
        bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                                 [](double v) { return std::isnan(v); });
        if (!x_nan && !y_nan) valid_pred.push_back(idx);
      }

      // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
      bool L1norm = (dist_metric == 1);

      // // Precompute neighbors (The earlier implementation based on a serial version)
      // auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, L1norm, true),lib_indices,max_num_neighbors);
      // auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, L1norm, true),lib_indices,max_num_neighbors);

      // Precompute neighbors (parallel computation)
      auto nx = CppMatKNNeighbors(embedding_x, lib_indices, max_num_neighbors, threads_sizet, L1norm);
      auto ny = CppMatKNNeighbors(embedding_y, lib_indices, max_num_neighbors, threads_sizet, L1norm);

      // Parameter initialization
      const size_t n_excluded_sizet = static_cast<size_t>(exclude);

      for (size_t j = 0; j < bs.size(); ++j){
        const size_t k = static_cast<size_t>(bs[j]);

        // run cross mapping
        std::vector<IntersectionRes> res = IntersectionalCardinalitySingle(
          nx,ny,lib_indices.size(),lib_indices,valid_pred,k,n_excluded_sizet,threads_sizet,0
        );

        std::vector<double> cs = {0,1};
        if (!res.empty())  cs = CppCMCTest(res[0].Intersection,">");

        result[j + bs.size() * i][0] = Es[i];  // E
        result[j + bs.size() * i][1] = bs[j];  // k
        result[j + bs.size() * i][2] = cs[0];  // AUC
        result[j + bs.size() * i][3] = cs[1];  // P value
      }
    }
  } else {
    for (size_t i = 0; i < Es.size(); ++i) {
      // Generate embeddings
      auto embedding_x = Embed(source, Es[i], tau);
      auto embedding_y = Embed(target, Es[i], tau);

      // Filter valid prediction points (exclude those with all NaN values)
      std::vector<size_t> valid_pred;
      for (size_t idx : pred_indices) {
        if (idx < 0 || idx >= embedding_x.size()) continue;

        bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                                 [](double v) { return std::isnan(v); });
        bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                                 [](double v) { return std::isnan(v); });
        if (!x_nan && !y_nan) valid_pred.push_back(idx);
      }

      // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
      bool L1norm = (dist_metric == 1);

      // // Precompute neighbors (The earlier implementation based on a serial version)
      // auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, L1norm, true),lib_indices,max_num_neighbors);
      // auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, L1norm, true),lib_indices,max_num_neighbors);

      // Precompute neighbors (parallel computation)
      auto nx = CppMatKNNeighbors(embedding_x, lib_indices, max_num_neighbors, threads_sizet, L1norm);
      auto ny = CppMatKNNeighbors(embedding_y, lib_indices, max_num_neighbors, threads_sizet, L1norm);

      // Parameter initialization
      const size_t n_excluded_sizet = static_cast<size_t>(exclude);

      RcppThread::parallelFor(0, bs.size(), [&](size_t j) {
        const size_t k = static_cast<size_t>(bs[j]);

        // run cross mapping
        std::vector<IntersectionRes> res = IntersectionalCardinalitySingle(
          nx,ny,lib_indices.size(),lib_indices,valid_pred,k,n_excluded_sizet,threads_sizet,1
        );

        std::vector<double> cs = {0,1};
        if (!res.empty())  cs = CppCMCTest(res[0].Intersection,">");

        result[j + bs.size() * i][0] = Es[i];  // E
        result[j + bs.size() * i][1] = bs[j];  // k
        result[j + bs.size() * i][2] = cs[0];  // AUC
        result[j + bs.size() * i][3] = cs[1];  // P value
      }, threads_sizet);
    }
  }

  return result;
}
