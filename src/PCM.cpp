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
#include "tEDMDataStruct.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

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
    bool cumulate
){
  int n_controls = controls.size();
  std::vector<double> rho(2, std::numeric_limits<double>::quiet_NaN());

  if (cumulate) {
    // Cumulative control embedding path
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0]);
      } else {
        temp_pred = SimplexProjectionPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors[i]);
      }
      temp_embedding = Embed(temp_pred,conEs[i],taus[i]);
    }

    std::vector<double> con_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[n_controls]);
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0]);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3){
      rho[0] = PearsonCor(target,target_pred,true);
      rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
    }
  } else {
    // Independent control paths
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SimplexProjectionPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0]);
      temp_embedding = Embed(temp_pred,conEs[i],taus[i]);
      temp_pred = SimplexProjectionPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[i+1]);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SimplexProjectionPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0]);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3){
      rho[0] = PearsonCor(target,target_pred,true);
      rho[1] = PartialCor(target,target_pred,con_pred,true,false);
    }
  }

  return rho;
}

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
    double theta,
    bool cumulate
){
  int n_controls = controls.size();
  std::vector<double> rho(2, std::numeric_limits<double>::quiet_NaN());

  if (cumulate){
    // Cumulative control embedding path
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      if (i == 0){
        temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0], theta);
      } else {
        temp_pred = SMapPrediction(temp_embedding, controls[i], lib_indices, pred_indices, num_neighbors[i], theta);
      }
      temp_embedding = Embed(temp_pred,conEs[i],taus[i]);
    }

    std::vector<double> con_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[n_controls], theta);
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0], theta);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3){
      rho[0] = PearsonCor(target,target_pred,true);
      rho[1] = PartialCorTrivar(target,target_pred,con_pred,true,false);
    }
  } else {
    // Independent control paths
    std::vector<std::vector<double>> con_pred(n_controls);
    std::vector<double> temp_pred;
    std::vector<std::vector<double>> temp_embedding;

    for (int i = 0; i < n_controls; ++i) {
      temp_pred = SMapPrediction(vectors, controls[i], lib_indices, pred_indices, num_neighbors[0], theta);
      temp_embedding = Embed(temp_pred,conEs[i],taus[i]);
      temp_pred = SMapPrediction(temp_embedding, target, lib_indices, pred_indices, num_neighbors[i+1], theta);
      con_pred[i] = temp_pred;
    }
    std::vector<double> target_pred = SMapPrediction(vectors, target, lib_indices, pred_indices, num_neighbors[0], theta);

    if (checkOneDimVectorNotNanNum(target_pred) >= 3){
      rho[0] = PearsonCor(target,target_pred,true);
      rho[1] = PartialCor(target,target_pred,con_pred,true,false);
    }
  }

  return rho;
}

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
    bool cumulate                                       // Whether to cumulate the partial correlations
) {
  int max_lib_size = lib_indices.size();

  // No possible library variation if using all vectors
  if (lib_size == max_lib_size) {
    std::vector<PartialCorRes> x_xmap_y;
    // Run partial cross map and store results
    std::vector<double> rho;
    if (simplex) {
      rho = PartialSimplex4TS(x_vectors, y, controls, lib_indices, pred_indices, conEs, taus, b, cumulate);
    } else {
      rho = PartialSMap4TS(x_vectors, y, controls, lib_indices, pred_indices, conEs, taus, b, theta, cumulate);
    }
    x_xmap_y.emplace_back(lib_size, rho[0], rho[1]);
    return x_xmap_y;
  } else if (parallel_level == 0){
    // Precompute valid indices for the library
    std::vector<std::vector<int>> valid_lib_indices;
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Loop around to beginning of lib indices
      if (start_lib + lib_size > max_lib_size) {
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        // // no wrapping around
        // int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        // for (int i = 0; i < num_vectors_remaining; ++i) {
        //   local_lib_indices.emplace_back(lib_indices[i]);
        // }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }
      valid_lib_indices.emplace_back(local_lib_indices);
    }

    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<PartialCorRes> x_xmap_y(valid_lib_indices.size());

    // Perform the operations using RcppThread
    RcppThread::parallelFor(0, valid_lib_indices.size(), [&](size_t i) {
      // Run partial cross map and store results
      std::vector<double> rho;
      if (simplex) {
        rho = PartialSimplex4TS(x_vectors, y, controls, valid_lib_indices[i], pred_indices, conEs, taus, b, cumulate);
      } else {
        rho = PartialSMap4TS(x_vectors, y, controls, valid_lib_indices[i], pred_indices, conEs, taus, b, theta, cumulate);
      }
      // Directly initialize a PartialCorRes struct with the three values
      PartialCorRes result(lib_size, rho[0], rho[1]);
      x_xmap_y[i] = result;
    }, threads);

    return x_xmap_y;
  } else {
    std::vector<PartialCorRes> x_xmap_y;

    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Setup changing library
      if (start_lib + lib_size > max_lib_size) { // Loop around to beginning of lib indices
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        // // no wrapping around
        // int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        // for (int i = 0; i < num_vectors_remaining; ++i) {
        //   local_lib_indices.emplace_back(lib_indices[i]);
        // }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }

      // Run partial cross map and store results
      std::vector<double> rho;
      if (simplex) {
        rho = PartialSimplex4TS(x_vectors, y, controls, local_lib_indices, pred_indices, conEs, taus, b, cumulate);
      } else {
        rho = PartialSMap4TS(x_vectors, y, controls, local_lib_indices, pred_indices, conEs, taus, b, theta, cumulate);
      }
      x_xmap_y.emplace_back(lib_size, rho[0], rho[1]);
    }

    return x_xmap_y;
  }
}

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
 * - progressbar: Boolean flag indicating whether to display a progress bar during computation.
 *
 * Returns:
 *    A 2D vector of results, where each row contains:
 *      - The library size.
 *      - The mean pearson cross-mapping correlation.
 *      - The statistical significance of the pearson correlation.
 *      - The upper bound of the pearson correlation confidence interval.
 *      - The lower bound of the pearson correlation confidence interval.
 *      - The mean partial cross-mapping partial correlation.
 *      - The statistical significance of the partial correlation.
 *      - The upper bound of the partial correlation confidence interval.
 *      - The lower bound of the partial correlation confidence interval.
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
    bool progressbar                                    // Whether to print the progress bar
) {
  // If b is not provided correctly, default it to E + 2
  std::vector<int> bs = b;
  for (size_t i = 0; i < bs.size(); ++i){
    if (bs[i] <= 0) {
      bs[i] = Es[i] + 2;
    }
  }

  int Ex = Es[0];
  std::vector<int> conEs = Es;
  conEs.erase(conEs.begin());

  int taux = taus[0];
  std::vector<int> contaus = taus;
  contaus.erase(contaus.begin());

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  std::vector<std::vector<double>> x_vectors = Embed(x,Ex,taux);
  size_t n = pred.size();

  size_t n_confounds;
  if (cumulate){
    n_confounds = 1;
  } else {
    n_confounds = controls.size();
  }

  int max_lib_size = static_cast<int>(lib.size()); // Maximum lib size

  std::vector<int> unique_lib_sizes(lib_sizes.begin(), lib_sizes.end());

  // Transform to ensure no size exceeds max_lib_size
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size); });

  // Ensure the minimum value in unique_lib_sizes is Ex + 2 (uncomment this section if required)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](int size) { return std::max(size, Ex + 2); });

  // Remove duplicates
  std::sort(unique_lib_sizes.begin(), unique_lib_sizes.end());
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Local results for each library
  std::vector<std::vector<PartialCorRes>> local_results(unique_lib_sizes.size());

  if (parallel_level == 0){
    // Iterate over each library size
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        local_results[i] = PCMSingle(
          x_vectors,
          y,
          controls,
          unique_lib_sizes[i],
          lib,
          pred,
          conEs,
          contaus,
          bs,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate
        );
        bar++;
      }
    } else {
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        local_results[i] = PCMSingle(
          x_vectors,
          y,
          controls,
          unique_lib_sizes[i],
          lib,
          pred,
          conEs,
          contaus,
          bs,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate
        );
      }
    }
  } else {
    // Perform the operations using RcppThread
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        local_results[i] = PCMSingle(
          x_vectors,
          y,
          controls,
          lib_size,
          lib,
          pred,
          conEs,
          contaus,
          bs,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate
        );
        bar++;
      }, threads_sizet);
    } else {
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        local_results[i] = PCMSingle(
          x_vectors,
          y,
          controls,
          lib_size,
          lib,
          pred,
          conEs,
          contaus,
          bs,
          simplex,
          theta,
          threads_sizet,
          parallel_level,
          cumulate
        );
      }, threads_sizet);
    }
  }

  // Initialize the result container
  std::vector<PartialCorRes> x_xmap_y;

  // Merge all local results into the final result
  for (const auto& local_result : local_results) {
    x_xmap_y.insert(x_xmap_y.end(), local_result.begin(), local_result.end());
  }

  // Group by the first int and store second and third values as pairs
  std::map<int, std::vector<std::pair<double, double>>> grouped_results;

  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].emplace_back(result.second, result.third);
  }

  // // Previous implementation calculated significance and confidence intervals using the mean of rho vector only.
  // // This approach is now deprecated and kept here for comparison purposes.
  // std::vector<std::vector<double>> final_results;
  // // Compute the mean of second and third values for each group
  // for (const auto& group : grouped_results) {
  //   std::vector<double> second_values, third_values;
  //
  //   for (const auto& val : group.second) {
  //     second_values.push_back(val.first);
  //     third_values.push_back(val.second);
  //   }
  //
  //   double mean_second = CppMean(second_values, true);
  //   double mean_third = CppMean(third_values, true);
  //
  //   final_results.push_back({static_cast<double>(group.first), mean_second, mean_third});
  // }
  //
  // // Compute significance and confidence intervals for each result
  // for (size_t i = 0; i < final_results.size(); ++i) {
  //   double rho_second = final_results[i][1];
  //   double rho_third = final_results[i][2];
  //
  //   // Compute significance and confidence interval for second value
  //   double significance_second = CppCorSignificance(rho_second, n);
  //   std::vector<double> confidence_interval_second = CppCorConfidence(rho_second, n);
  //
  //   // Compute significance and confidence interval for third value
  //   double significance_third = CppCorSignificance(rho_third, n, n_confounds);
  //   std::vector<double> confidence_interval_third = CppCorConfidence(rho_third, n, n_confounds);
  //
  //   // Append computed statistical values to the result
  //   final_results[i].push_back(significance_second);
  //   final_results[i].push_back(confidence_interval_second[0]);
  //   final_results[i].push_back(confidence_interval_second[1]);
  //
  //   final_results[i].push_back(significance_third);
  //   final_results[i].push_back(confidence_interval_third[0]);
  //   final_results[i].push_back(confidence_interval_third[1]);
  // }

  // For each group, compute the mean of second and third values and calculate significance and confidence intervals using the original vectors
  std::vector<std::vector<double>> final_results;
  for (const auto& group : grouped_results) {
    std::vector<double> second_values, third_values;

    // Collect all second and third values from current group
    for (const auto& val : group.second) {
      second_values.push_back(val.first);
      third_values.push_back(val.second);
    }

    // Compute mean values for reporting
    double mean_second = CppMean(second_values, true);
    double mean_third = CppMean(third_values, true);

    // Compute significance and confidence intervals using the full vector of values (not just mean)
    double significance_second = CppMeanCorSignificance(second_values, n);
    std::vector<double> confidence_interval_second = CppMeanCorConfidence(second_values, n);

    double significance_third = CppMeanCorSignificance(third_values, n, n_confounds);
    std::vector<double> confidence_interval_third = CppMeanCorConfidence(third_values, n, n_confounds);

    // Store group ID, mean values, and corresponding statistical results
    final_results.push_back({
      static_cast<double>(group.first),
      mean_second,
      mean_third,
      significance_second,
      confidence_interval_second[0],
      confidence_interval_second[1],
      significance_third,
      confidence_interval_third[0],
      confidence_interval_third[1]
    });
  }

  return final_results;
}
