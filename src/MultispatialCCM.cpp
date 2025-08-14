#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include "CppStats.h"
#include "Embed.h"
#include "SimplexProjection.h"
#include <RcppThread.h>

// [[Rcpp::depends(RcppThread)]]

/*
 * Perform bootstrapped simplex prediction from spatially replicated short time series.
 * Supports parallel execution using RcppThread if parallel_level == 0.
 *
 * Parameters:
 *   source         - A vector of vectors representing explanatory variables (plots × time).
 *   target         - A vector of vectors representing response variables (plots × time).
 *   E              - Embedding dimension.
 *   tau            - Time delay.
 *   lib            - Number of plots to sample in each bootstrap replicate.
 *   num_neighbors  - Number of nearest neighbors to use in simplex projection.
 *   threads        - Number of nearest neighbors to use in bootstrap replicates.
 *   seed           - Random seed for reproducibility.
 *   boot           - Number of bootstrap replicates.
 *   parallel_level - If 0, run in parallel using RcppThread. Otherwise, run sequentially.
 *   dist_metric    - Distance metric selector (1: Manhattan, 2: Euclidean).
 *   dist_average   - Whether to average distance by the number of valid vector components.
 *
 * Returns:
 *   A vector of length 5:
 *     [0] = library size (libsize)
 *     [1] = mean of Pearson correlation (mean rho)
 *     [2] = one-sided p-value (proportion of rho ≤ 0)
 *     [3] = 97.5% percentile (upper bound of 95% CI)
 *     [4] = 2.5% percentile (lower bound of 95% CI)
 */
std::vector<double> SimplexPredictionBoot(
    const std::vector<std::vector<double>>& source,
    const std::vector<std::vector<double>>& target,
    int E,
    int tau,
    int lib,
    int num_neighbors,
    int boot,
    size_t threads,
    unsigned int seed = 42,
    int parallel_level = 0,
    int dist_metric = 2,
    bool dist_average = true
) {
  int n_plot = source.size();
  std::vector<double> rho_list(boot, std::numeric_limits<double>::quiet_NaN());

  // Prebuild RNG pool with seed sequence
  std::vector<std::mt19937> rng_pool(boot);
  for (int i = 0; i < boot; ++i) {
    std::seed_seq seq{static_cast<uint32_t>(seed), static_cast<uint32_t>(i)};
    rng_pool[i] = std::mt19937(seq);
  }

  auto run_boot_once = [&](int b) {
    std::mt19937& rng = rng_pool[b];
    std::uniform_int_distribution<> plot_sampler(0, n_plot - 1);

    // 1. Bootstrap plots
    std::vector<int> lib_plots(lib);
    for (int i = 0; i < lib; ++i) {
      lib_plots[i] = plot_sampler(rng);
    }

    // 2. Embedding
    std::vector<std::vector<double>> library_vectors;
    std::vector<double> library_targets;

    for (int plot : lib_plots) {
      auto emb = Embed(source[plot], E, tau);
      int T = target[plot].size();
      for (int i = 0; i < static_cast<int>(emb.size()); ++i) {
        int ti = i + tau * (E - 1);
        if (ti + 1 < T) {
          library_vectors.push_back(emb[i]);
          library_targets.push_back(target[plot][ti + 1]);
        }
      }
    }

    int N = library_vectors.size();
    if (N < num_neighbors) {
      return;
    }

    std::vector<int> all_indices(N);
    std::iota(all_indices.begin(), all_indices.end(), 0);

    auto pred = SimplexProjectionPrediction(library_vectors, library_targets,
                                            all_indices, all_indices, num_neighbors,
                                            dist_metric, dist_average);

    rho_list[b] = PearsonCor(library_targets, pred, true);
  };

  // Parallel or sequential execution
  if (parallel_level == 0) {
    RcppThread::parallelFor(0, boot, run_boot_once, threads);
  } else {
    for (int b = 0; b < boot; ++b) run_boot_once(b);
  }

  // 3. Summary statistics
  double mean_rho = CppMean(rho_list, true);

  int n_non_na = 0, count_le_0 = 0;
  std::vector<double> clean_rho;
  for (double r : rho_list) {
    if (!std::isnan(r)) {
      ++n_non_na;
      clean_rho.push_back(r);
      if (r <= 0) ++count_le_0;
    }
  }

  double pval = (n_non_na > 0) ? (1.0 + count_le_0) / (1.0 + n_non_na)
    : std::numeric_limits<double>::quiet_NaN();

  std::sort(clean_rho.begin(), clean_rho.end());
  double ci_lower = std::numeric_limits<double>::quiet_NaN();
  double ci_upper = std::numeric_limits<double>::quiet_NaN();
  int n = clean_rho.size();
  if (n >= 1) {
    ci_lower = clean_rho[std::clamp(int(std::floor(0.025 * n)), 0, n - 1)];
    ci_upper = clean_rho[std::clamp(int(std::ceil(0.975 * n)) - 1, 0, n - 1)];
  }

  return {static_cast<double>(lib), mean_rho, pval, ci_upper, ci_lower};
}

/*
 * Conduct spatial convergent cross mapping (multispatialCCM) with bootstrapped simplex projection
 * across multiple library sizes using spatial replicates of short time series.
 *
 * The method evaluates how forecast skill (e.g., Pearson's ρ between observed and predicted) improves
 * with increasing library size (number of spatial replicates), which serves as a signal of causality
 * in the empirical dynamic modeling (EDM) framework. The function runs bootstrapped simplex projection
 * over a user-defined set of library sizes and aggregates the results to produce mean forecast skill,
 * p-values, and confidence intervals.
 *
 * Supports parallel execution via RcppThread with optional progress monitoring.
 *
 * Parameters:
 *   x              - A vector of vectors representing explanatory variables (plots × time).
 *   y              - A vector of vectors representing response variables (plots × time).
 *   lib_sizes      - A list of library sizes to evaluate (number of plots to sample).
 *   E              - Embedding dimension for state space reconstruction.
 *   tau            - Time delay between lags in the embedding.
 *   b              - Number of nearest neighbors used in simplex projection (defaults to E + 1).
 *   boot           - Number of bootstrap replicates for each library size.
 *   threads        - Number of threads to use for parallel processing.
 *   seed           - Random seed for reproducibility.
 *   parallel_level - 0 for sequential execution; >0 enables parallelization with RcppThread.
 *   dist_metric    - Distance metric selector (1: Manhattan, 2: Euclidean).
 *   dist_average   - Whether to average distance by the number of valid vector components.
 *   progressbar    - Logical flag to enable/disable a progress bar during execution.
 *
 * Returns:
 *   A vector of vectors, where each inner vector corresponds to a library size and contains:
 *     [0] = library size (libsize)
 *     [1] = mean of Pearson correlation (mean rho)
 *     [2] = one-sided p-value (proportion of rho ≤ 0 across bootstraps)
 *     [3] = 97.5% percentile of rho distribution (upper CI bound)
 *     [4] = 2.5% percentile of rho distribution (lower CI bound)
 */
std::vector<std::vector<double>> MultispatialCCM(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& y,
    const std::vector<int>& lib_sizes,
    int E,
    int tau,
    int b,
    int boot,
    int threads,
    unsigned int seed = 42,
    int parallel_level = 0,
    int dist_metric = 2,
    bool dist_average = true,
    bool progressbar = true
) {
  // If b is not provided correctly, default it to E + 1
  if (b <= 0) {
    b = E + 1;
  }

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Transform to ensure no size exceeds max library size
  int max_lib_size = x.size();

  std::vector<int> unique_lib_sizes(lib_sizes.begin(), lib_sizes.end());
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size); });

  // Ensure the minimum value in unique_lib_sizes is E + 1 (uncomment this section if required)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](int size) { return std::max(size, E + 1); });

  // Remove duplicates
  std::sort(unique_lib_sizes.begin(), unique_lib_sizes.end());
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Local results for each library
  std::vector<std::vector<double>> local_results(unique_lib_sizes.size());

  if (parallel_level == 0){
    // Iterate over each library size
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        local_results[i] = SimplexPredictionBoot(
          x,
          y,
          E,
          tau,
          unique_lib_sizes[i],
          b,
          boot,
          threads_sizet,
          seed,
          parallel_level,
          dist_metric,
          dist_average);
        bar++;
      }
    } else {
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        local_results[i] = SimplexPredictionBoot(
          x,
          y,
          E,
          tau,
          unique_lib_sizes[i],
          b,
          boot,
          threads_sizet,
          seed,
          parallel_level,
          dist_metric,
          dist_average);
      }
    }
  } else {
    // Perform the operations using RcppThread
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        local_results[i] = SimplexPredictionBoot(
          x,
          y,
          E,
          tau,
          lib_size,
          b,
          boot,
          threads_sizet,
          seed,
          parallel_level,
          dist_metric,
          dist_average);
        bar++;
      }, threads_sizet);
    } else {
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        local_results[i] = SimplexPredictionBoot(
          x,
          y,
          E,
          tau,
          lib_size,
          b,
          boot,
          threads_sizet,
          seed,
          parallel_level,
          dist_metric,
          dist_average);
      }, threads_sizet);
    }
  }

  return local_results;
}
