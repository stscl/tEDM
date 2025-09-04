#ifndef MultispatialCCM_H
#define MultispatialCCM_H

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
 *     [3] = 2.5% percentile (lower bound of 95% CI)
 *     [4] = 97.5% percentile (upper bound of 95% CI)
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
);

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
 *     [3] = 2.5% percentile of rho distribution (lower CI bound)
 *     [4] = 97.5% percentile of rho distribution (upper CI bound)
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
);

#endif // CCM_H
