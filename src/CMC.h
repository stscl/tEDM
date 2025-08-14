#ifndef CMC_H
#define CMC_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <utility>
#include <unordered_set>
#include "CppStats.h"
#include "CppDistances.h"
#include "tEDMDataStruct.h"
#include "IntersectionCardinality.h"
#include <RcppThread.h>

/**
 * @brief Computes the Cross Mapping Cardinality (CMC) causal strength score.
 *
 * This function evaluates the directional causal influence from one time series
 * to another using Cross Mapping Cardinality. It performs state-space reconstruction,
 * neighbor searching, and statistical evaluation across a range of library sizes.
 *
 * @param embedding_x State-space reconstructed time series of the potential cause.
 * @param embedding_y State-space reconstructed time series of the potential effect.
 * @param lib_sizes A vector of library sizes to use for subsampling during CMC analysis.
 * @param lib Indices of library points (0-based).
 * @param pred Indices of prediction points (0-based).
 * @param num_neighbors Number of neighbors used in cross mapping.
 * @param n_excluded Number of temporally excluded neighbors (Theiler window).
 * @param dist_metric Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param threads Number of threads for parallel processing.
 * @param parallel_level Level of parallelism to control nested parallel execution.
 * @param progressbar Boolean flag to show or hide a progress bar.
 *
 * @return CMCRes A struct containing:
 *   - cross_mapping: A vector of AUC values for the largest library size.
 *   - causal_strength: A 2D vector with rows [library size, mean AUC] across all lib sizes.
 */
CMCRes CMC(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<size_t>& lib_sizes,
    const std::vector<size_t>& lib,
    const std::vector<size_t>& pred,
    size_t num_neighbors = 4,
    size_t n_excluded = 0,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0,
    bool progressbar = true);

#endif // CMC_H
