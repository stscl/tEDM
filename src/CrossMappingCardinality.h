#ifndef CrossMappingCardinality_H
#define CrossMappingCardinality_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <utility>
#include <unordered_set>
#include "CppStats.h"
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
 * @param threads Number of threads for parallel processing.
 * @param parallel_level Level of parallelism to control nested parallel execution.
 * @param progressbar Boolean flag to show or hide a progress bar.
 *
 * @return CMCRes A struct containing:
 *   - cross_mapping: A vector of AUC values for the largest library size.
 *   - causal_strength: A 2D vector with rows [library size, mean AUC] across all lib sizes.
 */
CMCRes CrossMappingCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<size_t>& lib_sizes,
    const std::vector<size_t>& lib,
    const std::vector<size_t>& pred,
    size_t num_neighbors = 4,
    size_t n_excluded = 0,
    int threads = 8,
    int parallel_level = 0,
    bool progressbar = true);

// /**
//  * Computes the Cross Mapping Cardinality (CMC) causal strength score (adjusted based on Python logic).
//  *
//  * Parameters:
//  *   embedding_x: State-space reconstruction (embedded) of the potential cause variable.
//  *   embedding_y: State-space reconstruction (embedded) of the potential effect variable.
//  *   lib: Library index vector (1-based in R, converted to 0-based).
//  *   pred: Prediction index vector (1-based in R, converted to 0-based).
//  *   num_neighbors: Vector of numbers of neighbors used for cross mapping (corresponding to n_neighbor in python package crossmapy).
//  *   n_excluded: Vector of numbers of neighbors excluded from the distance matrix (corresponding to n_excluded in python package crossmapy).
//  *   threads: Number of parallel threads.
//  *   parallel_level: the level of parallelization
//  *   progressbar: Whether to display a progress bar.
//  *
//  * Returns:
//  *   A vector the results of the DeLong test for the AUC values: [number of neighbors, IC score, p-value, confidence interval upper bound, confidence interval lower bound] one for each entry in num_neighbors.
//  *   The result contains multiple rows, each corresponding to a different number of neighbors.
//  */
// std::vector<std::vector<double>> CrossMappingCardinality(
//     const std::vector<std::vector<double>>& embedding_x,
//     const std::vector<std::vector<double>>& embedding_y,
//     const std::vector<int>& lib,
//     const std::vector<int>& pred,
//     const std::vector<int>& num_neighbors,
//     const std::vector<int>& n_excluded,
//     int threads,
//     int parallel_level = 0,
//     bool progressbar = true);

#endif // CrossMappingCardinality_H
