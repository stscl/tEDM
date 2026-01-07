#ifndef IntersectionCardinality_H
#define IntersectionCardinality_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include "CppStats.h"
#include "CppDistances.h"
#include "DataStruct.h"
#include <RcppThread.h>

/**
 * Computes intersection-based mapping ratio sequences between two neighbor graphs
 * for use in Cross Mapping Cardinality (CMC) or similar causal inference frameworks.
 *
 * Parameters:
 *   neighborsX     - Precomputed sorted neighbor indices for embedding X
 *   neighborsY     - Precomputed sorted neighbor indices for embedding Y
 *   lib_size       - Size of the moving library used in mapping
 *   lib_indices    - Global indices from which to draw the sliding libraries
 *   pred_indices   - Indices at which to perform prediction (evaluation points)
 *   num_neighbors  - Number of neighbors used for mapping (after exclusion)
 *   n_excluded     - Number of nearest neighbors to exclude from the front
 *   threads        - Number of parallel threads for computation
 *   parallel_level - Whether to use multithreaded (0) or serial (1) mode
 *
 * Returns:
 *   A vector of IntersectionRes structures, each containing the average intersection
 *   ratio sequence (IC curve) for a different starting point of the moving library.
 *   If lib_size == lib_indices.size(), returns a single result using full library.
 *
 * Notes:
 *   - Neighbor lists must use std::numeric_limits<size_t>::max() to indicate invalid entries.
 *   - This function assumes that the neighbor vectors are sorted by ascending distance.
 *   - Use in combination with AUC computation to assess causal strength.
 */
std::vector<IntersectionRes> IntersectionCardinalitySingle(
    const std::vector<std::vector<size_t>>& neighborsX,
    const std::vector<std::vector<size_t>>& neighborsY,
    size_t lib_size,
    const std::vector<size_t>& lib_indices,
    const std::vector<size_t>& pred_indices,
    size_t num_neighbors,
    size_t n_excluded,
    size_t threads,
    int parallel_level
);

/**
 * Computes the Intersection Cardinality (IC) curve for causal inference via cross mapping.
 *
 * This function evaluates the extent to which neighbors of the effect variable Y
 * are preserved when mapped through the neighbors of cause variable X.
 * Specifically, for each number of neighbors from 1 to `num_neighbors`, it computes
 * the intersection count between the k nearest neighbors of Y and the k nearest neighbors of X,
 * for each prediction point.
 *
 * The output is an Intersection Cardinality (IC) curve, which can be further processed
 * (e.g., calculating AUC, statistical significance) outside this function.
 *
 * @param embedding_x     State-space embedding of the potential cause variable (NxE matrix).
 * @param embedding_y     State-space embedding of the potential effect variable (NxE matrix).
 * @param lib             Vector of library indices (shouble be 0-based in C++).
 * @param pred            Vector of prediction indices (shouble be 0-based in C++).
 * @param num_neighbors   Maximum number of neighbors to consider in intersection (e.g., from 1 to k).
 * @param n_excluded      Number of nearest neighbors to exclude (e.g., due to temporal proximity).
 * @param dist_metric     Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param threads         Number of threads used for parallel computation.
 * @param parallel_level  Parallel mode flag: 0 = parallel, 1 = serial.
 *
 * @return A vector of size `num_neighbors`:
 *         - Each element represents the average number of overlapping neighbors
 *           between X and Y across prediction points, for each neighbor count k = 1, 2, ..., num_neighbors.
 *
 *         If inputs are invalid or no valid prediction points exist, the returned vector
 *         is filled with `NaN` values.
 *
 * @note
 *   - This function returns only the raw intersection values. To compute AUC or p-values,
 *     use additional post-processing such as DeLong’s test.
 */
std::vector<double> IntersectionCardinality(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<size_t>& lib,
    const std::vector<size_t>& pred,
    size_t num_neighbors,
    size_t n_excluded,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0);

/**
 * Computes the Intersection Cardinality (IC) AUC-based causal strength score.
 *
 * This function evaluates the extent to which neighbors of the effect variable Y
 * are preserved when mapped through the neighbors of cause variable X, by calculating
 * the intersection ratio curve and evaluating its AUC (area under curve).
 * Statistical significance (p-value) and confidence interval are computed via DeLong's test.
 *
 * Parameters:
 *   embedding_x    - State-space reconstruction (embedding) of the potential cause variable.
 *   embedding_y    - State-space reconstruction (embedding) of the potential effect variable.
 *   lib            - Library index vector (shouble be 0-based in C++).
 *   pred           - Prediction index vector (shouble be 0-based in C++).
 *   num_neighbors  - Number of neighbors used for cross mapping (after exclusion).
 *   n_excluded     - Number of nearest neighbors to exclude (e.g. temporal).
 *   dist_metric    - Distance metric selector (1: Manhattan, 2: Euclidean).
 *   threads        - Number of threads used in parallel computation.
 *   parallel_level - Whether to use multithreaded (0) or serial (1) mode
 *
 * Returns:
 *   A vector of 4 values:
 *     [0] - AUC (Intersection Cardinality score, bounded [0, 1])
 *     [1] - p-value from DeLong test (testing whether AUC > 0.5)
 *     [2] - Confidence interval lower bound
 *     [3] - Confidence interval upper bound
 */
std::vector<double> IntersectionCardinalityScores(
    const std::vector<std::vector<double>>& embedding_x,
    const std::vector<std::vector<double>>& embedding_y,
    const std::vector<size_t>& lib,
    const std::vector<size_t>& pred,
    size_t num_neighbors,
    size_t n_excluded,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0);

#endif // IntersectionCardinality_H
