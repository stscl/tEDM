#ifndef FNN_H
#define FNN_H

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include "NumericUtils.h"
#include "CppStats.h"
#include "CppDistances.h"

// Note: <RcppThread.h> is intentionally excluded from this header to avoid
//       unnecessary Rcpp dependencies and potential header inclusion order
//       issues (e.g., R.h being included before Rcpp headers). It should only
//       be included in the corresponding .cpp implementation file.
// #include <RcppThread.h>

/*
 * Compute the False Nearest Neighbors (FNN) ratio for time series data.
 *
 * This function determines whether nearest neighbors identified in a lower-dimensional
 * embedded space (E1) remain close in a higher-dimensional space (E2).
 * If not, the neighbor is considered a "false" neighbor, indicating the need for
 * a higher embedding dimension to accurately capture spatial proximity.
 *
 * The FNN test is computed in two modes:
 * - parallel_level = 0: each prediction is processed in parallel using RcppThreads.
 * - parallel_level = 1: all pairwise distances are precomputed once in advance
 *   (better for repeated queries or small prediction sets).
 *
 * Parameters:
 * - embedding: A matrix (vector of vectors) representing the spatial embedding,
 *              where each row corresponds to a spatial unit's attributes.
 *              Must contain at least E2 columns.
 * - lib: Library index vector (1-based in R, converted to 0-based).
 * - pred: Prediction index vector (1-based in R, converted to 0-based).
 * - E1: The base embedding dimension used to identify the nearest neighbor (E1 < E2).
 * - E2: The full embedding dimension used to test false neighbors (usually E1 + 1).
 * - threads: Number of threads used when parallel_level = 0.
 * - parallel_level: 0 for per-pred parallelism (default), 1 for precomputed full distance matrix.
 * - Rtol: Relative threshold (default 10.0).
 * - Atol: Absolute threshold (default 2.0).
 * - L1norm: Whether to use Manhattan (L1) distance instead of Euclidean (L2).
 *
 * Returns:
 * - A double value indicating the proportion of false nearest neighbors (0–1).
 *   If no valid pairs are found, returns NaN.
 */
double CppSingleFNN(const std::vector<std::vector<double>>& embedding,
                    const std::vector<size_t>& lib,
                    const std::vector<size_t>& pred,
                    size_t E1,
                    size_t E2,
                    size_t threads,
                    int parallel_level = 0,
                    double Rtol = 10.0,
                    double Atol = 2.0,
                    bool L1norm = false);

/*
 * Compute False Nearest Neighbor (FNN) ratios across multiple embedding dimensions
 * for time series data.
 *
 * For a given embedding matrix (with each row representing a spatial unit and
 * each column an embedding dimension), this function evaluates the proportion
 * of false nearest neighbors (FNN) as the embedding dimension increases.
 *
 * It iteratively calls `CppSingleFNN` for each embedding dimension pair (E1, E2),
 * where E1 ranges from 1 to D - 1 (D = number of columns), and E2 = E1 + 1.
 * The FNN ratio measures how often a nearest neighbor in dimension E1 becomes
 * distant in dimension E2, suggesting that E1 is insufficient for reconstructing
 * the system.
 *
 * If `parallel_level == 0`, the function executes in serial;
 * otherwise, it uses multithreading to compute FNN ratios for each (E1, E2) pair
 * in parallel.
 *
 * Parameters:
 * - embedding: A vector of vectors where each row is a spatial unit’s embedding.
 *              Must have at least 2 columns (dimensions).
 * - lib: A vector of integer indices indicating the library set (0-based).
 * - pred: A vector of integer indices indicating the prediction set (0-based).
 * - Rtol: A vector of relative distance thresholds (one per E1).
 * - Atol: A vector of absolute distance thresholds (one per E1).
 * - L1norm: If true, use L1 (Manhattan) distance; otherwise, use L2 (Euclidean).
 * - threads: Number of threads to use for parallel computation.
 * - parallel_level: 0 for serial loop over E1, >0 for parallel loop over E1.
 *
 * Returns:
 * - A vector of FNN ratios corresponding to each E1 from 1 to D - 1.
 *   If not computable for a given E1, NaN is returned at that position.
 */
std::vector<double> CppFNN(const std::vector<std::vector<double>>& embedding,
                           const std::vector<size_t>& lib,
                           const std::vector<size_t>& pred,
                           const std::vector<double>& Rtol,
                           const std::vector<double>& Atol,
                           bool L1norm = false,
                           int threads = 8,
                           int parallel_level = 0);

#endif // FNN_H
