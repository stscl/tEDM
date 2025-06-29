#ifndef FNN_H
#define FNN_H

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "CppStats.h"
#include <RcppThread.h>

/*
 * Compute the False Nearest Neighbors (FNN) ratio for spatial cross-sectional data.
 *
 * This function determines whether nearest neighbors identified in a lower-dimensional
 * embedded space (E1) remain close in a higher-dimensional space (E2).
 * If not, the neighbor is considered a "false" neighbor, indicating the need for
 * a higher embedding dimension to accurately capture spatial proximity.
 *
 * Parameters:
 * - embedding: A matrix (vector of vectors) representing the spatial embedding,
 *              where each row corresponds to a spatial unit's attributes.
 *              Must contain at least E2 columns.
 * - lib: Library index vector (1-based in R, converted to 0-based).
 * - pred: Prediction index vector (1-based in R, converted to 0-based).
 * - E1: The base embedding dimension used to identify the nearest neighbor (E1 < E2).
 * - E2: The full embedding dimension used to test false neighbors (usually E1 + 1).
 * - Rtol: Relative threshold (default 10.0). If the change in the added dimension is
 *         large relative to the E1 distance, the neighbor is considered false.
 * - Atol: Absolute threshold (default 2.0). If the added dimension changes too much
 *         absolutely, the neighbor is also considered false.
 * - L1norm: Whether to use Manhattan (L1) distance instead of Euclidean (L2).
 *
 * Returns:
 * - A double value indicating the proportion of false nearest neighbors (0–1).
 *   If no valid pairs are found, returns NaN.
 */
double CppSingleFNN(const std::vector<std::vector<double>>& embedding,
                    const std::vector<int>& lib,
                    const std::vector<int>& pred,
                    size_t E1,
                    size_t E2,
                    double Rtol = 10.0,
                    double Atol = 2.0,
                    bool L1norm = false);

/*
 * Compute False Nearest Neighbor (FNN) ratios for a range of embedding dimensions
 * on spatial cross-sectional data.
 *
 * This function iteratively calls CppSingleFNN using pairs of embedding dimensions:
 * from 1 to D-1 (where D is the number of columns in the embedding matrix), such that
 * each E1 = d and E2 = d + 1.
 *
 * It is used to identify the optimal embedding dimension by observing how the proportion
 * of false nearest neighbors decreases with higher dimensions.
 *
 * Parameters:
 * - embedding: A matrix (vector of vectors) where each row is a spatial unit's
 *              multidimensional embedding. Should have at least 2 columns.
 * - lib: Library index vector (1-based in R, converted to 0-based).
 * - pred: Prediction index vector (1-based in R, converted to 0-based).
 * - Rtol: Vectors of relative distance threshold.
 * - Atol: Vectors of absolute distance threshold.
 * - L1norm: Whether to use L1 (Manhattan) distance instead of L2 (Euclidean).
 * - threads: Number of parallel threads.
 *
 * Returns:
 * - A vector of doubles, where each entry corresponds to the FNN ratio for a given
 *   E1 (from 1 to D - 1), with E2 = E1 + 1. If the result is invalid, NaN is returned.
 *   Returns an NaN vector if the embedding has fewer than 2 columns.
 */
std::vector<double> CppFNN(const std::vector<std::vector<double>>& embedding,
                           const std::vector<int>& lib,
                           const std::vector<int>& pred,
                           const std::vector<double>& Rtol,
                           const std::vector<double>& Atol,
                           bool L1norm, int threads);

#endif // FNN_H
