#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <utility>
#include <stdexcept>
#include "NumericUtils.h"
#include <RcppThread.h>

/**
 * @brief Compute L1 or L2 distance between two numeric vectors with NaN handling.
 *
 * @param vec1     First input vector
 * @param vec2     Second input vector
 * @param L1norm   If true → Manhattan (L1) distance, else → Euclidean (L2)
 * @param NA_rm    If true → ignore NaNs; if false → return NaN if any NaN present
 * @return double  Distance value (NaN if invalid or all-NaN)
 */
double CppDistance(const std::vector<double>& vec1,
                   const std::vector<double>& vec2,
                   bool L1norm = false,
                   bool NA_rm = false){
  const size_t n = vec1.size();
  // if (n != vec2.size()) throw std::invalid_argument("CppChebyshevDistance: Input vectors must have the same length.");

  double dist = 0.0;
  bool has_valid = false;

  if (L1norm) {
    // --- L1 (Manhattan) distance ---
    for (size_t i = 0; i < n; ++i) {
      double a = vec1[i], b = vec2[i];
      if (std::isnan(a) || std::isnan(b)) {
        if (!NA_rm) return std::numeric_limits<double>::quiet_NaN();
        else continue;
      }
      dist += std::abs(a - b);
      has_valid = true;
    }
  } else {
    // --- L2 (Euclidean) distance ---
    for (size_t i = 0; i < n; ++i) {
      double a = vec1[i], b = vec2[i];
      if (std::isnan(a) || std::isnan(b)) {
        if (!NA_rm) return std::numeric_limits<double>::quiet_NaN();
        else continue;
      }
      double diff = a - b;
      dist += diff * diff;
      has_valid = true;
    }
    dist = std::sqrt(dist);
  }

  return has_valid ? dist : std::numeric_limits<double>::quiet_NaN();
}

/**
 * @brief Compute the Chebyshev (L∞) distance between two numeric vectors.
 *
 * The Chebyshev distance is defined as:
 * \f[
 *    d(x, y) = \max_i |x_i - y_i|
 * \f]
 *
 * @param vec1   First numeric vector.
 * @param vec2   Second numeric vector (must have the same length as vec1).
 * @param NA_rm  Whether to remove NaN pairs before computing distance.
 *
 * @return Chebyshev distance (double). Returns NaN if no valid pairs or NA_rm=false with NaN present.
 */
double CppChebyshevDistance(const std::vector<double>& vec1,
                            const std::vector<double>& vec2,
                            bool NA_rm = false){
  // if (vec1.size() != vec2.size()) {
  //   throw std::invalid_argument("CppChebyshevDistance: Input vectors must have the same length.");
  // }

  double max_diff = 0.0;
  bool has_valid = false;

  const size_t n = vec1.size();
  for (size_t i = 0; i < n; ++i) {
    double a = vec1[i];
    double b = vec2[i];

    // Handle missing values
    if (std::isnan(a) || std::isnan(b)) {
      if (!NA_rm) return std::numeric_limits<double>::quiet_NaN();
      continue;
    }

    double diff = std::abs(a - b);
    if (diff > max_diff) max_diff = diff;
    has_valid = true;
  }

  // Return NaN if all pairs invalid
  if (!has_valid) return std::numeric_limits<double>::quiet_NaN();
  return max_diff;
}

// Function to compute the k-th nearest distance for a vector.
std::vector<double> CppKNearestDistance(const std::vector<double>& vec, size_t k,
                                        bool L1norm = false, bool NA_rm = false) {
  size_t n = vec.size();
  std::vector<double> result(n,std::numeric_limits<double>::quiet_NaN());  // Vector to store the k-th nearest distances

  for (size_t i = 0; i < n; ++i) {
    if (std::isnan(vec[i])) {
      continue;  // Skip if NA is encountered
    }

    std::vector<double> distances;
    distances.reserve(n);  // Reserve space to avoid repeated allocations

    for (size_t j = 0; j < n; ++j) {
      if (std::isnan(vec[j])) {
        if (!NA_rm) {
          distances.push_back(std::numeric_limits<double>::quiet_NaN());
          continue;  // Skip if NA is encountered and NA_rm is false
        } else {
          continue;  // Skip if NA is encountered and NA_rm is true
        }
      }

      double dist_res;
      if (L1norm) {
        dist_res = std::abs(vec[i] - vec[j]);  // Manhattan distance (L1)
      } else {
        double diff = vec[i] - vec[j];
        dist_res = diff * diff;  // Squared Euclidean distance (L2 squared)
      }
      distances.push_back(dist_res);
    }

    // Use nth_element to partially sort the distances up to the k-th element
    // This is more efficient than fully sorting the entire vector.
    if (k < distances.size()) {
      std::nth_element(distances.begin(), distances.begin() + k, distances.end());
      result[i] = distances[k];  // (k+1)-th smallest distance (exclude itself)
    } else {
      result[i] = *std::max_element(distances.begin(), distances.end());  // Handle case where k is out of bounds
    }

    // If using Euclidean distance, take the square root of the k-th nearest squared distance
    if (!L1norm) {
      result[i] = std::sqrt(result[i]);
    }
  }

  return result;
}

// Function to compute the k-th nearest Chebyshev distance for each sample in a matrix
std::vector<double> CppMatKNearestDistance(const std::vector<std::vector<double>>& mat,
                                           size_t k, bool NA_rm = false) {
  size_t n = mat.size();
  std::vector<double> result(n, std::numeric_limits<double>::quiet_NaN());

  for (size_t i = 0; i < n; ++i) {
    const auto& vec_i = mat[i];

    if (std::any_of(vec_i.begin(), vec_i.end(), [](double val) { return std::isnan(val); }) && !NA_rm) {
      continue;  // Skip if NA and NA_rm is false
    }

    std::vector<double> distances;
    distances.reserve(n - 1);

    for (size_t j = 0; j < n; ++j) {
      if (i == j) continue;

      double dist = CppChebyshevDistance(vec_i, mat[j], NA_rm);
      if (std::isnan(dist)) {
        if (!NA_rm) {
          distances.clear();
          break;
        } else {
          continue;
        }
      }

      distances.push_back(dist);
    }

    if (distances.empty()) continue;

    if (k < distances.size()) {
      std::nth_element(distances.begin(), distances.begin() + k, distances.end());
      result[i] = distances[k];
    } else {
      result[i] = *std::max_element(distances.begin(), distances.end());  // fallback if not enough neighbors
    }
  }

  return result;
}

// Function to compute distance for a matrix:
std::vector<std::vector<double>> CppMatDistance(
    const std::vector<std::vector<double>>& mat,
    bool L1norm = false,
    bool NA_rm = false){
  size_t n = mat.size();
  // std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n, std::numeric_limits<double>::quiet_NaN()));
  std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n, 0));

  // Compute distance between every pair of rows
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i+1; j < n; ++j) {  // <-- Corrected: increment j
      double distv = CppDistance(mat[i], mat[j], L1norm, NA_rm);
      distance_matrix[i][j] = distv;  // Correctly assign distance to upper triangle
      distance_matrix[j][i] = distv;  // Mirror the value to the lower triangle
      // distance_matrix[i][j] = distance_matrix[j][i] = CppDistance(mat[i],mat[j],L1norm,NA_rm);
    }
  }
  return distance_matrix;
}

// Function to compute chebyshev distance for a matrix:
std::vector<std::vector<double>> CppMatChebyshevDistance(
    const std::vector<std::vector<double>>& mat,
    bool NA_rm = false){
  size_t n = mat.size();
  // std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n, std::numeric_limits<double>::quiet_NaN()));
  std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n, 0));

  // Compute distance between every pair of rows
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i+1; j < n; ++j) {  // <-- Corrected: increment j
      double distv = CppChebyshevDistance(mat[i], mat[j], NA_rm);
      distance_matrix[i][j] = distv;  // Correctly assign distance to upper triangle
      distance_matrix[j][i] = distv;  // Mirror the value to the lower triangle
      // distance_matrix[i][j] = distance_matrix[j][i] = CppDistance(mat[i],mat[j],L1norm,NA_rm);
    }
  }
  return distance_matrix;
}

// Function to compute the number of neighbors for each point (in a vector) within a given radius.
std::vector<int> CppNeighborsNum(
    const std::vector<double>& vec,     // A vector of 1D points.
    const std::vector<double>& radius,  // A vector where radius[i] specifies the search radius for the i-th point.
    bool equal = false,                 // Flag to include points at exactly the radius distance (default: false).
    bool L1norm = false,                // Flag to use Manhattan distance or Euclidean distance
    bool NA_rm = false                  // Whether to remove the nan value in cpp
) {
  size_t N = vec.size();
  std::vector<int> NAx(N, 0); // Initialize neighbor counts to 0

  // Iterate over all pairs of points (i, j)
  for (size_t i = 0; i < N; ++i) {
    if (std::isnan(vec[i])) {
      continue;  // Skip if NA is encountered
    }

    for (size_t j = 0; j < N; ++j) {
      if (i != j) { // Skip self-comparison
        if (std::isnan(vec[j])) {
          continue;  // Skip if NA is encountered
        }

        double distance;
        if (L1norm) {
          distance = std::abs(vec[i] - vec[j]);  // Manhattan distance (L1)
        } else {
          double diff = vec[i] - vec[j];
          distance = std::sqrt(diff * diff);  // Euclidean distance (L2)
        }

        // Check neighbor condition based on the 'equal' flag
        if (!equal && distance < radius[i]) {
          NAx[i]++;
        } else if (equal && distance <= radius[i]) {
          NAx[i]++;
        }
      }
    }
  }

  return NAx;
}

// Function to compute the number of neighbors for each point (in a matrix) within a given radius.
// use the chebyshev distance
std::vector<int> CppMatNeighborsNum(
    const std::vector<std::vector<double>>& mat,     // A vector of 2D points.
    const std::vector<double>& radius,               // A vector where radius[i] specifies the search radius for the i-th point.
    bool equal = false,                              // Flag to include points at exactly the radius distance (default: false).
    bool NA_rm = false                               // Whether to remove the nan value in cpp
) {
  size_t N = mat.size();
  std::vector<int> NAx(N, 0); // Initialize neighbor counts to 0

  std::vector<std::vector<double>> dist = CppMatChebyshevDistance(mat,NA_rm);

  // Iterate over all pairs of points (i, j)
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (i != j) { // Skip self-comparison
        double distance = dist[i][j];
        // Check neighbor condition based on the 'equal' flag
        if (!equal && !doubleNearlyEqual(distance,radius[i]) && distance < radius[i]) {
          NAx[i]++;
        } else if (equal && doubleNearlyEqual(distance,radius[i]) && distance <= radius[i]) {
          NAx[i]++;
        }
      }
    }
  }

  return NAx;
}

// Function to find k-nearest neighbors of a given index in the embedding space
// The `lib` parameter specifies the indices from which the k-nearest neighbors should be selected
std::vector<size_t> CppKNNIndice(
    const std::vector<std::vector<double>>& embedding_space,  // Embedding space containing vectors
    size_t target_idx,                                        // Target index for which to find neighbors
    size_t k,                                                 // Number of nearest neighbors to find
    const std::vector<size_t>& lib,                           // Indices from which to select neighbors
    bool include_self = false)                                // Whether to include the point itself as a valid neighbor.
{
  std::vector<std::pair<double, size_t>> distances;

  // Iterate through the specified library indices to collect valid distances
  for (size_t i : lib) {
    if (!include_self && i == target_idx) continue;  // Skip the target index itself if include_self is false

    // Check if the entire embedding_space[i] is NaN
    if (std::all_of(embedding_space[i].begin(), embedding_space[i].end(),
                    [](double v) { return std::isnan(v); })) {
      continue;
    }

    // Compute the distance between the target and the current index
    double dist = CppDistance(embedding_space[target_idx], embedding_space[i], false, true);

    // Skip NaN distances
    if (!std::isnan(dist)) {
      distances.emplace_back(dist, i);
    }
  }

  // Partial sort to get k-nearest neighbors, excluding NaN distances
  std::partial_sort(distances.begin(), distances.begin() + std::min(k, distances.size()), distances.end());

  // Extract the indices of the k-nearest neighbors
  std::vector<std::size_t> neighbors;
  for (std::size_t i = 0; i < k && i < distances.size(); ++i) {
    neighbors.push_back(distances[i].second);
  }

  return neighbors;
}

// Function to find k-nearest neighbors of a given index using a precomputed distance matrix
// The `lib` parameter specifies the indices from which the k-nearest neighbors should be selected
std::vector<size_t> CppDistKNNIndice(
    const std::vector<std::vector<double>>& dist_mat,  // Precomputed n * n distance matrix
    size_t target_idx,                                 // Target index for which to find neighbors
    size_t k,                                          // Number of nearest neighbors to find
    const std::vector<size_t>& lib,                    // Indices from which to select neighbors
    bool include_self = false)                         // Whether to include the point itself as a valid neighbor.
{
  std::vector<std::pair<double, size_t>> distances;

  // Iterate through the specified library indices to collect valid distances
  for (size_t i : lib) {
    if (!include_self && i == target_idx) continue;  // Skip the target index itself if include_self is false

    double dist = dist_mat[target_idx][i];

    // Skip NaN distances
    if (!std::isnan(dist)) {
      distances.emplace_back(dist, i);
    }
  }

  // Partial sort to get k-nearest neighbors, excluding NaN distances
  std::partial_sort(distances.begin(), distances.begin() + std::min(k, distances.size()), distances.end());

  // Extract the indices of the k-nearest neighbors
  std::vector<std::size_t> neighbors;
  for (std::size_t i = 0; i < k && i < distances.size(); ++i) {
    neighbors.push_back(distances[i].second);
  }

  return neighbors;
}

/**
 * @brief Computes k-nearest neighbor indices for selected rows in a precomputed distance matrix.
 *
 * This function returns, for each specified row (in `lib`), the indices of its `k` nearest neighbors
 * (excluding or including self) among a given subset of indices (also `lib`). Distances must be
 * precomputed and stored in a full n x n matrix, which may include NaN values for invalid entries.
 *
 * @param dist_mat     Precomputed n x n distance matrix (may contain NaN for invalid distances).
 * @param lib          Subset of indices to restrict both the query points and candidate neighbors.
 * @param k            Number of neighbors to retain (if more than available, return all valid ones).
 * @param include_self Whether to include the point itself (i == j) as a valid neighbor.
 *
 * @return A vector of length n. Each element is a vector of up to k neighbor indices (from `lib`),
 *         sorted by increasing distance. For rows not in `lib`, a single entry with `invalid_index`.
 */
std::vector<std::vector<size_t>> CppDistSortedIndice(
    const std::vector<std::vector<double>>& dist_mat,
    const std::vector<size_t>& lib,
    size_t k,
    bool include_self = false){
  const size_t n = dist_mat.size();
  const size_t invalid_index = std::numeric_limits<size_t>::max();

  // Initialize result: fill each row with a single invalid index by default
  std::vector<std::vector<size_t>> sorted_indices(n, std::vector<size_t>{invalid_index});

  // Process each row specified in lib
  for (size_t i : lib) {
    if (i >= n || dist_mat[i].size() != n) continue;

    const auto& row = dist_mat[i];

    // Skip if self-distance is invalid (often indicates unusable row)
    if (std::isnan(row[i])) continue;

    std::vector<std::pair<double, size_t>> valid_neighbors;

    // Collect valid neighbors from lib
    for (size_t j : lib) {
      if (!include_self && i == j) continue; // Skip the target index itself if include_self is false
      double d = row[j];
      if (!std::isnan(d)) {
        valid_neighbors.emplace_back(d, j);
      }
    }

    // Take min(k, available neighbors)
    size_t num_neighbors = std::min(k, valid_neighbors.size());

    // Efficiently select top-k using partial_sort
    std::partial_sort(
      valid_neighbors.begin(),
      valid_neighbors.begin() + num_neighbors,
      valid_neighbors.end(),
      [](const std::pair<double, size_t>& a,
         const std::pair<double, size_t>& b) {
        if (!doubleNearlyEqual(a.first,b.first)) {
          return a.first < b.first;
        } else {
          return a.second < b.second;
        }
      }
    );

    // Extract indices of the top-k neighbors
    std::vector<size_t> indices;
    for (size_t m = 0; m < num_neighbors; ++m) {
      indices.push_back(valid_neighbors[m].second);
    }

    sorted_indices[i] = indices;
  }

  return sorted_indices;
}

/**
 * @brief Computes the k-nearest neighbors for a subset of points (lib) within the embedding space.
 *
 * For each index in 'lib', this function calculates the distance to all other indices in 'lib'
 * based on their corresponding vectors in 'embedding_space'. It then identifies the k closest
 * neighbors (excluding the point itself) and returns a vector of neighbors for each point.
 *
 * The returned vector has the same number of rows as embedding_space; rows not in 'lib' are
 * initialized with a vector containing a single invalid index.
 *
 * This function uses RcppThread to parallelize the outer loop over the 'lib' indices,
 * with the number of threads controlled by the 'threads' parameter.
 *
 * @param embedding_space - The full set of vectors representing the embedding space.
 * @param lib - A vector of indices representing the subset of points to consider for neighbor search.
 * @param k - The number of nearest neighbors to find for each point in 'lib'.
 * @param threads - The number of threads to use for parallel computation.
 * @param L1norm - Flag to use Manhattan distance (true) or Euclidean distance (false).
 * @param include_self - Whether to include the point itself (i == j) as a valid neighbor.
 * @return std::vector<std::vector<size_t>> A vector where each row corresponds to an embedding_space
 *         point and contains the indices of its k nearest neighbors from 'lib', or an invalid index if not in 'lib'.
 */
std::vector<std::vector<size_t>> CppMatKNNeighbors(
    const std::vector<std::vector<double>>& embedding_space,
    const std::vector<size_t>& lib,
    size_t k,
    size_t threads,
    bool L1norm = false,
    bool include_self = false) {

  const size_t n = embedding_space.size();
  const size_t invalid_index = std::numeric_limits<size_t>::max();

  // Initialize the result vector with a single invalid index per row
  std::vector<std::vector<size_t>> sorted_indices(n, std::vector<size_t>{invalid_index});

  const size_t lib_size = lib.size();

  // // Iterate over each element in lib
  // for (size_t idx_i = 0; idx_i < lib_size; ++idx_i) {
  //   size_t i = lib[idx_i];
  //
  //   // Vector to store pairs of (distance, neighbor_index)
  //   std::vector<std::pair<double, size_t>> dist_idx_pairs;
  //
  //   // Compute distances to all other points in lib (excluding self)
  //   for (size_t idx_j = 0; idx_j < lib_size; ++idx_j) {
  //     if (!include_self && idx_i == idx_j) continue; // Skip the target index itself if include_self is false
  //     size_t j = lib[idx_j];
  //
  //     double dist = CppDistance(embedding_space[i], embedding_space[j], L1norm, true);
  //
  //     if (std::isnan(dist)) continue; // Skip invalid distances
  //
  //     dist_idx_pairs.emplace_back(dist, j);
  //   }
  //
  //   // Determine the number of neighbors to keep (min(k, available neighbors))
  //   size_t knn = std::min(k, dist_idx_pairs.size());
  //
  //   // Partially sort distances to find the k smallest distances efficiently
  //   std::partial_sort(dist_idx_pairs.begin(), dist_idx_pairs.begin() + knn, dist_idx_pairs.end(),
  //                    [](const std::pair<double, size_t>& a,
  //                       const std::pair<double, size_t>& b) {
  //                      if (!doubleNearlyEqual(a.first,b.first)) {
  //                        return a.first < b.first;
  //                      } else {
  //                       return a.second < b.second;
  //                      }
  //                    });
  //
  //   // Resize the output vector for the current point and fill with neighbor indices
  //   sorted_indices[i].resize(knn);
  //   for (size_t m = 0; m < knn; ++m) {
  //     sorted_indices[i][m] = dist_idx_pairs[m].second;
  //   }
  // }

  // Parallel loop over lib indices using RcppThread
  RcppThread::parallelFor(0, lib_size, [&](size_t idx_i) {
    size_t i = lib[idx_i];

    // Vector to store pairs of (distance, neighbor_index)
    std::vector<std::pair<double, size_t>> dist_idx_pairs;

    // Compute distances to all other points in lib (excluding self)
    for (size_t idx_j = 0; idx_j < lib_size; ++idx_j) {
      if (!include_self && idx_i == idx_j) continue; // Skip the target index itself if include_self is false
      size_t j = lib[idx_j];

      double dist = CppDistance(embedding_space[i], embedding_space[j], L1norm, true);

      if (std::isnan(dist)) continue; // Skip invalid distances

      dist_idx_pairs.emplace_back(dist, j);
    }

    // Determine the number of neighbors to keep (min(k, available neighbors))
    size_t knn = std::min(k, dist_idx_pairs.size());

    // Partially sort distances to find the k smallest distances efficiently
    std::partial_sort(dist_idx_pairs.begin(), dist_idx_pairs.begin() + knn, dist_idx_pairs.end(),
                      [](const std::pair<double, size_t>& a,
                         const std::pair<double, size_t>& b) {
                        if (!doubleNearlyEqual(a.first,b.first)) {
                          return a.first < b.first;
                        } else {
                          return a.second < b.second;
                        }
                      });

    // Resize the output vector for the current point and fill with neighbor indices
    sorted_indices[i].resize(knn);
    for (size_t m = 0; m < knn; ++m) {
      sorted_indices[i][m] = dist_idx_pairs[m].second;
    }
  }, threads);

  return sorted_indices;
}
