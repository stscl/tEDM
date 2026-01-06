#include <vector>
#include <cstdint>
#include <algorithm>
#include <utility>
#include <string>
#include "DataStruct.h"

/**
 * @brief Computes DeLong placements for ROC analysis.
 *
 * This function implements the DeLong method to calculate placements for cases and controls
 * in the context of ROC (Receiver Operating Characteristic) curve analysis. The method is
 * used for non-parametric estimation of the area under the curve (AUC) and supports handling
 * tied values efficiently.
 *
 * The function accepts two sets of data points: `cases` and `controls`. Depending on the
 * specified `direction` (either ">" or "<"), it processes the data accordingly by inverting
 * the sign of the input values if necessary.
 *
 * @reference https://github.com/xrobin/pROC/blob/master/src/delong.cpp
 *
 * @param cases A vector of numeric values representing positive cases.
 * @param controls A vector of numeric values representing negative controls.
 * @param direction A string indicating the comparison direction. If set to ">", values are inverted.
 *
 * @return A structure containing:
 *   - theta: The estimated AUC value.
 *   - X: A vector of normalized placement values for cases.
 *   - Y: A vector of normalized placement values for controls.
 */
DeLongPlacementsRes CppDeLongPlacements(const std::vector<double>& cases,
                                        const std::vector<double>& controls,
                                        const std::string& direction) {
  // Initialize variables
  size_t m = cases.size();
  size_t n = controls.size();
  size_t L = m + n;

  // Create working copies to handle direction
  std::vector<double> proc_cases = cases;
  std::vector<double> proc_controls = controls;

  // Handle direction parameter
  if (direction == ">") {
    for (auto& val : proc_cases) val = -val;
    for (auto& val : proc_controls) val = -val;
  }

  // Create combined vector with indices and class labels
  // Use size_t instead of int for indices and Replace vector<bool> with vector<uint8_t>
  std::vector<std::pair<size_t, double>> Z;
  std::vector<uint8_t> labels(L, 0);

  // Populate case data
  for (size_t i = 0; i < m; ++i) {
    Z.emplace_back(i, proc_cases[i]);
    labels[i] = 1;
  }

  // Populate control data
  for (size_t j = 0; j < n; ++j) {
    Z.emplace_back(m + j, proc_controls[j]);
    labels[m + j] = 0;
  }

  // Sort combined data by value
  std::sort(Z.begin(), Z.end(), [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
    return a.second < b.second;
  });

  // Calculate XY placements
  std::vector<double> XY(L, 0.0);
  size_t current_m = 0, current_n = 0;
  size_t i = 0;

  while (i < L) {
    std::vector<size_t> case_indices, control_indices;
    size_t case_count = 0, control_count = 0;

    // Handle tied values
    while (true) {
      size_t index = Z[i].first;
      if (labels[index]) {
        ++case_count;
        case_indices.push_back(index);
      } else {
        ++control_count;
        control_indices.push_back(index);
      }

      if (i == L-1 || Z[i].second != Z[i+1].second) break;
      ++i;
    }

    // Update XY values for cases
    for (size_t idx : case_indices) {
      XY[idx] = current_n + control_count/2.0;
    }

    // Update XY values for controls
    for (size_t idx : control_indices) {
      XY[idx] = current_m + case_count/2.0;
    }

    // Accumulate counts
    current_m += case_count;
    current_n += control_count;
    ++i;
  }

  // Calculate final X and Y vectors
  DeLongPlacementsRes ret;
  // ret.X.reserve(n);
  // ret.Y.reserve(m);

  double sum = 0.0;
  const double norm_n = static_cast<double>(n);
  const double norm_m = static_cast<double>(m);

  for (size_t k = 0; k < L; ++k) {
    if (labels[k]) {
      sum += XY[k];
      ret.X.push_back(XY[k] / norm_n);
    } else {
      ret.Y.push_back(1.0 - XY[k] / norm_m);
    }
  }

  // Calculate theta (AUC estimate)
  ret.theta = sum / (norm_m * norm_n);

  return ret;
}
