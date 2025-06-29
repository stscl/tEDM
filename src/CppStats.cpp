#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric> // for std::accumulate
#include <limits>  // for std::numeric_limits
#include "DeLongPlacements.h"
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Function to check if a value is NA
bool isNA(double value) {
  return std::isnan(value);
}

// Function to check if a indice of int type is NA
bool checkIntNA(int value) {
  return value == std::numeric_limits<int>::min();
}

// Function to check whether a one-dimensional vector contains NaN values
bool checkOneDimVectorHasNaN(const std::vector<double>& vec) {
  for (double val : vec) {
    if (std::isnan(val)) {
      return true; // Return false if any NaN is found
    }
  }
  return false;
}

// Function to count the number of non-NaN values in a one-dimensional vector
int checkOneDimVectorNotNanNum(const std::vector<double>& vec) {
  int count = 0; // Initialize the counter for non-NaN values
  for (double val : vec) {
    if (!std::isnan(val)) {
      ++count; // Increment the counter if the value is not NaN
    }
  }
  return count; // Return the count of non-NaN values
}

/**
 * Computes the factorial of a non-negative integer using an iterative approach.
 *
 * @param n The non-negative integer for which factorial is to be computed.
 * @return The factorial of n as unsigned long long. Returns 0 if n > 20 due to overflow.
 *
 * @note Efficiently computes factorial in O(n) time with O(1) space complexity.
 *       Maximum representable value is 20! (2432902008176640000) due to 64-bit limits.
 */
unsigned long long CppFactorial(unsigned int n) {
  if (n > 20) {
    return 0;  // Overflow occurs beyond 20! for 64-bit unsigned integers
  }
  unsigned long long result = 1;
  for (unsigned int i = 2; i <= n; ++i) {
    result *= i;
  }
  return result;
}

/**
 * Computes the binomial coefficient C(n, k) using multiplicative decomposition.
 *
 * @param n The total number of items.
 * @param k The number of items to choose.
 * @return The binomial coefficient C(n, k) as unsigned long long. Returns 0 if k > n.
 *
 * @note Optimizes by selecting smaller k (or n-k) and using stepwise multiply-divide
 *       operations to reduce intermediate overflows. Time complexity O(k), space O(1).
 *       Result accuracy depends on final value fitting within 64-bit limits.
 */
unsigned long long CppCombine(unsigned int n, unsigned int k) {
  if (k > n) {
    return 0;  // Invalid case when selecting more items than available
  }
  // Use smaller k to minimize iterations (C(n, k) == C(n, n-k))
  k = (k < n - k) ? k : n - k;
  if (k == 0) {
    return 1;  // Directly return 1 when no selection is required
  }

  unsigned long long result = 1;
  for (unsigned int i = 1; i <= k; ++i) {
    // Multiplicative formula: result = result * (n - k + i) / i
    result *= (n - k + i);  // Accumulate numerator term
    result /= i;            // Divide by denominator term stepwise
  }
  return result;
}

/**
 * @brief Computes the digamma (ψ) function approximation for a given input x.
 *
 * This function provides an approximation of the digamma function, which is
 * the first derivative of the logarithm of the gamma function. The implementation
 * follows an asymptotic expansion for large values of x and applies a recurrence
 * relation to shift smaller values upward.
 *
 * The approach consists of:
 * - A loop to increase x if it is less than or equal to 5, adjusting the result accordingly.
 * - An asymptotic series expansion for the digamma function when x is sufficiently large.
 * - The final result is computed as a sum of the logarithm term, correction terms,
 *   and accumulated adjustments from the recurrence relation.
 *
 * @param x The input value for which the digamma function is evaluated.
 * @return The approximate digamma function value ψ(x).
 */
double CppDigamma(double x) {
  double a = 0;
  while (x <= 5) {
    a -= 1 / x;
    x += 1;
  }

  double b = 1 / (x * x);
  double c = b * (-1/12.0 +
             b * (1/120.0 +
             b * (-1/252.0 +
             b * (1/240.0 +
             b * (-1/132.0 +
             b * (691/32760.0 +
             b * (-1/12.0 +
             b * 3617/8160.0)))))));

  return (a + std::log(x) - 0.5 / x + c);
}

// Function to calculate the logarithm of x with the specified base
double CppLog(double x, double base = 10) {
  return std::log(x) / std::log(base);
}

// Function to calculate the median of a vector.
double CppMedian(const std::vector<double>& vec, bool NA_rm = false) {
  std::vector<double> filtered_vec;
  for (const double& value : vec) {
    if (std::isnan(value)) {
      if (!NA_rm) {
        // Return NaN immediately if not removing NaNs
        return std::numeric_limits<double>::quiet_NaN();
      }
      // else skip the NaN
    } else {
      filtered_vec.push_back(value);
    }
  }

  // If no valid values remain, return NaN
  if (filtered_vec.empty()) return std::numeric_limits<double>::quiet_NaN();

  // Sort the filtered vector
  std::sort(filtered_vec.begin(), filtered_vec.end());

  // Compute the median
  size_t n = filtered_vec.size();
  if (n % 2 == 0) {
    // If even, return the average of the two middle values
    return (filtered_vec[n / 2 - 1] + filtered_vec[n / 2]) / 2.0;
  } else {
    // If odd, return the middle value
    return filtered_vec[n / 2];
  }
}

// Function to calculate the mean of a vector, ignoring NA values
double CppMean(const std::vector<double>& vec, bool NA_rm = false) {
  double sum = 0.0;
  size_t count = 0;
  for (const auto& value : vec) {
    if (!NA_rm || !isNA(value)) {
      sum += value;
      ++count;
    }
  }
  return count > 0 ? sum / count : std::numeric_limits<double>::quiet_NaN();
}

// Function to calculate the minimum of a vector, ignoring NA values if NA_rm is true
double CppMin(const std::vector<double>& vec, bool NA_rm = false) {
  double min_val = std::numeric_limits<double>::infinity();
  bool found_valid = false;

  for (const auto& value : vec) {
    if (isNA(value)) {
      if (!NA_rm) {
        // Return NaN immediately if NA_rm is false and we hit a NaN
        return std::numeric_limits<double>::quiet_NaN();
      }
      continue;  // skip if NA_rm is true
    }

    if (!found_valid || value < min_val) {
      min_val = value;
      found_valid = true;
    }
  }

  // Return NaN if no valid value was found
  return found_valid ? min_val : std::numeric_limits<double>::quiet_NaN();
}


// Function to calculate the maximum of a vector, ignoring NA values if NA_rm is true
double CppMax(const std::vector<double>& vec, bool NA_rm = false) {
  double max_val = -std::numeric_limits<double>::infinity();
  bool found_valid = false;

  for (const auto& value : vec) {
    if (isNA(value)) {
      if (!NA_rm) return std::numeric_limits<double>::quiet_NaN();
      continue;
    }

    if (!found_valid || value > max_val) {
      max_val = value;
      found_valid = true;
    }
  }

  // Return NaN if no valid value was found
  return found_valid ? max_val : std::numeric_limits<double>::quiet_NaN();
}

// Function to calculate the sum of a vector, ignoring NA values if NA_rm is true
double CppSum(const std::vector<double>& vec,
              bool NA_rm = false) {
  double sum = 0.0;
  for (const auto& value : vec) {
    if (!NA_rm || !isNA(value)) {
      sum += value;
    }
  }
  return sum;
}

// Function to compute Mean Absolute Error (MAE) between two vectors
double CppMAE(const std::vector<double>& vec1,
              const std::vector<double>& vec2,
              bool NA_rm = false) {
  // Check if input vectors have the same size
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  // Initialize variables for MAE calculation
  double sum_abs_diff = 0.0; // Sum of absolute differences
  size_t valid_count = 0;    // Count of valid (non-NaN) pairs

  // Iterate through the vectors
  for (size_t i = 0; i < vec1.size(); ++i) {
    // Check if either vec1[i] or vec2[i] is NaN
    if (isNA(vec1[i]) || isNA(vec2[i])) {
      if (!NA_rm) {
        // If NA_rm is false and NaN is encountered, return NaN
        return std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      // If both values are valid, compute absolute difference and add to sum
      sum_abs_diff += std::fabs(vec1[i] - vec2[i]);
      valid_count++;
    }
  }

  // If no valid pairs are found, return NaN
  if (valid_count == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Compute and return MAE
  return sum_abs_diff / static_cast<double>(valid_count);
}

// Function to compute Root Mean Squared Error (RMSE) between two vectors
double CppRMSE(const std::vector<double>& vec1,
               const std::vector<double>& vec2,
               bool NA_rm = false) {
  // Check if input vectors have the same size
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  // Initialize variables for RMSE calculation
  double sum_squared_diff = 0.0; // Sum of squared differences
  size_t valid_count = 0;        // Count of valid (non-NaN) pairs

  // Iterate through the vectors
  for (size_t i = 0; i < vec1.size(); ++i) {
    // Check if either vec1[i] or vec2[i] is NaN
    if (isNA(vec1[i]) || isNA(vec2[i])) {
      if (!NA_rm) {
        // If NA_rm is false and NaN is encountered, return NaN
        return std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      // If both values are valid, compute squared difference and add to sum
      double diff = vec1[i] - vec2[i];
      sum_squared_diff += std::pow(diff, 2);
      valid_count++;
    }
  }

  // If no valid pairs are found, return NaN
  if (valid_count == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Compute and return RMSE
  return std::sqrt(sum_squared_diff / static_cast<double>(valid_count));
}

// Function to calculate the cumulative sum of a vector of doubles
// Returns a new vector where each element at index i is the sum of elements from index 0 to i in the input vector
std::vector<double> CppCumSum(const std::vector<double>& vec) {
  // Create a vector to hold the cumulative sums with the same size as the input vector
  std::vector<double> cumSum(vec.size());

  // Initialize the first element of the cumulative sum
  if (!vec.empty()) {
    cumSum[0] = vec[0];
  }

  // Calculate the cumulative sum
  for (size_t i = 1; i < vec.size(); ++i) {
    cumSum[i] = cumSum[i - 1] + vec[i]; // Add the current element to the previous cumulative sum
  }

  return cumSum; // Return the vector containing cumulative sums
}

// Function to calculate the absolute difference between two vectors
std::vector<double> CppAbsDiff(const std::vector<double>& vec1,
                               const std::vector<double>& vec2) {
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  std::vector<double> result(vec1.size());
  for (size_t i = 0; i < vec1.size(); ++i) {
    result[i] = std::abs(vec1[i] - vec2[i]);
  }
  return result;
}

// Function to normalize a vector by dividing each element by the sum of all elements
std::vector<double> CppSumNormalize(const std::vector<double>& vec,
                                    bool NA_rm = false) {
  double sum = CppSum(vec, NA_rm);
  if (sum == 0.0) {
    throw std::invalid_argument("Sum of vector elements is zero, cannot normalize.");
  }

  std::vector<double> normalizedVec(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    if (!isNA(vec[i])) {
      normalizedVec[i] = vec[i] / sum;
    } else {
      normalizedVec[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  return normalizedVec;
}

// Generates an arithmetic sequence of numbers starting from `from` and ending at `to`,
// with a total of `length_out` elements. The sequence is evenly spaced.
std::vector<double> CppArithmeticSeq(double from, double to, size_t length_out) {
  // Check for invalid input
  if (length_out < 1) {
    throw std::invalid_argument("length_out must be at least 1.");
  }

  // Initialize the result vector
  std::vector<double> res;
  res.reserve(length_out); // Reserve space for efficiency

  // If length_out is 1, return a vector containing only `from`
  if (length_out == 1) {
    res.push_back(from);
    return res;
  }

  // Calculate the step size for evenly spaced elements
  double step = (to - from) / (length_out - 1);

  // Generate the sequence
  for (size_t i = 0; i < length_out; ++i) {
    res.push_back(from + i * step);
  }

  return res;
}

// Function to calculate the variance of a vector, ignoring NA values
double CppVariance(const std::vector<double>& vec, bool NA_rm = false) {
  double mean_val = CppMean(vec, NA_rm);
  double var = 0.0;
  size_t count = 0;
  for (const auto& value : vec) {
    if (!NA_rm || !isNA(value)) {
      var += (value - mean_val) * (value - mean_val);
      ++count;
    }
  }
  return count > 1 ? var / (count - 1) : std::numeric_limits<double>::quiet_NaN();
}

// Function to calculate the covariance of two vectors, ignoring NA values
double CppCovariance(const std::vector<double>& vec1,
                     const std::vector<double>& vec2,
                     bool NA_rm = false) {
  if (vec1.size() != vec2.size()) {
    throw std::invalid_argument("Vectors must have the same size");
  }

  double mean1 = CppMean(vec1, NA_rm);
  double mean2 = CppMean(vec2, NA_rm);
  double cov = 0.0;
  size_t count = 0;
  for (size_t i = 0; i < vec1.size(); ++i) {
    if ((!NA_rm || !isNA(vec1[i])) && (!NA_rm || !isNA(vec2[i]))) {
      cov += (vec1[i] - mean1) * (vec2[i] - mean2);
      ++count;
    }
  }
  return count > 1 ? cov / (count - 1) : std::numeric_limits<double>::quiet_NaN();
}

// Function to compute Pearson correlation using Armadillo
double PearsonCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  bool NA_rm = false) {
  // Check input sizes
  if (y.size() != y_hat.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  // Handle NA values
  std::vector<double> clean_y, clean_y_hat;
  for (size_t i = 0; i < y.size(); ++i) {
    bool is_na = isNA(y[i]) || isNA(y_hat[i]);
    if (is_na) {
      if (!NA_rm) {
        return std::numeric_limits<double>::quiet_NaN(); // Return NaN if NA_rm is false
      }
    } else {
      clean_y.push_back(y[i]);
      clean_y_hat.push_back(y_hat[i]);
    }
  }

  // If no valid data, return NaN
  if (clean_y.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Convert cleaned vectors to Armadillo vectors
  arma::vec arma_y(clean_y);
  arma::vec arma_y_hat(clean_y_hat);

  // Compute Pearson correlation using Armadillo
  double corr = arma::as_scalar(arma::cor(arma_y, arma_y_hat));

  // Ensure correlation is within valid range [-1, 1]
  if (corr < -1.0) corr = -1.0;
  if (corr > 1.0) corr = 1.0;

  return corr;
}

// // Function to calculate the Pearson correlation coefficient, ignoring NA values
// double PearsonCor(const std::vector<double>& y,
//                   const std::vector<double>& y_hat,
//                   bool NA_rm = false) {
//   // // Check input sizes
//   // if (y.size() != y_hat.size()) {
//   //   throw std::invalid_argument("Input vectors must have the same size.");
//   // }
//
//   // Handle NA values
//   std::vector<double> clean_y, clean_y_hat;
//   for (size_t i = 0; i < y.size(); ++i) {
//     bool is_na = isNA(y[i]) || isNA(y_hat[i]);
//     if (is_na) {
//       if (!NA_rm) {
//         return std::numeric_limits<double>::quiet_NaN(); // Return NaN if NA_rm is false
//       }
//     } else {
//       clean_y.push_back(y[i]);
//       clean_y_hat.push_back(y_hat[i]);
//     }
//   }
//
//   // If no valid data, return NaN
//   if (clean_y.empty()) {
//     return std::numeric_limits<double>::quiet_NaN();
//   }
//
//   double cov_yy_hat = CppCovariance(clean_y, clean_y_hat, true);
//   double var_y = CppVariance(clean_y, true);
//   double var_y_hat = CppVariance(clean_y_hat, true);
//
//   // If any of the values is NaN, return NaN
//   if (isNA(cov_yy_hat) || isNA(var_y) || isNA(var_y_hat)) {
//     return std::numeric_limits<double>::quiet_NaN();
//   }
//
//   // Check if variances are zero
//   if (var_y == 0.0 || var_y_hat == 0.0) {
//     return std::numeric_limits<double>::quiet_NaN(); // Return NaN if variance is zero
//   }
//
//   // Calculate Pearson correlation coefficient
//   double corr = cov_yy_hat / std::sqrt(var_y * var_y_hat);
//
//   // Ensure correlation is within valid range [-1, 1]
//   if (corr < -1.0) corr = -1.0;
//   if (corr > 1.0) corr = 1.0;
//
//   return corr;
// }

// Function to compute Spearman correlation using Armadillo
double SpearmanCor(const std::vector<double>& y,
                   const std::vector<double>& y_hat,
                   bool NA_rm = false) {
  if (y.size() != y_hat.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  std::vector<double> clean_y, clean_y_hat;
  for (size_t i = 0; i < y.size(); ++i) {
    bool is_na = isNA(y[i]) || isNA(y_hat[i]);
    if (is_na) {
      if (!NA_rm) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      clean_y.push_back(y[i]);
      clean_y_hat.push_back(y_hat[i]);
    }
  }

  if (clean_y.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  arma::vec arma_y(clean_y);
  arma::vec arma_y_hat(clean_y_hat);

  // Rank-transform both vectors
  arma::vec rank_y = arma::conv_to<arma::vec>::from(arma::sort_index(arma::sort_index(arma_y))) + 1;
  arma::vec rank_y_hat = arma::conv_to<arma::vec>::from(arma::sort_index(arma::sort_index(arma_y_hat))) + 1;

  // Compute Pearson correlation of ranks
  double corr = arma::as_scalar(arma::cor(rank_y, rank_y_hat));
  return std::max(-1.0, std::min(1.0, corr));
}

// Function to compute Kendall's tau correlation coefficient
double KendallCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  bool NA_rm = false) {
  if (y.size() != y_hat.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  std::vector<double> clean_y, clean_y_hat;
  for (size_t i = 0; i < y.size(); ++i) {
    bool is_na = isNA(y[i]) || isNA(y_hat[i]);
    if (is_na) {
      if (!NA_rm) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    } else {
      clean_y.push_back(y[i]);
      clean_y_hat.push_back(y_hat[i]);
    }
  }

  size_t n = clean_y.size();
  if (n < 2) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Count concordant and discordant pairs
  int concordant = 0, discordant = 0;
  for (size_t i = 0; i < n - 1; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      double dy = clean_y[i] - clean_y[j];
      double dyh = clean_y_hat[i] - clean_y_hat[j];
      double sign = dy * dyh;
      if (sign > 0) ++concordant;
      else if (sign < 0) ++discordant;
      // ties are ignored (as in Kendall's tau-a)
    }
  }

  double denom = 0.5 * n * (n - 1);
  double tau = (concordant - discordant) / denom;
  return std::max(-1.0, std::min(1.0, tau));
}

/*
 * Function to compute Partial Correlation using Armadillo
 *
 * Computes the partial correlation between the dependent variable 'y' and the predicted variable 'y_hat',
 * after controlling for the variables specified in the 'controls' matrix. The partial correlation can be computed
 * either through linear regression or by using the correlation matrix, depending on the 'linear' flag.
 * Optionally, missing values (NA) can be removed if 'NA_rm' is set to true.
 *
 * Parameters:
 *   y          - A vector representing the dependent variable.
 *   y_hat      - A vector representing the predicted variable.
 *   controls   - A matrix where each row corresponds to a control variable to adjust for in the correlation.
 *   NA_rm      - A boolean flag to indicate whether to remove missing values (default is false).
 *   linear     - A boolean flag to specify whether to use linear regression (true) or correlation matrix (false)
 *                for computing the partial correlation (default is false).
 *
 * Returns:
 *   A double representing the partial correlation coefficient between 'y' and 'y_hat' after controlling for
 *   the variables in 'controls'.
 */
double PartialCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  const std::vector<std::vector<double>>& controls,
                  bool NA_rm = false,
                  bool linear = false) {
  // Check input sizes
  if (y.size() != y_hat.size()) {
    throw std::invalid_argument("Input vectors y and y_hat must have the same size.");
  }
  if (!controls.empty()) {
    bool all_controls_valid = std::all_of(controls.begin(), controls.end(),
                                          [&](const std::vector<double>& ctrl) { return ctrl.size() == y.size(); });
    if (!all_controls_valid) {
      throw std::invalid_argument("All control variables must have the same size as y.");
    }
  }

  // Handle NA values
  std::vector<double> clean_y, clean_y_hat;
  std::vector<std::vector<double>> clean_controls(controls.size());

  for (size_t i = 0; i < y.size(); ++i) {
    bool is_na = isNA(y[i]) || isNA(y_hat[i]);
    for (const auto& control : controls) {
      if (isNA(control[i])) {
        is_na = true;
        break;
      }
    }
    if (is_na) {
      if (!NA_rm) {
        return std::numeric_limits<double>::quiet_NaN(); // Return NaN if NA_rm is false
      }
    } else {
      clean_y.push_back(y[i]);
      clean_y_hat.push_back(y_hat[i]);
      for (size_t j = 0; j < controls.size(); ++j) {
        clean_controls[j].push_back(controls[j][i]);
      }
    }
  }

  // If no valid data, return NaN
  if (clean_y.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Check sample adequacy
  if (clean_y.size() <= (clean_controls.size() + 2)) {  // n_samples > n_vars
    return std::numeric_limits<double>::quiet_NaN();
  }

  double partial_corr;
  if (linear){
    // Convert cleaned vectors to Armadillo vectors/matrices
    arma::vec arma_y(clean_y);
    arma::vec arma_y_hat(clean_y_hat);
    arma::mat arma_controls(clean_y.size(), controls.size());
    for (size_t i = 0; i < controls.size(); ++i) {
      arma_controls.col(i) = arma::vec(clean_controls[i]);
    }

    // // Compute residuals of y and y_hat after regressing on controls
    // arma::vec residuals_y = arma_y - arma_controls * arma::solve(arma_controls, arma_y);
    // arma::vec residuals_y_hat = arma_y_hat - arma_controls * arma::solve(arma_controls, arma_y_hat);

    // Use a more robust method for solving the linear system, such as arma::pinv (pseudo-inverse):
    arma::vec residuals_y = arma_y - arma_controls * arma::pinv(arma_controls) * arma_y;
    arma::vec residuals_y_hat = arma_y_hat - arma_controls * arma::pinv(arma_controls) * arma_y_hat;

    // Compute Pearson correlation of the residuals
    partial_corr = arma::as_scalar(arma::cor(residuals_y, residuals_y_hat));

  } else {
    int i = controls.size();
    int j = controls.size() + 1;
    arma::mat data(clean_y.size(), i + 2);
    for (size_t k = 0; k < controls.size(); ++k) {
      data.col(k) = arma::vec(clean_controls[k]);
    }
    data.col(i) = arma::vec(clean_y);
    data.col(j) = arma::vec(clean_y_hat);

    if (data.n_rows < 2 || data.n_cols < 1) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    // Compute the correlation matrix of the data
    arma::mat corrm = arma::cor(data);

    // // Compute the precision matrix (inverse of the correlation matrix)
    // arma::mat precm = arma::inv(corrm);

    // Moore-Penrose pseudo-inverse
    // arma::mat precm = arma::pinv(corrm);
    arma::mat precm;
    try {
      precm = arma::pinv(corrm, 1e-10);
    } catch (...) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    // Get the correlation between y and y_hat after controlling for the others
    partial_corr = -precm(i, j) / std::sqrt(precm(i, i) * precm(j, j));
  }

  // Ensure partial correlation is within valid range [-1, 1]
  if (partial_corr < -1.0) partial_corr = -1.0;
  if (partial_corr > 1.0) partial_corr = 1.0;

  return partial_corr;
}

double PartialCorTrivar(const std::vector<double>& y,
                        const std::vector<double>& y_hat,
                        const std::vector<double>& control,
                        bool NA_rm = false,
                        bool linear = false){
  std::vector<std::vector<double>> conmat;
  conmat.push_back(control);

  double res = PartialCor(y,y_hat,conmat,NA_rm,linear);
  return res;
}

/**
 * Calculates the significance (two-sided p-value) of a (partial) correlation coefficient.
 *
 * This function computes the t-statistic for a given (partial) correlation coefficient `r`
 * and returns the corresponding two-tailed p-value under the null hypothesis that the true
 * correlation is zero.
 *
 * The t-statistic is calculated using:
 *     t = r * sqrt((n - k - 2) / (1 - r^2))
 * where:
 *     - r is the correlation coefficient
 *     - n is the sample size
 *     - k is the number of control variables (0 for simple correlation)
 *
 * The degrees of freedom used is (n - k - 2). The resulting two-sided p-value is computed
 * using the cumulative distribution function of the t-distribution.
 *
 * @param r The (partial) correlation coefficient.
 * @param n The number of observations.
 * @param k The number of control variables (default = 0).
 * @return The two-sided p-value.
 */
double CppCorSignificance(double r, size_t n, size_t k = 0) {
  // Check if degrees of freedom are valid: df = n - k - 2 >= 1
  if (n <= k + 2) {
    return std::numeric_limits<double>::quiet_NaN();  // Invalid degrees of freedom, return non-significant result
  }

  double df = static_cast<double>(n - k - 2);
  double t = r * std::sqrt(df / (1.0 - r * r));

  double pvalue = (1 - R::pt(t, df, true, false)) * 2;

  // Ensure p value is within valid range [-1, 1]
  if (pvalue < 0) pvalue = 0;
  if (pvalue > 1.0) pvalue = 1.0;

  return pvalue;
}

/**
 * Calculates the confidence interval for a (partial) correlation coefficient.
 *
 * This function uses Fisher's z-transformation to compute the confidence interval
 * for a correlation or partial correlation coefficient `r`. The transformation
 * stabilizes the variance of `r` for more accurate interval estimation.
 *
 * The steps include:
 *   1. Transforming r using Fisher's z.
 *   2. Computing the standard error of z.
 *   3. Determining the critical z-value for the specified confidence level.
 *   4. Calculating the confidence interval in the z-domain.
 *   5. Back-transforming to get the interval in the correlation domain.
 *
 * The degrees of freedom are adjusted for partial correlation with `k` control variables.
 *
 * @param r The (partial) correlation coefficient.
 * @param n The number of observations.
 * @param k The number of control variables (default = 0; use 0 for simple correlation).
 * @param level The significance level α for the confidence interval (default = 0.05).
 * @return A vector containing the upper and lower bounds of the confidence interval.
 */
std::vector<double> CppCorConfidence(double r, size_t n, size_t k = 0,
                                     double level = 0.05) {
  if (n <= k + 3) {
    // Not enough degrees of freedom to compute confidence interval
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  }

  // Calculate the Fisher's z-transformation
  double z = 0.5 * std::log((1 + r) / (1 - r));

  // Calculate the standard error of z
  double ztheta = 1 / std::sqrt(static_cast<double>(n - k - 3));

  // Calculate the z-value for the given confidence level
  double qZ = R::qnorm(1 - level / 2, 0.0, 1.0, true, false);

  // Calculate the upper and lower bounds of the confidence interval
  double upper = z + qZ * ztheta;
  double lower = z - qZ * ztheta;

  // Convert the bounds back to correlation coefficients
  double r_upper = (std::exp(2 * upper) - 1) / (std::exp(2 * upper) + 1);
  double r_lower = (std::exp(2 * lower) - 1) / (std::exp(2 * lower) + 1);

  // Clamp results to [-1,1] to prevent numerical issues
  if (r_lower < -1.0) r_lower = -1.0;
  if (r_lower > 1.0) r_lower = 1.0;
  if (r_upper < -1.0) r_upper = -1.0;
  if (r_upper > 1.0) r_upper = 1.0;

  // Return the result as a std::vector<double>
  return {r_upper, r_lower};
}

/**
 * Computes a two-sided p-value for correlation or partial correlation coefficients.
 *
 * This function handles both single and multiple correlation coefficients:
 *
 * - If only one valid correlation coefficient is provided, it uses a t-distribution
 *   to compute the p-value.
 *
 * - If multiple valid correlation coefficients are provided (e.g., from bootstrapping),
 *   it applies Fisher's z-transformation, computes the mean z, and evaluates significance
 *   using the standard normal distribution.
 *
 * Correlation coefficients must be within the valid range (-1, 1); NaN values are ignored.
 *
 * @param rho_vec A vector of correlation or partial correlation coefficients.
 * @param n The number of observations.
 * @param k The number of control variables (0 for simple correlation, >0 for partial correlation).
 * @return The two-sided p-value; returns NaN if no valid rho is found.
 */
double CppMeanCorSignificance(const std::vector<double>& rho_vec,
                              size_t n, size_t k = 0) {
  if (n <= k + 2) return std::numeric_limits<double>::quiet_NaN();  // Invalid degrees of freedom

  std::vector<double> z_values;
  for (double rho : rho_vec) {
    if (std::isnan(rho) || std::fabs(rho) >= 1.0) continue;
    double z = 0.5 * std::log((1 + rho) / (1 - rho));
    z_values.push_back(z);
  }

  size_t m = z_values.size();
  if (m == 0) return std::numeric_limits<double>::quiet_NaN();  // No valid correlations

  if (m == 1) {
    // Use t-statistic for single correlation
    double r = rho_vec[0];
    double df = static_cast<double>(n - k - 2);
    double t = r * std::sqrt(df / (1.0 - r * r));
    double pvalue = 2.0 * R::pt(-std::fabs(t), df, true, false);

    // Clamp p-value to [0,1] to prevent numerical issues
    if (pvalue < 0.0) pvalue = 0.0;
    if (pvalue > 1.0) pvalue = 1.0;

    return pvalue;
  }

  // Use Fisher z transformation for multiple correlations
  double z_mean = std::accumulate(z_values.begin(), z_values.end(), 0.0) / m;
  double se = 1.0 / std::sqrt(static_cast<double>(n - k - 3));
  double z_stat = z_mean / se;

  double pvalue = 2.0 * R::pnorm(-std::fabs(z_stat), 0.0, 1.0, true, false);

  // Clamp p-value to [0,1]
  if (pvalue < 0.0) pvalue = 0.0;
  if (pvalue > 1.0) pvalue = 1.0;

  return pvalue;
}

/**
 * Computes the confidence interval for correlation or partial correlation coefficients.
 *
 * This function handles both single and multiple correlation coefficients:
 *
 * - If only one valid correlation coefficient is provided, it uses Fisher's z-transformation
 *   to compute the confidence interval using the standard normal quantile.
 *
 * - If multiple valid correlation coefficients are provided (e.g., from bootstrapping),
 *   it applies Fisher's z-transformation to each, averages them, and computes the confidence interval
 *   on the transformed scale, then transforms back to correlation scale.
 *
 * Correlation coefficients must be within the valid range (-1, 1); NaN values are ignored.
 *
 * @param rho_vec A vector of correlation or partial correlation coefficients.
 * @param n The number of observations.
 * @param k The number of control variables (0 for simple correlation).
 * @param level The significance level alpha (default = 0.05 for 95% confidence).
 * @return A vector of two elements: lower and upper bounds of the confidence interval;
 *         returns {NaN, NaN} if no valid rho is found.
 */
std::vector<double> CppMeanCorConfidence(const std::vector<double>& rho_vec,
                                         size_t n, size_t k = 0,
                                         double level = 0.05) {
  if (n <= k + 3) {
    // Not enough degrees of freedom
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  }

  std::vector<double> z_values;
  for (double rho : rho_vec) {
    if (std::isnan(rho) || std::fabs(rho) >= 1.0) continue;
    double z = 0.5 * std::log((1 + rho) / (1 - rho));
    z_values.push_back(z);
  }

  size_t m = z_values.size();
  if (m == 0) return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

  double z_mean = std::accumulate(z_values.begin(), z_values.end(), 0.0) / m;
  double se = 1.0 / std::sqrt(static_cast<double>(n - k - 3));
  double z_crit = R::qnorm(1.0 - level / 2.0, 0.0, 1.0, true, false);

  double z_lower = z_mean - z_crit * se;
  double z_upper = z_mean + z_crit * se;

  double r_lower = (std::exp(2 * z_lower) - 1) / (std::exp(2 * z_lower) + 1);
  double r_upper = (std::exp(2 * z_upper) - 1) / (std::exp(2 * z_upper) + 1);

  // Clamp results to [-1,1] to prevent numerical issues
  if (r_lower < -1.0) r_lower = -1.0;
  if (r_lower > 1.0) r_lower = 1.0;
  if (r_upper < -1.0) r_upper = -1.0;
  if (r_upper > 1.0) r_upper = 1.0;

  return {r_lower, r_upper};
}

/**
 * Computes the AUC (theta) and confidence interval using the DeLong method.
 *
 * @reference https://github.com/xrobin/pROC/blob/master/R/delong.R ci_auc_delong function
 *
 * @param cases A vector of scores for the cases (positive class).
 * @param controls A vector of scores for the controls (negative class).
 * @param direction A string indicating the direction of comparison (">" for greater, "<" for less).
 * @param level The confidence level, default is 0.05.
 *
 * @return A vector of four elements:
 *   - theta: The computed AUC value.
 *   - ci_lower: The lower bound of the confidence interval.
 *   - ci_upper: The upper bound of the confidence interval.
 */
std::vector<double> CppDeLongAUCConfidence(const std::vector<double>& cases,
                                           const std::vector<double>& controls,
                                           const std::string& direction,
                                           double level = 0.05) {
  // Get sizes of cases and controls
  size_t m = cases.size();
  size_t n = controls.size();

  // Compute DeLong placements
  DeLongPlacementsRes ret = CppDeLongPlacements(cases, controls, direction);
  double theta = ret.theta;
  std::vector<double> X = ret.X;
  std::vector<double> Y = ret.Y;

  // If there are too few cases or controls, return default values
  if (m <= 1 || n <= 1) {
    return {theta, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()}; // Default values for invalid input
  }

  // Compute variances SX and SY
  double SX = 0.0, SY = 0.0;
  for (size_t i = 0; i < m; ++i) {
    SX += (X[i] - theta) * (X[i] - theta);
  }
  SX /= (m - 1);

  for (size_t i = 0; i < n; ++i) {
    SY += (Y[i] - theta) * (Y[i] - theta);
  }
  SY /= (n - 1);

  // Compute the overall variance S
  double S = SX / m + SY / n;

  // Compute the confidence interval using R::qnorm
  double ci_lower = R::qnorm(level / 2, theta, std::sqrt(S), true, false);
  double ci_upper = R::qnorm(1 - level / 2, theta, std::sqrt(S), true, false);

  // Ensure the confidence interval is within [0, 1]
  ci_lower = std::max(0.0, ci_lower);
  ci_upper = std::min(1.0, ci_upper);

  // Return the results as a three-element vector
  return {theta, ci_upper, ci_lower};
}

/**
 * Computes the AUC (theta), p-value, and confidence interval using the DeLong method.
 *
 * @param cases A vector of scores for the cases (positive class).
 * @param direction A string indicating the direction of comparison (">" for greater, "<" for less).
 * @param level The confidence level, default is 0.05.
 * @param num_samples Number of effective sample size
 *
 * @return A vector of four elements:
 *   - theta: The computed AUC value.
 *   - p_value: The p-value for testing the null hypothesis that AUC = 0.5.
 *   - ci_lower: The lower bound of the confidence interval.
 *   - ci_upper: The upper bound of the confidence interval.
 */
std::vector<double> CppCMCTest(const std::vector<double>& cases,
                               const std::string& direction,
                               double level = 0.05,
                               size_t num_samples = 0) {
  size_t m = cases.size(), n = cases.size();  // Both m and n are set to cases.size()

  if (num_samples == 0){
    num_samples = m;
  }

  std::vector<double> controls;
  for (size_t i = 0; i < cases.size(); ++i) {
    // controls.push_back(static_cast<double>(i) / num_samples);
    controls.push_back(static_cast<double>(i + 1) / num_samples);
  }

  // Compute DeLong placements
  DeLongPlacementsRes ret = CppDeLongPlacements(cases, controls, direction);
  double theta = ret.theta;
  std::vector<double> X = ret.X;
  std::vector<double> Y = ret.Y;

  // If there are too few cases or controls, return default values
  if (m <= 1 || n <= 1) {
    return {theta, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()}; // Default values for invalid input
  }

  // Compute variances SX and SY
  double SX = 0.0, SY = 0.0;
  for (size_t i = 0; i < m; ++i) {
    SX += (X[i] - theta) * (X[i] - theta);
  }
  SX /= (m - 1);

  for (size_t i = 0; i < n; ++i) {
    SY += (Y[i] - theta) * (Y[i] - theta);
  }
  SY /= (n - 1);

  // Compute the overall variance S
  double S = SX / m + SY / n;

  // Compute the Z-score for the p-value
  double z = (theta - 0.5) / std::sqrt(S);

  // Compute the two-tailed p-value (AUC ≠ 0.5)
  double p_value = 2 * R::pnorm(-std::abs(z), 0.0, 1.0, true, false);

  // // Compute the one-sided test (right-tailed) p-value (AUC > 0.5)
  // // double p_value = R::pnorm(z, 0.0, 1.0, true, false);
  // // Set theta to negative when sample size is small to mitigate sample size effect`
  // double p_value = R::pnorm(-std::abs(z), 0.0, 1.0, true, false);

  // Compute the confidence interval using R::qnorm
  double ci_lower = R::qnorm(level / 2, theta, std::sqrt(S), true, false);
  double ci_upper = R::qnorm(1 - level / 2, theta, std::sqrt(S), true, false);

  // Ensure the confidence interval is within [0, 1]
  ci_lower = std::max(0.0, ci_lower);
  ci_upper = std::min(1.0, ci_upper);

  // Return the results as a four-element vector
  return {theta, p_value, ci_upper, ci_lower};
}

// /**
//  * Computes the AUC (via area method), p-value, and confidence interval using the DeLong method.
//  *
//  * @param cases A vector of scores for the cases (positive class).
//  * @param direction A string indicating the direction of comparison (">" or "<").
//  * @param level The confidence level, default is 0.05.
//  * @param num_samples Number of effective sample size for the control group.
//  *
//  * @return A vector of four elements:
//  *   - area_auc: AUC computed by trapezoidal rule (area method).
//  *   - p_value: The p-value testing H0: AUC = 0.5 (via DeLong).
//  *   - ci_upper: Upper bound of 100*(1-level)% confidence interval (via DeLong).
//  *   - ci_lower: Lower bound of 100*(1-level)% confidence interval (via DeLong).
//  */
// std::vector<double> CppCMCTest(const std::vector<double>& cases,
//                                const std::string& direction,
//                                double level = 0.05,
//                                size_t num_samples = 0) {
//   size_t m = cases.size();
//   if (m <= 1) {
//     return {std::numeric_limits<double>::quiet_NaN(), 1.0,
//             std::numeric_limits<double>::quiet_NaN(),
//             std::numeric_limits<double>::quiet_NaN()};
//   }
//
//   if (num_samples == 0) {
//     num_samples = m;
//   }
//
//   // === AUC via area under curve (trapezoidal rule) ===
//   std::vector<double> x(m);
//   for (size_t i = 0; i < m; ++i) {
//     x[i] = static_cast<double>(i) / (m - 1);  // normalized x-axis [0,1]
//   }
//
//   // Compute area AUC using trapezoidal integration
//   double area_auc = 0.0;
//   for (size_t i = 1; i < m; ++i) {
//     area_auc += 0.5 * (cases[i] + cases[i - 1]) * (x[i] - x[i - 1]);
//   }
//
//   // === DeLong p-value and CI ===
//   std::vector<double> controls;
//   // for (size_t i = 0; i < cases.size(); ++i) {
//   //   controls.push_back(static_cast<double>(i) / num_samples);
//   // }
//   for (size_t i = 0; i < m; ++i) {
//     controls.push_back(0.5 * static_cast<double>(i + 1) / num_samples);
//   }
//
//   // DeLong placements
//   DeLongPlacementsRes ret = CppDeLongPlacements(cases, controls, direction);
//   double theta = ret.theta;
//   std::vector<double> X = ret.X;
//   std::vector<double> Y = ret.Y;
//
//   if (X.size() <= 1 || Y.size() <= 1) {
//     return {area_auc, 1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
//   }
//
//   // Compute variances
//   double SX = 0.0, SY = 0.0;
//   for (size_t i = 0; i < m; ++i) {
//     SX += (X[i] - theta) * (X[i] - theta);
//     SY += (Y[i] - theta) * (Y[i] - theta);
//   }
//   SX /= (m - 1);
//   SY /= (m - 1);  // controls.size() == cases.size()
//
//   // Compute the overall variance S
//   double S = SX / m + SY / m;
//
//   // Compute the Z-score for the p-value
//   double z = (theta - 0.5) / std::sqrt(S);
//
//   // Compute the two-tailed p-value (AUC ≠ 0.5)
//   double p_value = 2 * R::pnorm(-std::abs(z), 0.0, 1.0, true, false);
//
//   // // Compute the one-sided test (right-tailed) p-value (AUC > 0.5)
//   // // double p_value = R::pnorm(z, 0.0, 1.0, true, false);
//   // // Set theta to negative when sample size is small to mitigate sample size effect`
//   // double p_value = R::pnorm(-std::abs(z), 0.0, 1.0, true, false);
//
//   // Confidence interval using DeLong variance
//   double ci_lower = R::qnorm(level / 2, theta, std::sqrt(S), true, false);
//   double ci_upper = R::qnorm(1 - level / 2, theta, std::sqrt(S), true, false);
//   ci_lower = std::max(0.0, ci_lower);
//   ci_upper = std::min(1.0, ci_upper);
//
//   return {area_auc, p_value, ci_upper, ci_lower};
// }

// Function to compute distance between two vectors:
double CppDistance(const std::vector<double>& vec1,
                   const std::vector<double>& vec2,
                   bool L1norm = false,
                   bool NA_rm = false){
  // Handle NA values
  std::vector<double> clean_v1, clean_v2;
  for (size_t i = 0; i < vec1.size(); ++i) {
    bool is_na = isNA(vec1[i]) || isNA(vec2[i]);
    if (is_na) {
      if (!NA_rm) {
        return std::numeric_limits<double>::quiet_NaN(); // Return NaN if NA_rm is false
      }
    } else {
      clean_v1.push_back(vec1[i]);
      clean_v2.push_back(vec2[i]);
    }
  }

  // If no valid data, return NaN
  if (clean_v1.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double dist_res = 0.0;
  if (L1norm) {
    for (size_t i = 0; i < clean_v1.size(); ++i) {
      dist_res += std::abs(clean_v1[i] - clean_v2[i]);
    }
  } else {
    for (size_t i = 0; i < clean_v1.size(); ++i) {
      dist_res += (clean_v1[i] - clean_v2[i]) * (clean_v1[i] - clean_v2[i]);
    }
    dist_res = std::sqrt(dist_res);
  }

  return dist_res;
}

// Function to compute the chebyshev distance between two vectors:
double CppChebyshevDistance(const std::vector<double>& vec1,
                            const std::vector<double>& vec2,
                            bool NA_rm = false){
  // Handle NA values
  std::vector<double> clean_v1, clean_v2;
  for (size_t i = 0; i < vec1.size(); ++i) {
    bool is_na = isNA(vec1[i]) || isNA(vec2[i]);
    if (is_na) {
      if (!NA_rm) {
        return std::numeric_limits<double>::quiet_NaN(); // Return NaN if NA_rm is false
      }
    } else {
      clean_v1.push_back(vec1[i]);
      clean_v2.push_back(vec2[i]);
    }
  }

  // If no valid data, return NaN
  if (clean_v1.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double dist_res = 0.0;
  for (size_t i = 0; i < clean_v1.size(); ++i){
    if (dist_res < std::abs(clean_v1[i] - clean_v2[i])){
      dist_res = std::abs(clean_v1[i] - clean_v2[i]);
    }
  }

  return dist_res;
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

// Function to find k-nearest neighbors of a given index in the embedding space
// The `lib` parameter specifies the indices from which the k-nearest neighbors should be selected
std::vector<size_t> CppKNNIndice(
    const std::vector<std::vector<double>>& embedding_space,  // Embedding space containing vectors
    size_t target_idx,                                        // Target index for which to find neighbors
    size_t k,                                                 // Number of nearest neighbors to find
    const std::vector<int>& lib)                              // Indices from which to select neighbors
{
  std::vector<std::pair<double, size_t>> distances;

  // Iterate through the specified library indices to collect valid distances
  for (std::size_t i : lib) {
    if (i == target_idx) continue;  // Skip the target index itself

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
    const std::vector<int>& lib)                       // Indices from which to select neighbors
{
  std::vector<std::pair<double, size_t>> distances;

  // Iterate through the specified library indices to collect valid distances
  for (size_t i : lib) {
    if (i == target_idx) continue;  // Skip the target index itself

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
 * Computes sorted neighbor indices for selected rows in a precomputed distance matrix,
 * with neighbors restricted to a specified subset of indices (lib).
 *
 * Parameters:
 *   dist_mat      - Precomputed n x n distance matrix (may include NaN).
 *   lib           - Indices to compute neighbors for (and restrict neighbors to).
 *   include_self  - Whether to include self (i == j) as a neighbor.
 *
 * Returns:
 *   A vector of n vectors. For rows in lib, each subvector contains the indices of valid neighbors
 *   (from lib, sorted by increasing distance). Other rows are filled with `invalid_index`.
 */
std::vector<std::vector<size_t>> CppDistSortedIndice(
    const std::vector<std::vector<double>>& dist_mat,
    const std::vector<size_t>& lib,
    bool include_self = false)
{
  const size_t n = dist_mat.size();
  const size_t invalid_index = std::numeric_limits<size_t>::max();

  // Initialize all rows as invalid by default
  std::vector<std::vector<size_t>> sorted_indices(n, std::vector<size_t>{invalid_index});

  for (size_t i : lib) {
    if (i >= n || dist_mat[i].size() != n) continue;

    const auto& row = dist_mat[i];

    if (std::isnan(row[i])) continue;

    std::vector<std::pair<double, size_t>> valid_neighbors;

    for (size_t j : lib) {
      if (!include_self && i == j) continue;

      double d = row[j];
      if (!std::isnan(d)) {
        valid_neighbors.emplace_back(d, j);
      }
    }

    std::sort(valid_neighbors.begin(), valid_neighbors.end(),
              [](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) {
                return (a.first < b.first) || (a.first == b.first && a.second < b.second);
              });

    std::vector<size_t> indices;
    for (const auto& pair : valid_neighbors) {
      indices.push_back(pair.second);
    }

    sorted_indices[i] = indices;
  }

  return sorted_indices;
}

/*
 * Function to compute Singular Value Decomposition (SVD) similar to R's svd()
 *
 * This function computes the Singular Value Decomposition (SVD) of the input matrix 'X'.
 * The decomposition breaks the matrix 'X' into three components: singular values, left singular vectors,
 * and right singular vectors. These components are returned in a nested vector structure.
 *
 * Parameters:
 *   - X: A matrix represented as a std::vector<std::vector<double>>.
 *        This matrix is decomposed into its singular values and vectors.
 *
 * Returns:
 *   A std::vector containing three components:
 *     1. d: A vector of singular values (std::vector<double>).
 *     2. u: A matrix of left singular vectors (std::vector<std::vector<double>>).
 *     3. v: A matrix of right singular vectors (std::vector<std::vector<double>>).
 */
std::vector<std::vector<std::vector<double>>> CppSVD(const std::vector<std::vector<double>>& X) {
  // Convert input matrix to Armadillo matrix
  size_t m = X.size();
  size_t n = X[0].size();
  arma::mat A(m, n);
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      A(i, j) = X[i][j];
    }
  }

  // Perform SVD using Armadillo
  arma::mat U; // Left singular vectors
  arma::vec S; // Singular values
  arma::mat V; // Right singular vectors
  arma::svd(U, S, V, A);

  // Convert Armadillo objects back to std::vector
  std::vector<std::vector<double>> u(m, std::vector<double>(m));
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      u[i][j] = U(i, j);
    }
  }

  std::vector<double> d(S.n_elem);
  for (size_t i = 0; i < S.n_elem; ++i) {
    d[i] = S(i);
  }

  std::vector<std::vector<double>> v(V.n_rows, std::vector<double>(V.n_cols));
  for (size_t i = 0; i < V.n_rows; ++i) {
    for (size_t j = 0; j < V.n_cols; ++j) {
      v[i][j] = V(i, j);
    }
  }

  // Return as a std::vector to match R's svd() output
  return {u, {d}, v};
}

// Function to remove linear trend using Armadillo for internal calculations
std::vector<double> LinearTrendRM(const std::vector<double>& vec,
                                  const std::vector<double>& xcoord,
                                  const std::vector<double>& ycoord,
                                  bool NA_rm = false) {
  // Check input sizes
  if (vec.size() != xcoord.size() || vec.size() != ycoord.size()) {
    throw std::invalid_argument("Input vectors must have the same size.");
  }

  // Create a vector to store the result
  std::vector<double> result(vec.size(), std::numeric_limits<double>::quiet_NaN()); // Initialize with NA

  // Find indices where all three vectors are not NA
  std::vector<size_t> valid_indices;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (!isNA(vec[i]) && !isNA(xcoord[i]) && !isNA(ycoord[i])) {
      valid_indices.push_back(i);
    } else if (!NA_rm) {
      throw std::invalid_argument("Input contains NA values and NA_rm is false.");
    }
  }

  // If no valid data, return the result filled with NA
  if (valid_indices.empty()) {
    return result;
  }

  // Extract non-NA values
  std::vector<double> clean_vec, clean_xcoord, clean_ycoord;
  for (size_t i : valid_indices) {
    clean_vec.push_back(vec[i]);
    clean_xcoord.push_back(xcoord[i]);
    clean_ycoord.push_back(ycoord[i]);
  }

  // Convert cleaned vectors to Armadillo vectors
  arma::vec arma_vec(clean_vec);
  arma::vec arma_xcoord(clean_xcoord);
  arma::vec arma_ycoord(clean_ycoord);

  // Create design matrix for linear regression
  arma::mat X(clean_vec.size(), 3);
  X.col(0) = arma::ones<arma::vec>(clean_vec.size()); // Intercept
  X.col(1) = arma_xcoord; // x1
  X.col(2) = arma_ycoord; // x2

  // Perform linear regression using Armadillo
  arma::vec coefficients;
  if (!arma::solve(coefficients, X, arma_vec)) {
    throw std::invalid_argument("Linear regression failed due to singular matrix.");
  }

  // Predict vec_hat using the linear regression model
  arma::vec arma_vec_hat = X * coefficients;

  // Fill the result vector
  for (size_t i = 0; i < valid_indices.size(); ++i) {
    result[valid_indices[i]] = clean_vec[i] - arma_vec_hat(i);
  }

  return result;
}

// // Function to perform Linear Trend Removal
// std::vector<double> LinearTrendRM(const std::vector<double>& vec,
//                                   const std::vector<double>& xcoord,
//                                   const std::vector<double>& ycoord,
//                                   bool NA_rm = false) {
//   if (vec.size() != xcoord.size() || vec.size() != ycoord.size()) {
//     throw std::invalid_argument("Input vectors must have the same size.");
//   }
//
//   // Perform linear regression
//   double x1_mean = CppMean(xcoord, NA_rm);
//   double x2_mean = CppMean(ycoord, NA_rm);
//   double y_mean = CppMean(vec, NA_rm);
//
//   double x1_var = CppVariance(xcoord, NA_rm);
//   double x2_var = CppVariance(ycoord, NA_rm);
//   double x1_x2_cov = CppCovariance(xcoord, ycoord, NA_rm);
//
//   double x1_y_cov = CppCovariance(xcoord, vec, NA_rm);
//   double x2_y_cov = CppCovariance(ycoord, vec, NA_rm);
//
//   double denom = x1_var * x2_var - x1_x2_cov * x1_x2_cov;
//   if (denom == 0.0) {
//     throw std::invalid_argument("Linear regression cannot be performed due to collinearity.");
//   }
//
//   double b1 = (x2_var * x1_y_cov - x1_x2_cov * x2_y_cov) / denom;
//   double b2 = (x1_var * x2_y_cov - x1_x2_cov * x1_y_cov) / denom;
//   double b0 = y_mean - b1 * x1_mean - b2 * x2_mean;
//
//   // Predict vec_hat using the linear regression model
//   std::vector<double> vec_hat(vec.size());
//   for (size_t i = 0; i < vec.size(); ++i) {
//     vec_hat[i] = b0 + b1 * xcoord[i] + b2 * ycoord[i];
//   }
//
//   // Calculate vec - vec_hat
//   std::vector<double> result(vec.size());
//   for (size_t i = 0; i < vec.size(); ++i) {
//     result[i] = vec[i] - vec_hat[i];
//   }
//
//   return result;
// }
