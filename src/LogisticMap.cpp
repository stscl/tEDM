#include <cmath>
#include <limits>
#include <vector>
#include <numeric>

/**
 * @brief Simulate a univariate logistic map.
 *
 * @param val                 Initial value.
 * @param step                Number of simulation time steps to run.
 * @param alpha               Growth/interaction parameter in the logistic update rule.
 * @param escape_threshold    Threshold to treat divergent values as invalid (default: 1e10).
 *
 * @return A vector of simulation results.
 */
std::vector<double> LogisticMapUni(
    double val,
    int step,
    double alpha,
    double escape_threshold = 1e10
){
  // Initialize result vector with NaNs
  std::vector<double> res(step + 1,
                          std::numeric_limits<double>::quiet_NaN());

  // Set initial values at time step 0
  res[0] = val;

  // Time-stepped simulation
  for (int s = 1; s <= step; ++s){
    // Skip if the current value is invalid (NaN)
    if (std::isnan(res[s - 1])) continue;

    // Apply the logistic map update if no neighbors exist
    double v_next = res[s - 1] * (alpha - alpha * res[s - 1]);

    // Update result only if the value is within the escape threshold
    if (!std::isinf(v_next) && std::abs(v_next) <= escape_threshold){
      res[s] = v_next;
    }
  }

  return res;
}

/**
 * @brief Simulate a bivariate logistic map.
 *
 * @param v1                 Initial values of the first variable (e.g., species A density).
 * @param v2                 Initial values of the second variable (e.g., species B density).
 * @param step               Number of simulation time steps to run.
 * @param alpha1             Growth/interaction parameter for the first variable.
 * @param alpha2             Growth/interaction parameter for the second variable.
 * @param beta12             Cross-inhibition coefficient from variable 1 to variable 2.
 * @param beta21             Cross-inhibition coefficient from variable 2 to variable 1.
 * @param escape_threshold   Threshold to treat divergent values as invalid (default: 1e10).
 *
 * @return A 2D vector of simulation results:
 *         - First dimension: variable index (0 for v1, 1 for v2),
 *         - Second dimension: time steps (0 to step).
 */
std::vector<std::vector<double>> LogisticMapBi(
    double v1,
    double v2,
    int step,
    double alpha1,
    double alpha2,
    double beta12,
    double beta21,
    double escape_threshold = 1e10
){
  // Initialize result vector with NaNs
  std::vector<std::vector<double>> res (2,
                                        std::vector<double>(step + 1,
                                                            std::numeric_limits<double>::quiet_NaN()));

  // Set initial values at time step 0
  res[0][0] = v1;
  res[1][0] = v2;

  // Time-stepped simulation
  for (int s = 1; s <= step; ++s){
    // Skip if the current value is invalid (NaN)
    if (std::isnan(res[0][s - 1]) && std::isnan(res[1][s - 1])) continue;

    // Apply the logistic map update if no neighbors exist
    double v_next_1 = res[0][s - 1] * (alpha1 - alpha1 * res[0][s - 1] - beta21 * res[1][s - 1]);
    double v_next_2 = res[1][s - 1] * (alpha2 - alpha2 * res[1][s - 1] - beta12 * res[0][s - 1]);

    // Update result only if the value is within the escape threshold
    if (!std::isinf(v_next_1) && std::abs(v_next_1) <= escape_threshold){
      res[0][s] = v_next_1;
    }
    if (!std::isinf(v_next_2) && std::abs(v_next_2) <= escape_threshold){
      res[1][s] = v_next_2;
    }
  }

  return res;
}

/**
 * @brief Simulate a trivariate logistic map.
 *
 * @param v1                 Initial values of the first variable (e.g., species A density).
 * @param v2                 Initial values of the second variable (e.g., species B density).
 * @param v3                 Initial values of the third variable (e.g., species C density).
 * @param step               Number of simulation time steps to perform.
 * @param alpha1             Growth/interaction parameter for variable 1.
 * @param alpha2             Growth/interaction parameter for variable 2.
 * @param alpha3             Growth/interaction parameter for variable 3.
 * @param beta12             Cross-inhibition from variable 1 to variable 2.
 * @param beta13             Cross-inhibition from variable 1 to variable 3.
 * @param beta21             Cross-inhibition from variable 2 to variable 1.
 * @param beta23             Cross-inhibition from variable 2 to variable 3.
 * @param beta31             Cross-inhibition from variable 3 to variable 1.
 * @param beta32             Cross-inhibition from variable 3 to variable 2.
 * @param escape_threshold   Threshold beyond which values are treated as divergent (default: 1e10).
 *
 * @return A 2D vector of simulation results:
 *         - First dimension: variable index (0 for v1, 1 for v2, 2 for v3),
 *         - Second dimension: time steps (0 to step).
 */
std::vector<std::vector<double>> LogisticMapTri(
    double v1,
    double v2,
    double v3,
    int step,
    double alpha1,
    double alpha2,
    double alpha3,
    double beta12,
    double beta13,
    double beta21,
    double beta23,
    double beta31,
    double beta32,
    double escape_threshold = 1e10
){
  // Initialize result vector with NaNs
  std::vector<std::vector<double>> res (3,
                                        std::vector<double>(step + 1,
                                                            std::numeric_limits<double>::quiet_NaN()));

  // Set initial values at time step 0
  res[0][0] = v1;
  res[1][0] = v2;
  res[2][0] = v3;

  // Time-stepped simulation
  for (int s = 1; s <= step; ++s){
    // Skip if the current value is invalid (NaN)
    if (std::isnan(res[0][s - 1]) &&
        std::isnan(res[1][s - 1]) &&
        std::isnan(res[2][s - 1])) continue;

    // Apply the logistic map update if no neighbors exist
    double v_next_1 = res[0][s - 1] * (alpha1 - alpha1 * res[0][s - 1] - beta21 * res[1][s - 1] - beta31 * res[2][s - 1]);
    double v_next_2 = res[1][s - 1] * (alpha2 - alpha2 * res[1][s - 1] - beta12 * res[0][s - 1] - beta32 * res[2][s - 1]);
    double v_next_3 = res[2][s - 1] * (alpha3 - alpha3 * res[2][s - 1] - beta13 * res[0][s - 1] - beta23 * res[1][s - 1]);

    // Update result only if the value is within the escape threshold
    if (!std::isinf(v_next_1) && std::abs(v_next_1) <= escape_threshold){
      res[0][s] = v_next_1;
    }
    if (!std::isinf(v_next_2) && std::abs(v_next_2) <= escape_threshold){
      res[1][s] = v_next_2;
    }
    if (!std::isinf(v_next_3) && std::abs(v_next_3) <= escape_threshold){
      res[2][s] = v_next_3;
    }
  }

  return res;
}
