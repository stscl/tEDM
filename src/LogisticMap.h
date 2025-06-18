#ifndef LogisticMap_H
#define LogisticMap_H

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
);

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
);

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
);

#endif // LogisticMap_H
