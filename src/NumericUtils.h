#ifndef NUMERIC_UTILS_H
#define NUMERIC_UTILS_H

#include <cmath>
#include <algorithm>
#include <limits>
#include <initializer_list>

/**
 * @file NumericUtils.h
 * @brief Utility functions for safe and consistent floating-point operations.
 *
 * Provides helper functions for:
 *   - Floating-point comparison with combined relative and absolute tolerance.
 *   - Portable numeric constants (epsilon and tolerance).
 *
 * Intended for scientific computation where double precision stability matters.
 *
 * Author: Wenbo Lv
 * Created: 2025-11-11
 * License: GPL-3
 */

// ==============================
// Common numeric constants
// ==============================
constexpr double DOUBLE_EPS = std::numeric_limits<double>::epsilon();   // ≈ 2.22e-16
constexpr double DOUBLE_TOL_ABS = 1.5e-16;  // Absolute tolerance
constexpr double DOUBLE_TOL_REL = 1.5e-8;   // Relative tolerance

// ==============================
// Floating-point comparison
// ==============================
/**
 * @brief Compare two double values with combined relative and absolute tolerance.
 *
 * Implements a numerically stable test for "near equality":
 * |x - y| <= max(rel_tol * max(|x|, |y|, 1.0), abs_tol)
 *
 * @param x First value
 * @param y Second value
 * @param rel_tol Relative tolerance (default DOUBLE_TOL_REL)
 * @param abs_tol Absolute tolerance (default DOUBLE_TOL_ABS)
 * @return true if x and y are considered equal within tolerance
 */
inline bool doubleNearlyEqual(double x, double y,
                              double rel_tol = DOUBLE_TOL_REL,
                              double abs_tol = DOUBLE_TOL_ABS) noexcept {
  double diff = std::fabs(x - y);
  double scale = std::max({1.0, std::fabs(x), std::fabs(y)});
  return diff <= std::max(rel_tol * scale, abs_tol);
}

#endif // NUMERIC_UTILS_H
