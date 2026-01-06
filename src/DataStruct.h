#ifndef DataStruct_H
#define DataStruct_H

#include <vector>
#include <cstdint>
#include <limits>
#include <string>
#include <utility> // for std::move

struct PartialCorRes {
  int first;
  double second;
  double third;

  /// Default constructor: initializes all members to zero.
  PartialCorRes() : first(0), second(0.0), third(0.0) {}

  /// Parameterized constructor: initializes all members.
  PartialCorRes(int f, double s, double t) : first(f), second(s), third(t) {}

  /// Defaulted copy/move constructors and assignment operators.
  PartialCorRes(const PartialCorRes&) = default;
  PartialCorRes& operator=(const PartialCorRes&) = default;
  PartialCorRes(PartialCorRes&&) noexcept = default;
  PartialCorRes& operator=(PartialCorRes&&) noexcept = default;
};

struct CMCRes {
  std::vector<double> cross_mapping;                 ///< Cross mapping values.
  std::vector<std::vector<double>> causal_strength;  ///< Causal strength matrix.

  /// Default constructor.
  CMCRes() = default;

  /// Copy constructor (for lvalues).
  CMCRes(const std::vector<double>& cross_mapping,
         const std::vector<std::vector<double>>& causal_strength)
    : cross_mapping(cross_mapping), causal_strength(causal_strength) {}

  /// Move constructor (for rvalues).
  CMCRes(std::vector<double>&& cross_mapping,
         std::vector<std::vector<double>>&& causal_strength) noexcept
    : cross_mapping(std::move(cross_mapping)), causal_strength(std::move(causal_strength)) {}

  /// Defaulted copy/move assignment operators.
  CMCRes(const CMCRes&) = default;
  CMCRes& operator=(const CMCRes&) = default;
  CMCRes(CMCRes&&) noexcept = default;
  CMCRes& operator=(CMCRes&&) noexcept = default;
};

/**
 * @brief Represents the result of an intersection calculation.
 *
 * This structure stores:
 * - `libsize`: The library size.
 * - `Intersection`: The intersection values as a vector of doubles.
 *
 */
struct IntersectionRes {
  size_t libsize;                   ///< Library size.
  std::vector<double> Intersection; ///< Intersection values.

  /// Default constructor: initializes libsize to 0 and empty Intersection vector.
  IntersectionRes() : libsize(0), Intersection() {}

  /// Copy constructor (for lvalues): deep copies Intersection vector.
  IntersectionRes(size_t t, const std::vector<double>& x)
    : libsize(t), Intersection(x) {}

  /// Move constructor (for rvalues): moves Intersection vector.
  IntersectionRes(size_t t, std::vector<double>&& x) noexcept
    : libsize(t), Intersection(std::move(x)) {}

  /// Defaulted copy/move assignment operators.
  IntersectionRes(const IntersectionRes&) = default;
  IntersectionRes& operator=(const IntersectionRes&) = default;
  IntersectionRes(IntersectionRes&&) noexcept = default;
  IntersectionRes& operator=(IntersectionRes&&) noexcept = default;
};

/**
 * @brief Holds the results of DeLong placement calculations.
 *
 * This structure stores:
 * - `theta`: The AUC or placement statistic.
 * - `X`: Placement values for the positive class.
 * - `Y`: Placement values for the negative class.
 *
 * The move constructor and move assignment are marked `noexcept`
 * to enable efficient transfers when stored in standard containers
 * like `std::vector`.
 */
struct DeLongPlacementsRes {
  double theta;                  ///< The AUC or placement statistic.
  std::vector<double> X;         ///< Placement values for the positive class.
  std::vector<double> Y;         ///< Placement values for the negative class.

  /**
   * @brief Default constructor.
   * Initializes theta = 0.0 and empty vectors.
   */
  DeLongPlacementsRes() : theta(0.0), X(), Y() {}

  /**
   * @brief Copy constructor (for lvalues).
   * Performs deep copies of input vectors.
   *
   * @param t The statistic value (e.g., theta).
   * @param x The placement values for the positive class.
   * @param y The placement values for the negative class.
   */
  DeLongPlacementsRes(double t, const std::vector<double>& x, const std::vector<double>& y)
    : theta(t), X(x), Y(y) {}

  /**
   * @brief Move constructor (for rvalues).
   * Transfers ownership of vectors without copying.
   */
  DeLongPlacementsRes(double t, std::vector<double>&& x, std::vector<double>&& y) noexcept
    : theta(t), X(std::move(x)), Y(std::move(y)) {}

  /// Defaulted copy/move assignment operators.
  DeLongPlacementsRes(const DeLongPlacementsRes&) = default;
  DeLongPlacementsRes& operator=(const DeLongPlacementsRes&) = default;
  DeLongPlacementsRes(DeLongPlacementsRes&&) noexcept = default;
  DeLongPlacementsRes& operator=(DeLongPlacementsRes&&) noexcept = default;
};

#endif // DataStruct_H
