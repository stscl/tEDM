#ifndef CPP_COMBN_H
#define CPP_COMBN_H

#include <vector>
#include <string>
#include <functional>

/**
 * @brief Generate all combinations of m elements from a given vector vec.
 *
 * @tparam T The type of elements in the vector.
 * @param vec The input vector to generate combinations from.
 * @param m The number of elements in each combination.
 * @return std::vector<std::vector<T>> A vector containing all combinations.
 */
template <typename T>
std::vector<std::vector<T>> CppCombn(const std::vector<T>& vec, int m) {
  std::vector<std::vector<T>> result;
  std::vector<T> current;

  int vec_size = static_cast<int>(vec.size());

  std::function<void(int)> combnHelper = [&](int start) {
    if (static_cast<int>(current.size()) == m) {
      result.push_back(current);
      return;
    }
    int remaining = m - static_cast<int>(current.size());
    for (int i = start; i <= vec_size - remaining; ++i) {
      current.push_back(vec[i]);
      combnHelper(i + 1);
      current.pop_back();
    }
  };

  combnHelper(0);
  return result;
}

/**
 * @brief Generate all non-empty subsets of a given vector.
 *
 * This function generates all subsets of the input vector with sizes from 1 to vec.size().
 * Internally it calls CppCombn repeatedly for all sizes.
 *
 * @tparam T The type of elements in the vector.
 * @param set The input vector to generate subsets from.
 * @return std::vector<std::vector<T>> A vector containing all non-empty subsets.
 */
template <typename T>
std::vector<std::vector<T>> CppGenSubsets(const std::vector<T>& set) {
  std::vector<std::vector<T>> allSubsets;
  for (int m = 1; m <= static_cast<int>(set.size()); ++m) {
    std::vector<std::vector<T>> combs = CppCombn(set, m);
    allSubsets.insert(allSubsets.end(), combs.begin(), combs.end());
  }
  return allSubsets;
}

#endif // CPP_COMBN_H
