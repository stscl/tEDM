#ifndef CppDistances_H
#define CppDistances_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <utility>
#include <stdexcept>
#include "NumericUtils.h"
#include <RcppThread.h>

double CppDistance(const std::vector<double>& vec1,
                   const std::vector<double>& vec2,
                   bool L1norm = false,
                   bool NA_rm = false);

double CppChebyshevDistance(const std::vector<double>& vec1,
                            const std::vector<double>& vec2,
                            bool NA_rm = false);

std::vector<double> CppKNearestDistance(const std::vector<double>& vec, size_t k,
                                        bool L1norm = false, bool NA_rm = false);

std::vector<double> CppMatKNearestDistance(const std::vector<std::vector<double>>& mat,
                                           size_t k, bool NA_rm = false);

std::vector<std::vector<double>> CppMatDistance(
    const std::vector<std::vector<double>>& mat,
    bool L1norm = false,
    bool NA_rm = false);

std::vector<std::vector<double>> CppMatChebyshevDistance(
    const std::vector<std::vector<double>>& mat,
    bool NA_rm = false);

std::vector<int> CppNeighborsNum(
    const std::vector<double>& vec,
    const std::vector<double>& radius,
    bool equal = false,
    bool L1norm = false,
    bool NA_rm = false);

std::vector<int> CppMatNeighborsNum(
    const std::vector<std::vector<double>>& mat,
    const std::vector<double>& radius,
    bool equal = false,
    bool NA_rm = false);

std::vector<size_t> CppKNNIndice(
    const std::vector<std::vector<double>>& embedding_space,
    size_t target_idx,
    size_t k,
    const std::vector<int>& lib,
    bool include_self = false);

std::vector<size_t> CppDistKNNIndice(
    const std::vector<std::vector<double>>& dist_mat,
    size_t target_idx,
    size_t k,
    const std::vector<int>& lib,
    bool include_self = false);

std::vector<std::vector<size_t>> CppDistSortedIndice(
    const std::vector<std::vector<double>>& dist_mat,
    const std::vector<size_t>& lib,
    size_t k,
    bool include_self = false);

std::vector<std::vector<size_t>> CppMatKNNeighbors(
    const std::vector<std::vector<double>>& embedding_space,
    const std::vector<size_t>& lib,
    size_t k,
    size_t threads,
    bool L1norm = false,
    bool include_self = false);

#endif // CppDistances_H
