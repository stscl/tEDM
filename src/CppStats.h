#ifndef CppStats_H
#define CppStats_H

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

bool isNA(double value);

bool checkIntNA(int value);

bool checkOneDimVectorHasNaN(const std::vector<double>& vec);

int checkOneDimVectorNotNanNum(const std::vector<double>& vec);

unsigned long long CppFactorial(unsigned int n);

unsigned long long CppCombine(unsigned int n, unsigned int k);

double CppDigamma(double x);

double CppLog(double x, double base = 10);

double CppMedian(const std::vector<double>& vec,
                 bool NA_rm = false);

double CppMean(const std::vector<double>& vec,
               bool NA_rm = false);

double CppMin(const std::vector<double>& vec,
              bool NA_rm = false);

double CppMax(const std::vector<double>& vec,
              bool NA_rm = false);

double CppSum(const std::vector<double>& vec,
              bool NA_rm = false);

double CppMAE(const std::vector<double>& x1,
              const std::vector<double>& x2,
              bool NA_rm = false);

double CppRMSE(const std::vector<double>& x1,
               const std::vector<double>& x2,
               bool NA_rm = false);

std::vector<double> CppCumSum(const std::vector<double>& vec);

std::vector<double> CppAbsDiff(const std::vector<double>& vec1,
                               const std::vector<double>& vec2);

std::vector<double> CppSumNormalize(const std::vector<double>& vec,
                                    bool NA_rm = false);

std::vector<double> CppArithmeticSeq(double from, double to, size_t length_out);

double CppVariance(const std::vector<double>& vec, bool NA_rm = false);

double CppCovariance(const std::vector<double>& vec1,
                     const std::vector<double>& vec2,
                     bool NA_rm = false);

double PearsonCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  bool NA_rm = false);

double SpearmanCor(const std::vector<double>& y,
                   const std::vector<double>& y_hat,
                   bool NA_rm = false);

double KendallCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  bool NA_rm = false);

double PartialCor(const std::vector<double>& y,
                  const std::vector<double>& y_hat,
                  const std::vector<std::vector<double>>& controls,
                  bool NA_rm = false,
                  bool linear = false);

double PartialCorTrivar(const std::vector<double>& y,
                        const std::vector<double>& y_hat,
                        const std::vector<double>& control,
                        bool NA_rm = false,
                        bool linear = false);

double CppCorSignificance(double r, size_t n, size_t k = 0);

std::vector<double> CppCorConfidence(double r, size_t n, size_t k = 0,
                                     double level = 0.05);

double CppMeanCorSignificance(const std::vector<double>& rho_vec,
                              size_t n, size_t k = 0);

std::vector<double> CppMeanCorConfidence(const std::vector<double>& rho_vec,
                                         size_t n, size_t k = 0,
                                         double level = 0.05);

std::vector<double> CppDeLongAUCConfidence(const std::vector<double>& cases,
                                           const std::vector<double>& controls,
                                           const std::string& direction,
                                           double level = 0.05);

std::vector<double> CppCMCTest(const std::vector<double>& cases,
                               const std::string& direction,
                               double level = 0.05,
                               size_t num_samples = 0);

std::vector<double> CppDeLongTest(const std::vector<double>& cases,
                                  const std::vector<double>& controls,
                                  const std::string& direction,
                                  double level = 0.05);

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
    const std::vector<int>& lib);

std::vector<size_t> CppDistKNNIndice(
    const std::vector<std::vector<double>>& dist_mat,
    size_t target_idx,
    size_t k,
    const std::vector<int>& lib);

std::vector<std::vector<size_t>> CppDistSortedIndice(
    const std::vector<std::vector<double>>& dist_mat,
    const std::vector<size_t>& lib,
    bool include_self = false);

std::vector<std::vector<std::vector<double>>> CppSVD(const std::vector<std::vector<double>>& X);

std::vector<double> LinearTrendRM(const std::vector<double>& vec,
                                  const std::vector<double>& xcoord,
                                  const std::vector<double>& ycoord,
                                  bool NA_rm = false);

#endif // CppStats_H
