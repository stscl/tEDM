#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "CppStats.h"
#include "Embed.h"
#include "SimplexProjection.h"
#include "SMap.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Wrapper function to generate time-delay embeddings for a univariate time series
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppEmbed(const Rcpp::NumericVector& vec,
                              int E,
                              int tau) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Generate embeddings
  std::vector<std::vector<double>> embeddings = Embed(vec_std, E, tau);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  int rows = embeddings.size();
  int cols = embeddings[0].size();
  Rcpp::NumericMatrix result(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      result(i, j) = embeddings[i][j];
    }
  }

  return result;
}

// Wrapper function to generate simplex projection forecasting
// [[Rcpp::export]]
Rcpp::NumericVector RcppSimplexForecast(
    const Rcpp::NumericVector& source,
    const Rcpp::NumericVector& target,
    int E,
    int tau,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> source_std = Rcpp::as<std::vector<double>>(source);
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);

  // Generate time-delay embeddings for a univariate time series
  std::vector<std::vector<double>> embeddings = Embed(source_std, E, tau);

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target_std.size();
  // Convert lib and pred (1-based in R) to 0-based indices and check validity
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 0 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 0 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
  }

  // Call the SimplexProjectionPrediction function
  std::vector<double> pred_res = SimplexProjectionPrediction(
    embeddings,
    target_std,
    lib_indices,
    pred_indices,
    num_neighbors
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(pred_res);
}

// Wrapper function to generate s-mapping forecasting
// [[Rcpp::export]]
Rcpp::NumericVector RcppSMapForecast(
    const Rcpp::NumericVector& source,
    const Rcpp::NumericVector& target,
    int E,
    int tau,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors,
    const double& theta){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> source_std = Rcpp::as<std::vector<double>>(source);
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);

  // Generate time-delay embeddings for a univariate time series
  std::vector<std::vector<double>> embeddings = Embed(source_std, E, tau);

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target_std.size();
  // Convert lib and pred (1-based in R) to 0-based indices and check validity
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 0 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 0 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
  }

  // Call the SMapPrediction function
  std::vector<double> pred_res = SMapPrediction(
    embeddings,
    target_std,
    lib_indices,
    pred_indices,
    num_neighbors,
    theta
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(pred_res);
}
