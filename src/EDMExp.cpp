#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "CppStats.h"
#include "Embed.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "CCM.h"
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

// Wrapper function to perform convergent cross mapping for time serise data
// predict y based on x ====> x xmap y ====> y causes x
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppCCM(const Rcpp::NumericVector& x,
                            const Rcpp::NumericVector& y,
                            const Rcpp::IntegerVector& libsizes,
                            const Rcpp::IntegerVector& lib,
                            const Rcpp::IntegerVector& pred,
                            int E,
                            int tau,
                            int b,
                            bool simplex,
                            double theta,
                            int threads,
                            int parallel_level,
                            bool progressbar) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);
  std::vector<int> lib_std;
  std::vector<int> pred_std;

  // Check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  int n = y_std.size();
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(y_std[lib[i] - 1]) && (lib[i] < n - E * tau)) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(y_std[pred[i] - 1])&& (pred[i] < n - E * tau)) {
      pred_std.push_back(pred[i] - 1);
    }
  }

  // Perform GCCM Lattice
  std::vector<std::vector<double>> result = CCM(
    x_std,
    y_std,
    libsizes_std,
    lib_std,
    pred_std,
    E,
    tau,
    b,
    simplex,
    theta,
    threads,
    parallel_level,
    progressbar);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  Rcpp::NumericMatrix resultMatrix(result.size(), 5);
  for (size_t i = 0; i < result.size(); ++i) {
    resultMatrix(i, 0) = result[i][0];
    resultMatrix(i, 1) = result[i][1];
    resultMatrix(i, 2) = result[i][2];
    resultMatrix(i, 3) = result[i][3];
    resultMatrix(i, 4) = result[i][4];
  }

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create("libsizes",
                 "x_xmap_y_mean","x_xmap_y_sig",
                 "x_xmap_y_upper","x_xmap_y_lower");
  return resultMatrix;
}
