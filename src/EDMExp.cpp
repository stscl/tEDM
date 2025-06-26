#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "CppStats.h"
#include "LogisticMap.h"
#include "Embed.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "Forecast4TS.h"
#include "IntersectionCardinality.h"
#include "CCM.h"
#include "PCM.h"
#include "CMC.h"
#include "MultispatialCCM.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Wrapper function to perform trivariate logistic map
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppLogisticMap(
    double x = 3.6,
    double y = 3.72,
    double z = 3.68,
    int step = 20,
    double alpha_x = 0.625,
    double alpha_y = 0.77,
    double alpha_z = 0.55,
    double beta_xy = 0.05,
    double beta_xz = 0.05,
    double beta_yx = 0.4,
    double beta_yz = 0.4,
    double beta_zx = 0.65,
    double beta_zy = 0.65,
    double escape_threshold = 1e10
) {
  // Call the core function
  std::vector<std::vector<double>> result = LogisticMapTri(
    x, y, z, step, alpha_x, alpha_y, alpha_z,
    beta_xy, beta_xz, beta_yx, beta_yz, beta_zx, beta_zy,
    escape_threshold
  );

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_cols = step + 1;
  Rcpp::NumericVector out_x(n_cols);
  Rcpp::NumericVector out_y(n_cols);
  Rcpp::NumericVector out_z(n_cols);

  // Copy data into NumericVector
  for (int j = 0; j < n_cols; ++j) {
    out_x(j) = result[0][j];
    out_y(j) = result[1][j];
    out_z(j) = result[2][j];
  }

  // Wrap results into an Rcpp::List
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("x") = out_x,
    Rcpp::Named("y") = out_y,
    Rcpp::Named("z") = out_z
  );

  return out;
}

// Wrapper function to generate time-delay embeddings for a univariate time series
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppEmbed(const Rcpp::NumericVector& vec,
                              int E,
                              int tau = 0) {
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

// Wrapper function to perform simplex projection forecasting
// [[Rcpp::export(rng = false)]]
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
  int max_lag = (tau == 0) ? (E - 1) : (E * tau);
  // Convert lib and pred (1-based in R) to 0-based indices and check validity
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 0 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(source_std[lib[i] - 1]) &&
        !std::isnan(target_std[lib[i] - 1]) &&
        (lib[i] > max_lag + 1)) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 0 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(source_std[pred[i] - 1]) &&
        !std::isnan(target_std[pred[i] - 1]) &&
        (pred[i] > max_lag + 1)) {
      pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
    }
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

// Wrapper function to perform s-mapping forecasting
// [[Rcpp::export(rng = false)]]
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
  int max_lag = (tau == 0) ? (E - 1) : (E * tau);
  // Convert lib and pred (1-based in R) to 0-based indices and check validity
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 0 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(source_std[lib[i] - 1]) &&
        !std::isnan(target_std[lib[i] - 1]) &&
        (lib[i] > max_lag + 1)) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 0 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(source_std[pred[i] - 1]) &&
        !std::isnan(target_std[pred[i] - 1]) &&
        (pred[i] > max_lag + 1)) {
      pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
    }
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

// Wrapper function to compute the intersection cardinality curve
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppIntersectionCardinality(
    const Rcpp::NumericVector& source,
    const Rcpp::NumericVector& target,
    int E,
    int tau,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const int& num_neighbors = 4,
    const int& n_excluded = 0,
    const int& threads = 8,
    const int& parallel_level = 0){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> source_std = Rcpp::as<std::vector<double>>(source);
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);

  // Generate embeddings
  std::vector<std::vector<double>> embedding_x = Embed(source_std, E, tau);
  std::vector<std::vector<double>> embedding_y = Embed(target_std, E, tau);

  // Initialize lib_indices and pred_indices
  std::vector<size_t> lib_indices;
  std::vector<size_t> pred_indices;

  int target_len = target_std.size();
  int max_lag = (tau == 0) ? (E - 1) : (E * tau);
  // Convert lib and pred (1-based in R) to 0-based indices and check validity
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 0 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(source_std[lib[i] - 1]) &&
        !std::isnan(target_std[lib[i] - 1]) &&
        (lib[i] > max_lag + 1)){
        lib_indices.push_back(static_cast<size_t>(lib[i] - 1)); // Convert to 0-based index
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 0 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(source_std[pred[i] - 1]) &&
        !std::isnan(target_std[pred[i] - 1]) &&
        (pred[i] > max_lag + 1)){
        pred_indices.push_back(static_cast<size_t>(pred[i] - 1)); // Convert to 0-based index
    }
  }

  if (lib_indices.size() < static_cast<size_t>(num_neighbors)){
    Rcpp::stop("Library size must not exceed the number of nearest neighbors used for mapping.");
  }

  // Call the IntersectionCardinality function
  std::vector<double> res = IntersectionCardinality(
    embedding_x,
    embedding_y,
    lib_indices,
    pred_indices,
    static_cast<size_t>(num_neighbors),
    static_cast<size_t>(n_excluded),
    threads,
    parallel_level
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(res);
}

//  Wrapper function to help determining embedding dimension `E` and numbers of neighbors `k` parameters
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSimplex4TS(const Rcpp::NumericVector& source,
                                   const Rcpp::NumericVector& target,
                                   const Rcpp::IntegerVector& lib,
                                   const Rcpp::IntegerVector& pred,
                                   const Rcpp::IntegerVector& E,
                                   const Rcpp::IntegerVector& b,
                                   int tau,
                                   int threads) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> source_std = Rcpp::as<std::vector<double>>(source);
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target_std.size();
  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    if (lib[i] < 1 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(source_std[lib[i] - 1]) &&
        !std::isnan(target_std[lib[i] - 1])) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(source_std[pred[i] - 1]) &&
        !std::isnan(target_std[pred[i] - 1])) {
      pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
    }
  }

  std::vector<std::vector<double>> res_std = Simplex4TS(
    source_std,
    target_std,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau,
    threads);

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "rho", "mae", "rmse");
  return result;
}

//  Wrapper function to help determining theta parameters
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSMap4TS(const Rcpp::NumericVector& source,
                                const Rcpp::NumericVector& target,
                                const Rcpp::IntegerVector& lib,
                                const Rcpp::IntegerVector& pred,
                                const Rcpp::NumericVector& theta,
                                int E,
                                int tau,
                                int b,
                                int threads) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> source_std = Rcpp::as<std::vector<double>>(source);
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);
  std::vector<double> theta_std = Rcpp::as<std::vector<double>>(theta);

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target_std.size();
  int max_lag = (tau == 0) ? (E - 1) : (E * tau);
  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    if (lib[i] < 1 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(source_std[lib[i] - 1]) &&
        !std::isnan(target_std[lib[i] - 1]) &&
        (lib[i] > max_lag + 1)) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(source_std[pred[i] - 1]) &&
        !std::isnan(target_std[pred[i] - 1]) &&
        (pred[i] > max_lag + 1)) {
      pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
    }
  }

  std::vector<std::vector<double>> res_std = SMap4TS(
    source_std,
    target_std,
    lib_indices,
    pred_indices,
    theta_std,
    E,
    tau,
    b,
    threads);

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("theta", "rho", "mae", "rmse");
  return result;
}

//  Wrapper function to help determining embedding dimension `E` and numbers of neighbors `k` parameters
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppMultiSimplex4TS(const Rcpp::NumericMatrix& source,
                                        const Rcpp::NumericMatrix& target,
                                        const Rcpp::IntegerVector& lib,
                                        const Rcpp::IntegerVector& pred,
                                        const Rcpp::IntegerVector& E,
                                        const Rcpp::IntegerVector& b,
                                        int tau,
                                        int threads) {
  // Convert Rcpp NumericMatrix to std::vector of std::vectors
  std::vector<std::vector<double>> source_std(source.ncol());
  std::vector<std::vector<double>> target_std(target.ncol());
  for (int i = 0; i < source.ncol(); ++i) {
    Rcpp::NumericVector covvar = source.column(i);
    source_std[i] = Rcpp::as<std::vector<double>>(covvar);
  }
  for (int i = 0; i < target.ncol(); ++i) {
    Rcpp::NumericVector covvar = target.column(i);
    target_std[i] = Rcpp::as<std::vector<double>>(covvar);
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  std::vector<int> pred_indices;

  int target_len = target.nrow();
  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    if (lib[i] < 1 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }

    bool allnotnan = true;
    for (size_t j = 0; j < target_std.size(); ++j){
      if (std::isnan(source_std[j][lib[i] - 1]) ||
          std::isnan(target_std[j][lib[i] - 1])){
        allnotnan = false;
      }
    }

    if (allnotnan) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }

  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }

    bool allnotnan = true;
    for (size_t j = 0; j < target_std.size(); ++j){
      if (std::isnan(source_std[j][pred[i] - 1]) ||
          std::isnan(target_std[j][pred[i] - 1])){
        allnotnan = false;
      }
    }

    if (allnotnan) {
      pred_indices.push_back(pred[i] - 1); // Convert to 0-based index
    }
  }

  std::vector<std::vector<double>> res_std = MultiSimplex4TS(
    source_std,
    target_std,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau,
    threads);

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "rho", "mae", "rmse");
  return result;
}

// Wrapper function to compute intersection cardinality for time series data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppIC4TS(const Rcpp::NumericVector& source,
                              const Rcpp::NumericVector& target,
                              const Rcpp::IntegerVector& lib,
                              const Rcpp::IntegerVector& pred,
                              const Rcpp::IntegerVector& E,
                              const Rcpp::IntegerVector& b,
                              int tau,
                              int exclude = 0,
                              int threads = 8,
                              int parallel_level = 0) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> source_std = Rcpp::as<std::vector<double>>(source);
  std::vector<double> target_std = Rcpp::as<std::vector<double>>(target);

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);

  // Initialize lib_indices and pred_indices
  std::vector<size_t> lib_indices;
  std::vector<size_t> pred_indices;

  int target_len = target_std.size();
  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  size_t n_libsize = lib.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_libsize; ++i) {
    if (lib[i] < 1 || lib[i] > target_len) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(source_std[lib[i] - 1]) && !std::isnan(target_std[lib[i] - 1])) {
      lib_indices.push_back(static_cast<size_t>(lib[i] - 1)); // Convert to 0-based index
    }
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(source_std[pred[i] - 1]) && !std::isnan(target_std[pred[i] - 1])) {
      pred_indices.push_back(static_cast<size_t>(pred[i] - 1)); // Convert to 0-based index
    }
  }

  // Check the validity of the neignbor numbers
  std::vector<int> b_std;
  for (int i = 0; i < b.size(); ++i){
    if (b[i] > static_cast<int>(lib_indices.size())) {
      Rcpp::stop("Neighbor numbers count out of acceptable range at position %d (value: %d)", i + 1, b[i]);
    }
    b_std.push_back(b[i]);
  }

  std::vector<std::vector<double>> res_std = IC4TS(
    source_std,
    target_std,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau,
    exclude,
    threads,
    parallel_level);

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "CausalScore", "Significance");
  return result;
}

// Wrapper function to perform convergent cross mapping for time series data
// predict y based on x ====> x xmap y ====> y causes x
// [[Rcpp::export(rng = false)]]
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
  int max_lag = (tau == 0) ? (E - 1) : (E * tau);
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(x_std[lib[i] - 1]) &&
        !std::isnan(y_std[lib[i] - 1]) &&
        (lib[i] > max_lag + 1)) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(x_std[pred[i] - 1]) &&
        !std::isnan(y_std[pred[i] - 1]) &&
        (pred[i] > max_lag + 1)) {
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

// Wrapper function to perform partial cross mapping for time series data
// predict y based on x ====> x xmap y ====> y causes x (account for controls)
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppPCM(const Rcpp::NumericVector& x,
                            const Rcpp::NumericVector& y,
                            const Rcpp::NumericMatrix& z,
                            const Rcpp::IntegerVector& libsizes,
                            const Rcpp::IntegerVector& lib,
                            const Rcpp::IntegerVector& pred,
                            const Rcpp::IntegerVector& E,
                            const Rcpp::IntegerVector& tau,
                            const Rcpp::IntegerVector& b,
                            bool simplex,
                            double theta,
                            int threads,
                            int parallel_level,
                            bool cumulate,
                            bool progressbar) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp NumericMatrix to std::vector of std::vectors
  std::vector<std::vector<double>> z_std(z.ncol());
  for (int i = 0; i < z.ncol(); ++i) {
    Rcpp::NumericVector covvar = z.column(i);
    z_std[i] = Rcpp::as<std::vector<double>>(covvar);
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);

  // Convert and check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  std::vector<int> lib_std;
  std::vector<int> pred_std;
  int max_E = *std::max_element(E_std.begin(), E_std.end());
  int max_tau = *std::max_element(tau_std.begin(), tau_std.end());
  int max_lag = (max_tau == 0) ? (max_E - 1) : (max_E * max_tau);
  int n = y_std.size();
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > n) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(x_std[lib[i] - 1]) &&
        !std::isnan(y_std[lib[i] - 1]) &&
        (lib[i] > max_lag + 1)) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(x_std[pred[i] - 1]) &&
        !std::isnan(y_std[pred[i] - 1]) &&
        (pred[i] > max_lag + 1)) {
      pred_std.push_back(pred[i] - 1);
    }
  }

  // Perform partial cross mapping
  std::vector<std::vector<double>> result = PCM(
    x_std,
    y_std,
    z_std,
    libsizes_std,
    lib_std,
    pred_std,
    E_std,
    tau_std,
    b_std,
    simplex,
    theta,
    threads,
    parallel_level,
    cumulate,
    progressbar);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  Rcpp::NumericMatrix resultMatrix(result.size(), 9);
  for (size_t i = 0; i < result.size(); ++i) {
    resultMatrix(i, 0) = result[i][0];
    resultMatrix(i, 1) = result[i][1];
    resultMatrix(i, 2) = result[i][2];
    resultMatrix(i, 3) = result[i][3];
    resultMatrix(i, 4) = result[i][4];
    resultMatrix(i, 5) = result[i][5];
    resultMatrix(i, 6) = result[i][6];
    resultMatrix(i, 7) = result[i][7];
    resultMatrix(i, 8) = result[i][8];
  }

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create(
    "libsizes","T_mean","D_mean",
    "T_sig","T_upper","T_lower",
    "D_sig","D_upper","D_lower");
  return resultMatrix;
}

// Wrapper function to perform cross mapping cardinality for time series data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppCMC(
    const Rcpp::NumericVector& x,
    const Rcpp::NumericVector& y,
    const Rcpp::IntegerVector& libsizes,
    const Rcpp::IntegerVector& lib,
    const Rcpp::IntegerVector& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    int b,
    int r,
    int threads,
    int parallel_level,
    bool progressbar){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x_std = Rcpp::as<std::vector<double>>(x);
  std::vector<double> y_std = Rcpp::as<std::vector<double>>(y);

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<size_t> libsizes_std = Rcpp::as<std::vector<size_t>>(libsizes);
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);

  int validSampleNum = x_std.size();
  // Convert and check that lib and pred indices are within bounds & convert R based 1 index to C++ based 0 index
  std::vector<size_t> lib_std;
  std::vector<size_t> pred_std;
  for (int i = 0; i < lib.size(); ++i) {
    if (lib[i] < 1 || lib[i] > validSampleNum) {
      Rcpp::stop("lib contains out-of-bounds index at position %d (value: %d)", i + 1, lib[i]);
    }
    if (!std::isnan(x_std[lib[i] - 1]) && !std::isnan(y_std[lib[i] - 1])) {
      lib_std.push_back(static_cast<size_t>(lib[i] - 1));
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > validSampleNum) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(x_std[pred[i] - 1]) && !std::isnan(y_std[pred[i] - 1])) {
      pred_std.push_back(static_cast<size_t>(pred[i] - 1));
    }
  }

  // check b that are greater than validSampleNum or less than or equal to 3
  if (b < 3 || b > validSampleNum) {
    Rcpp::stop("k cannot be less than or equal to 3 or greater than the number of non-NA values.");
  } else if (b + 1 > static_cast<int>(lib_std.size())){
    Rcpp::stop("Please check `libsizes` or `lib`; no valid libraries available for running GCMC.");
  }

  // Generate embeddings
  std::vector<std::vector<double>> e1 = Embed(x_std, E[0], tau_std[0]);
  std::vector<std::vector<double>> e2 = Embed(y_std, E[1], tau_std[1]);

  // Perform CMC for time series data
  CMCRes res = CMC(e1,e2,libsizes_std,lib_std,pred_std,
                   static_cast<size_t>(b),static_cast<size_t>(r),
                   threads,parallel_level,progressbar);

  // Convert mean_aucs to Rcpp::DataFrame
  std::vector<double> libs, aucs;
  for (const auto& cm : res.causal_strength) {
    libs.push_back(cm[0]);
    aucs.push_back(cm[1]);
  }

  Rcpp::DataFrame cs_df = Rcpp::DataFrame::create(
    Rcpp::Named("libsizes") = libs,
    Rcpp::Named("x_xmap_y_mean") = aucs
  );

  // Wrap causal_strength with names
  Rcpp::DataFrame xmap_df = Rcpp::DataFrame::create(
    Rcpp::Named("neighbors") = res.cross_mapping[0],
    Rcpp::Named("x_xmap_y_mean") = res.cross_mapping[1],
    Rcpp::Named("x_xmap_y_sig") = res.cross_mapping[2],
    Rcpp::Named("x_xmap_y_upper") = res.cross_mapping[3],
    Rcpp::Named("x_xmap_y_lower")  = res.cross_mapping[4]
  );

  return Rcpp::List::create(
    Rcpp::Named("xmap") = xmap_df,
    Rcpp::Named("cs") = cs_df
  );
}

// Wrapper function to perform multispatial convergent cross mapping for time series data
// predict y based on x ====> x xmap y ====> y causes x
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppMultispatialCCM(const Rcpp::NumericMatrix& x,
                                        const Rcpp::NumericMatrix& y,
                                        const Rcpp::IntegerVector& libsizes,
                                        int E,
                                        int tau,
                                        int b,
                                        int boot,
                                        int seed,
                                        int threads,
                                        int parallel_level,
                                        bool progressbar) {
  // Convert Rcpp NumericMatrix to std::vector of std::vectors
  std::vector<std::vector<double>> x_std(x.ncol());
  std::vector<std::vector<double>> y_std(y.ncol());
  for (int i = 0; i < x.ncol(); ++i) {
    Rcpp::NumericVector covvar = x.column(i);
    x_std[i] = Rcpp::as<std::vector<double>>(covvar);
  }
  for (int i = 0; i < y.ncol(); ++i) {
    Rcpp::NumericVector covvar = y.column(i);
    y_std[i] = Rcpp::as<std::vector<double>>(covvar);
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> libsizes_std = Rcpp::as<std::vector<int>>(libsizes);


  // Perform GCCM Lattice
  std::vector<std::vector<double>> result = MultispatialCCM(
    x_std,
    y_std,
    libsizes_std,
    E,
    tau,
    b,
    boot,
    threads,
    std::abs(seed),
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
