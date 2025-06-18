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
#include "CCM.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Wrapper function to perform trivariate logistic map
// [[Rcpp::export(rng=false)]]
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
// [[Rcpp::export(rng=false)]]
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
// [[Rcpp::export(rng=false)]]
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

// Wrapper function to perform s-mapping forecasting
// [[Rcpp::export(rng=false)]]
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

//  Wrapper function to help determining embedding dimension `E` and numbers of neighbors `k` parameters
// [[Rcpp::export(rng=false)]]
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
    if (!std::isnan(target_std[lib[i] - 1])) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(target_std[pred[i] - 1])) {
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
// [[Rcpp::export(rng=false)]]
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
    if (!std::isnan(target_std[lib[i] - 1]) && (lib[i] > max_lag)) {
      lib_indices.push_back(lib[i] - 1); // Convert to 0-based index
    }
  }
  size_t n_predsize = pred.size();   // convert R R_xlen_t to C++ size_t
  for (size_t i = 0; i < n_predsize; ++i) {
    if (pred[i] < 1 || pred[i] > target_len) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(target_std[pred[i] - 1]) && (pred[i] > max_lag)) {
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


// Wrapper function to perform convergent cross mapping for time series data
// predict y based on x ====> x xmap y ====> y causes x
// [[Rcpp::export(rng=false)]]
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
    if (!std::isnan(y_std[lib[i] - 1]) && (lib[i] > max_lag)) {
      lib_std.push_back(lib[i] - 1);
    }
  }
  for (int i = 0; i < pred.size(); ++i) {
    if (pred[i] < 1 || pred[i] > n) {
      Rcpp::stop("pred contains out-of-bounds index at position %d (value: %d)", i + 1, pred[i]);
    }
    if (!std::isnan(y_std[pred[i] - 1])&& (pred[i] > max_lag)) {
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
