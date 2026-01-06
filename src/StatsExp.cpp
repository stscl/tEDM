#include <cmath>
#include <vector>
#include <string>
#include "CppStats.h"
#include "CppCombn.h"
#include "CppDistances.h"
#include "DeLongPlacements.h"
#include "SpatialBlockBootstrap.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(rng = false)]]
int RcppFactorial(int n){
  return(CppFactorial(n));
};

// [[Rcpp::export(rng = false)]]
double RcppCombine(int n,int k){
  return(CppCombine(n,k));
};

// [[Rcpp::export(rng = false)]]
Rcpp::List RcppCombn(const Rcpp::RObject& vec, int m) {
  if (TYPEOF(vec) == REALSXP) {
    std::vector<double> input = Rcpp::as<std::vector<double>>(vec);
    return Rcpp::wrap(CppCombn(input, m));  // Calls the double version of the template
  } else if (TYPEOF(vec) == INTSXP) {
    std::vector<int> input = Rcpp::as<std::vector<int>>(vec);
    return Rcpp::wrap(CppCombn(input, m));  // Calls the int version of the template
  } else if (TYPEOF(vec) == STRSXP) {
    std::vector<std::string> input = Rcpp::as<std::vector<std::string>>(vec);
    return Rcpp::wrap(CppCombn(input, m));  // Calls the string version of the template
  } else {
    Rcpp::stop("Unsupported vector type. Must be numeric, integer, or character.");
  }
}

// [[Rcpp::export(rng = false)]]
Rcpp::List RcppGenSubsets(const Rcpp::RObject& vec) {
  if (TYPEOF(vec) == REALSXP) {
    std::vector<double> input = Rcpp::as<std::vector<double>>(vec);
    return Rcpp::wrap(CppGenSubsets(input));  // Calls the double version of the template
  } else if (TYPEOF(vec) == INTSXP) {
    std::vector<int> input = Rcpp::as<std::vector<int>>(vec);
    return Rcpp::wrap(CppGenSubsets(input));  // Calls the int version of the template
  } else if (TYPEOF(vec) == STRSXP) {
    std::vector<std::string> input = Rcpp::as<std::vector<std::string>>(vec);
    return Rcpp::wrap(CppGenSubsets(input));  // Calls the string version of the template
  } else {
    Rcpp::stop("Unsupported vector type. Must be numeric, integer, or character.");
  }
}

// [[Rcpp::export(rng = false)]]
double RcppDigamma(double x){
  return(CppDigamma(x));
};

// [[Rcpp::export(rng = false)]]
double RcppLog(double x, double base = 10){
  return(CppLog(x, base));
};

// [[Rcpp::export(rng = false)]]
double RcppMedian(const Rcpp::NumericVector& vec,
                  bool NA_rm = false) {
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);
  return CppMedian(y, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppMean(const Rcpp::NumericVector& vec,
                bool NA_rm = false) {
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);
  return CppMean(y, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppMin(const Rcpp::NumericVector& vec,
               bool NA_rm = false) {
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);
  return CppMin(y, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppMax(const Rcpp::NumericVector& vec,
               bool NA_rm = false) {
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);
  return CppMax(y, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppSum(const Rcpp::NumericVector& vec,
               bool NA_rm = false) {
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);
  return CppSum(y, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppVariance(const Rcpp::NumericVector& vec,
                    bool NA_rm = false) {
  std::vector<double> y = Rcpp::as<std::vector<double>>(vec);
  return CppVariance(y, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppCovariance(const Rcpp::NumericVector& vec1,
                      const Rcpp::NumericVector& vec2,
                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x1_vec = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> x2_vec = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppCovariance function
  return CppCovariance(x1_vec, x2_vec, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppMAE(const Rcpp::NumericVector& vec1,
               const Rcpp::NumericVector& vec2,
               bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x1_vec = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> x2_vec = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppMAE function
  return CppMAE(x1_vec, x2_vec, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppRMSE(const Rcpp::NumericVector& vec1,
                const Rcpp::NumericVector& vec2,
                bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> x1_vec = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> x2_vec = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppRMSE function
  return CppRMSE(x1_vec, x2_vec, NA_rm);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppCumSum(const Rcpp::NumericVector& vec) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);

  // Call the CppCumSum function
  std::vector<double> result = CppCumSum(vec_std);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppAbsDiff(const Rcpp::NumericVector& vec1,
                                const Rcpp::NumericVector& vec2) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec1_cpp = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> vec2_cpp = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppAbsDiff function
  std::vector<double> result = CppAbsDiff(vec1_cpp, vec2_cpp);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSumNormalize(const Rcpp::NumericVector& vec,
                                     bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_cpp = Rcpp::as<std::vector<double>>(vec);

  // Call the CppSumNormalize function
  std::vector<double> result = CppSumNormalize(vec_cpp, NA_rm);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppArithmeticSeq(double from, double to, int length_out) {
  // Call the CppArithmeticSeq function
  std::vector<double> result = CppArithmeticSeq(from, to, static_cast<size_t>(length_out));
  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppQuantile(const Rcpp::NumericVector& vec,
                                 const Rcpp::NumericVector& probs,
                                 bool NA_rm = true) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);
  std::vector<double> probs_std = Rcpp::as<std::vector<double>>(probs);

  // Call the CppQuantilefunction
  std::vector<double> result = CppQuantile(vec_std, probs_std, NA_rm);

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// [[Rcpp::export(rng = false)]]
double RcppPearsonCor(const Rcpp::NumericVector& y,
                      const Rcpp::NumericVector& y_hat,
                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y_vec = Rcpp::as<std::vector<double>>(y);
  std::vector<double> y_hat_vec = Rcpp::as<std::vector<double>>(y_hat);

  // Call the PearsonCor function
  return PearsonCor(y_vec, y_hat_vec, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppSpearmanCor(const Rcpp::NumericVector& y,
                       const Rcpp::NumericVector& y_hat,
                       bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y_vec = Rcpp::as<std::vector<double>>(y);
  std::vector<double> y_hat_vec = Rcpp::as<std::vector<double>>(y_hat);

  // Call the SpearmanCorfunction
  return SpearmanCor(y_vec, y_hat_vec, NA_rm);
}

// [[Rcpp::export(rng = false)]]
double RcppKendallCor(const Rcpp::NumericVector& y,
                      const Rcpp::NumericVector& y_hat,
                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> y_vec = Rcpp::as<std::vector<double>>(y);
  std::vector<double> y_hat_vec = Rcpp::as<std::vector<double>>(y_hat);

  // Call the KendallCor function
  return KendallCor(y_vec, y_hat_vec, NA_rm);
}

// Rcpp wrapper for PartialCor function
// [[Rcpp::export(rng = false)]]
double RcppPartialCor(const Rcpp::NumericVector& y,
                      const Rcpp::NumericVector& y_hat,
                      const Rcpp::NumericMatrix& controls,
                      bool NA_rm = false,
                      bool linear = false,
                      double pinv_tol = 1e-10) {

  // Convert Rcpp NumericVector to std::vector
  std::vector<double> std_y = Rcpp::as<std::vector<double>>(y);
  std::vector<double> std_y_hat = Rcpp::as<std::vector<double>>(y_hat);

  // Convert Rcpp NumericMatrix to std::vector of std::vectors
  std::vector<std::vector<double>> std_controls(controls.ncol());
  for (int i = 0; i < controls.ncol(); ++i) {
    Rcpp::NumericVector covvar = controls.column(i);
    std_controls[i] = Rcpp::as<std::vector<double>>(covvar);
  }

  // Call the PartialCor function
  return PartialCor(std_y, std_y_hat, std_controls, NA_rm, linear);
}

// [[Rcpp::export(rng = false)]]
double RcppPartialCorTrivar(const Rcpp::NumericVector& y,
                            const Rcpp::NumericVector& y_hat,
                            const Rcpp::NumericVector& control,
                            bool NA_rm = false,
                            bool linear = false,
                            double pinv_tol = 1e-10) {
  // Convert Rcpp NumericVector to std::vector
  std::vector<double> std_y = Rcpp::as<std::vector<double>>(y);
  std::vector<double> std_y_hat = Rcpp::as<std::vector<double>>(y_hat);
  std::vector<double> std_control = Rcpp::as<std::vector<double>>(control);

  // Call the PartialCorTrivar function
  return PartialCorTrivar(std_y, std_y_hat, std_control, NA_rm, linear, pinv_tol);
}

// Wrapper function to calculate the significance of a (partial) correlation coefficient
// [[Rcpp::export(rng = false)]]
double RcppCorSignificance(double r, int n, int k = 0){
  return CppCorSignificance(r, static_cast<size_t>(n), static_cast<size_t>(k));
}

// Wrapper function to calculate the confidence interval for a (partial) correlation coefficient
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppCorConfidence(double r, int n, int k = 0,
                                      double level = 0.05) {
  // Calculate the confidence interval
  std::vector<double> result = CppCorConfidence(r,
                                                static_cast<size_t>(n),
                                                static_cast<size_t>(k),
                                                level);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Wrapper function to calculate the significance of a vector of (partial) correlation coefficients
// [[Rcpp::export(rng = false)]]
double RcppMeanCorSignificance(const Rcpp::NumericVector& r, int n, int k = 0){
  // Convert Rcpp inputs to standard C++ types
  std::vector<double> r_std = Rcpp::as<std::vector<double>>(r);
  return CppMeanCorSignificance(r_std, static_cast<size_t>(n), static_cast<size_t>(k));
}

// Wrapper function to calculate the confidence interval for a vector of (partial) correlation coefficients
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppMeanCorConfidence(const Rcpp::NumericVector& r,
                                          int n, int k = 0,
                                          double level = 0.05) {
  // Convert Rcpp inputs to standard C++ types
  std::vector<double> r_std = Rcpp::as<std::vector<double>>(r);
  // Calculate the confidence interval
  std::vector<double> result = CppMeanCorConfidence(r_std,
                                                    static_cast<size_t>(n),
                                                    static_cast<size_t>(k),
                                                    level);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Wrapper function to performs delong's test for ROC AUC comparison and return a NumericVector
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppDeLongAUCConfidence(const Rcpp::NumericVector& cases,
                                            const Rcpp::NumericVector& controls,
                                            const Rcpp::CharacterVector& direction,
                                            double level = 0.05) {
  // Convert Rcpp inputs to standard C++ types
  std::vector<double> cpp_cases = Rcpp::as<std::vector<double>>(cases);
  std::vector<double> cpp_controls = Rcpp::as<std::vector<double>>(controls);
  std::string cpp_direction = Rcpp::as<std::string>(direction[0]);

  // Call the CppDeLongAUCConfidence function
  std::vector<double> result = CppDeLongAUCConfidence(cpp_cases, cpp_controls, cpp_direction, level);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Wrapper function to performs delong's test for CMC causal score and return a NumericVector
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppCMCTest(const Rcpp::NumericVector& cases,
                                const Rcpp::CharacterVector& direction,
                                double level = 0.05,
                                int num_samples = 0) {
  // Convert Rcpp inputs to standard C++ types
  std::vector<double> cpp_cases = Rcpp::as<std::vector<double>>(cases);
  std::string cpp_direction = Rcpp::as<std::string>(direction[0]);

  // Call the CppCMCTest function
  std::vector<double> result = CppCMCTest(cpp_cases, cpp_direction, level, static_cast<size_t>(num_samples));

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Wrapper function to compute distance between two vectors
// [[Rcpp::export(rng = false)]]
double RcppDistance(const Rcpp::NumericVector& vec1,
                    const Rcpp::NumericVector& vec2,
                    bool L1norm = false,
                    bool NA_rm = false){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> v1 = Rcpp::as<std::vector<double>>(vec1);
  std::vector<double> v2 = Rcpp::as<std::vector<double>>(vec2);

  // Call the CppDistance function
  return CppDistance(v1, v2, L1norm ,NA_rm);
}

// Wrapper function to compute the k-th nearest distance for a vector.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppKNearestDistance(const Rcpp::NumericVector& vec1,
                                         int k,
                                         bool L1norm = false,
                                         bool NA_rm = false){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> v1 = Rcpp::as<std::vector<double>>(vec1);

  // Call the CppKNearestDistance function
  std::vector<double> res = CppKNearestDistance(v1, static_cast<size_t>(std::abs(k)), L1norm ,NA_rm);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(res);
}

// Wrapper function to compute the distance matrix of a given matrix 'mat'
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppMatDistance(const Rcpp::NumericMatrix& mat,
                                    bool L1norm = false,
                                    bool NA_rm = false) {

  // Convert the Rcpp::NumericMatrix to a C++ vector of vectors (std::vector)
  size_t rownum = mat.nrow();
  size_t colnum = mat.ncol();
  std::vector<std::vector<double>> cppMat(rownum, std::vector<double>(colnum));

  // Fill cppMat with values from the R matrix
  for (size_t i = 0; i < rownum; ++i) {
    for (size_t j = 0; j < colnum; ++j) {
      cppMat[i][j] = mat(i, j);
    }
  }

  // Call the C++ function to compute the distance matrix
  std::vector<std::vector<double>> distanceMatrix = CppMatDistance(cppMat, L1norm, NA_rm);

  // Convert the resulting distance matrix back into an Rcpp::NumericMatrix
  Rcpp::NumericMatrix result(rownum, rownum);
  for (size_t i = 0; i < rownum; ++i) {
    for (size_t j = 0; j < rownum; ++j) {
      result(i, j) = distanceMatrix[i][j];
    }
  }

  return result;
}

// Wrapper function to compute the number of neighbors for each point within a given radius.
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector RcppNeighborsNum(
    const Rcpp::NumericVector& vec,
    const Rcpp::NumericVector& radius,
    bool equal = false,
    bool L1norm = false,
    bool NA_rm = false){
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> v = Rcpp::as<std::vector<double>>(vec);
  std::vector<double> r = Rcpp::as<std::vector<double>>(radius);

  // Call the CppKNearestDistance function
  std::vector<int> res = CppNeighborsNum(v,r,equal,L1norm,NA_rm);

  // Convert std::vector<int> to Rcpp::IntegerVector
  return Rcpp::wrap(res);
}

// Wrapper function to find k-nearest neighbors of a given index in the embedding space
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector RcppKNNIndice(const Rcpp::NumericMatrix& embedding_space,
                                  int target_idx,
                                  int k,
                                  const Rcpp::IntegerVector& lib) {
  // Get the number of rows and columns
  size_t n_rows = embedding_space.nrow();
  size_t n_cols = embedding_space.ncol();

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> embedding_vec(n_rows, std::vector<double>(n_cols));
  for (std::size_t i = 0; i < n_rows; ++i) {
    for (std::size_t j = 0; j < n_cols; ++j) {
      embedding_vec[i][j] = embedding_space(i, j);
    }
  }

  // Ensure target index is within valid range
  if (target_idx < 0 || static_cast<size_t>(target_idx) >= n_rows) {
    Rcpp::stop("target_idx is out of range.");
  }

  // Ensure k is positive
  if (k <= 0) {
    Rcpp::stop("k must be greater than 0.");
  }

  // Convert lib(1-based R index) to lib_std (0-based C++ index)
  std::vector<int> lib_std;
  size_t n_libsize = lib.size();
  for (size_t i = 0; i < n_libsize; ++i) {
    lib_std.push_back(lib[i] - 1); // Convert to 0-based index
  }

  // Call the C++ function
  std::vector<size_t> knn_indices = CppKNNIndice(embedding_vec,
                                                 static_cast<size_t>(target_idx - 1),
                                                 static_cast<size_t>(k),
                                                 lib_std);

  // Convert result to Rcpp::IntegerVector (R uses 1-based indexing)
  Rcpp::IntegerVector result(knn_indices.size());
  for (size_t i = 0; i < knn_indices.size(); ++i) {
    result[i] = static_cast<int>(knn_indices[i]) + 1;  // Convert to 1-based index
  }

  return result;
}

// Wrapper function to find k-nearest neighbors of a given index using a precomputed distance matrix
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector RcppDistKNNIndice(const Rcpp::NumericMatrix& dist_mat,
                                      int target_idx,
                                      int k,
                                      const Rcpp::IntegerVector& lib) {
  // Get the number of rows and columns
  size_t n_rows = dist_mat.nrow();
  size_t n_cols = dist_mat.ncol();

  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> distmat(n_rows, std::vector<double>(n_cols));
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      distmat[i][j] = dist_mat(i, j);
    }
  }

  // Ensure target index is within valid range
  if (target_idx < 0 || static_cast<size_t>(target_idx) >= n_rows) {
    Rcpp::stop("target_idx is out of range.");
  }

  // Ensure k is positive
  if (k <= 0) {
    Rcpp::stop("k must be greater than 0.");
  }

  // Convert lib(1-based R index) to lib_std (0-based C++ index)
  std::vector<int> lib_std;
  size_t n_libsize = lib.size();
  for (size_t i = 0; i < n_libsize; ++i) {
    lib_std.push_back(lib[i] - 1); // Convert to 0-based index
  }

  // Call the C++ function
  std::vector<size_t> knn_indices = CppDistKNNIndice(distmat,
                                                     static_cast<size_t>(target_idx - 1),
                                                     static_cast<size_t>(k),
                                                     lib_std);

  // Convert result to Rcpp::IntegerVector (R uses 1-based indexing)
  Rcpp::IntegerVector result(knn_indices.size());
  for (size_t i = 0; i < knn_indices.size(); ++i) {
    result[i] = static_cast<int>(knn_indices[i]) + 1;  // Convert to 1-based index
  }

  return result;
}

// Wrapper function to generate sorted neighbor indices
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppDistSortedIndice(const Rcpp::NumericMatrix& dist_mat,
                                const Rcpp::IntegerVector& lib,
                                int k, bool include_self = false) {
  // Get number of rows and columns
  const int n = dist_mat.nrow();
  const int m = dist_mat.ncol();

  // Convert Rcpp data structure to std::vector<>
  std::vector<std::vector<double>> dist_vec(n, std::vector<double>(m));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      dist_vec[i][j] = dist_mat(i, j);
    }
  }

  std::vector<size_t> lib_std(lib.size());
  for (int i = 0; i < lib.size(); ++i) {
    lib_std[i] = static_cast<size_t>(i);
  }

  // Call the existing C++ function to compute sorted neighbor indices
  std::vector<std::vector<size_t>> result = CppDistSortedIndice(dist_vec, lib_std, static_cast<size_t>(k), include_self);

  // Convert the result to an R list of integer vectors
  Rcpp::List out(n);
  for (int i = 0; i < n; ++i) {
    const auto& row = result[i];
    Rcpp::IntegerVector indices(row.size());
    for (size_t j = 0; j < row.size(); ++j) {
      if (row[j] == std::numeric_limits<size_t>::max()) {
        indices[j] = NA_INTEGER;
      } else {
        indices[j] = static_cast<int>(row[j]);
      }
    }
    out[i] = indices;
  }

  // Return the list where each element contains sorted neighbor indices for that row
  return out;
}

// Wrapper function to generate k-nearest neighbors within the embedding space.
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppMatKNNeighbors(const Rcpp::NumericMatrix& embeddings,
                              const Rcpp::IntegerVector& lib,
                              int k, int threads = 8, bool L1norm = false) {
  // Get number of rows and columns
  const int n = embeddings.nrow();
  const int m = embeddings.ncol();

  // Convert Rcpp data structure to std::vector<>
  std::vector<std::vector<double>> emb(n, std::vector<double>(m));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      emb[i][j] = embeddings(i, j);
    }
  }

  std::vector<size_t> lib_std(lib.size());
  for (int i = 0; i < lib.size(); ++i) {
    lib_std[i] = static_cast<size_t>(i);
  }

  // Call the existing C++ function to compute sorted neighbor indices
  std::vector<std::vector<size_t>> result = CppMatKNNeighbors(emb, lib_std, static_cast<size_t>(k), static_cast<size_t>(threads), L1norm);

  // Convert the result to an R list of integer vectors
  Rcpp::List out(n);
  for (int i = 0; i < n; ++i) {
    const auto& row = result[i];
    Rcpp::IntegerVector indices(row.size());
    for (size_t j = 0; j < row.size(); ++j) {
      if (row[j] == std::numeric_limits<size_t>::max()) {
        indices[j] = NA_INTEGER;
      } else {
        indices[j] = static_cast<int>(row[j]);
      }
    }
    out[i] = indices;
  }

  // Return the list where each element contains sorted neighbor indices for that row
  return out;
}

// Wrapper function to perform Linear Trend Removal and return a NumericVector
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppLinearTrendRM(const Rcpp::NumericVector& vec,
                                      const Rcpp::NumericVector& xcoord,
                                      const Rcpp::NumericVector& ycoord,
                                      bool NA_rm = false) {
  // Convert Rcpp::NumericVector to std::vector<double>
  std::vector<double> vec_std = Rcpp::as<std::vector<double>>(vec);
  std::vector<double> xcoord_std = Rcpp::as<std::vector<double>>(xcoord);
  std::vector<double> ycoord_std = Rcpp::as<std::vector<double>>(ycoord);

  // Perform Linear Trend Removal
  std::vector<double> result = LinearTrendRM(vec_std, xcoord_std, ycoord_std, NA_rm);

  // Convert std::vector<double> to Rcpp::NumericVector
  return Rcpp::wrap(result);
}

// Rcpp wrapper function for CppSVD
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppSVD(const Rcpp::NumericMatrix& X) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  size_t m = X.nrow();
  size_t n = X.ncol();
  std::vector<std::vector<double>> X_vec(m, std::vector<double>(n));
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      X_vec[i][j] = X(i, j);
    }
  }

  // Call the original CppSVD function
  std::vector<std::vector<std::vector<double>>> result = CppSVD(X_vec);

  // Extract results from CppSVD output
  std::vector<std::vector<double>> u = result[0]; // Left singular vectors
  std::vector<double> d = result[1][0];           // Singular values
  std::vector<std::vector<double>> v = result[2]; // Right singular vectors

  // Convert std::vector results to Rcpp objects
  Rcpp::NumericMatrix u_rcpp(m, m);
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < m; ++j) {
      u_rcpp(i, j) = u[i][j];
    }
  }

  Rcpp::NumericVector d_rcpp(d.size());
  for (size_t i = 0; i < d.size(); ++i) {
    d_rcpp(i) = d[i];
  }

  Rcpp::NumericMatrix v_rcpp(v.size(), v[0].size());
  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = 0; j < v[0].size(); ++j) {
      v_rcpp(i, j) = v[i][j];
    }
  }

  // Return results as an Rcpp::List to match R's svd() output
  return Rcpp::List::create(
    Rcpp::Named("u") = u_rcpp, // Left singular vectors
    Rcpp::Named("d") = d_rcpp, // Singular values
    Rcpp::Named("v") = v_rcpp  // Right singular vectors
  );
}

// Rcpp wrapper function for CppDeLongPlacements
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppDeLongPlacements(const Rcpp::NumericVector& cases,
                                const Rcpp::NumericVector& controls,
                                const Rcpp::CharacterVector& direction) {
  // Convert Rcpp inputs to standard C++ types
  std::vector<double> cpp_cases = Rcpp::as<std::vector<double>>(cases);
  std::vector<double> cpp_controls = Rcpp::as<std::vector<double>>(controls);
  std::string cpp_direction = Rcpp::as<std::string>(direction[0]);

  // Call the CppDeLongPlacements function
  DeLongPlacementsRes result = CppDeLongPlacements(cpp_cases, cpp_controls, cpp_direction);

  // Return the result as an Rcpp List
  return Rcpp::List::create(
    Rcpp::Named("theta") = result.theta,
    Rcpp::Named("X") = result.X,
    Rcpp::Named("Y") = result.Y
  );
}

// Rcpp wrapper function for SpatialBlockBootstrap
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector RcppSpatialBlockBootstrap(
    const Rcpp::IntegerVector& block,
    unsigned int seed = 42){
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(block);
  std::vector<int> result = SpatialBlockBootstrap(b_std,seed);
  return Rcpp::wrap(result);
}
