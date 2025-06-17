#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include "CppStats.h"
#include "Embed.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>

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
