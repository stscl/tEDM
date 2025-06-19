#include <vector>
#include <RcppThread.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(rng = false)]]
unsigned int DetectMaxNumThreads(){
  unsigned int max_threads = std::thread::hardware_concurrency();
  return max_threads;
}

/**
 * Determine the optimal embedding dimension (E) and number of nearest neighbors (k).
 *
 * This function takes a matrix `Emat` with columns "E", "k", "rho", "mae", and "rmse".
 * It selects the optialmal embedding dimension (E) and number of nearest neighbors (k)
 * by first maximizing "rho", then minimizing "rmse", and finally minimizing "mae" if necessary.
 *
 * @param Emat A NumericMatrix with five columns: "E", k", "rho", "mae", and "rmse".
 * @return The optimal embedding dimension (E) and number of nearest neighbors (k) as an integer vector.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptEmbedDim(Rcpp::NumericMatrix Emat) {
  // Check if the input matrix has exactly 5 columns
  if (Emat.ncol() != 5) {
    Rcpp::stop("Input matrix must have exactly 5 columns: E, k, rho, mae, and rmse.");
  }

  // Initialize variables to store the optialmal row index and its metrics
  int optialmal_row = 0;
  double optialmal_rho = Emat(0, 2); // Initialize with the first row's rho
  double optialmal_rmse = Emat(0, 4); // Initialize with the first row's rmse
  double optialmal_mae = Emat(0, 3); // Initialize with the first row's mae

  // Iterate through each row of the matrix
  for (int i = 1; i < Emat.nrow(); ++i) {
    double current_rho = Emat(i, 2); // Current row's rho
    double current_rmse = Emat(i, 4); // Current row's rmse
    double current_mae = Emat(i, 3); // Current row's mae

    // Compare rho values first
    if (current_rho > optialmal_rho) {
      optialmal_row = i;
      optialmal_rho = current_rho;
      optialmal_rmse = current_rmse;
      optialmal_mae = current_mae;
    } else if (current_rho == optialmal_rho) {
      // If rho is equal, compare rmse values
      if (current_rmse < optialmal_rmse) {
        optialmal_row = i;
        optialmal_rho = current_rho;
        optialmal_rmse = current_rmse;
        optialmal_mae = current_mae;
      } else if (current_rmse == optialmal_rmse) {
        // If rmse is also equal, compare mae values
        if (current_mae < optialmal_mae) {
          optialmal_row = i;
          optialmal_rho = current_rho;
          optialmal_rmse = current_rmse;
          optialmal_mae = current_mae;
        }
      }
    }
  }

  // Return the optimal E and b from the optialmal row
  Rcpp::IntegerVector result(2);
  result[0] = static_cast<int>(Emat(optialmal_row, 0));
  result[1] = static_cast<int>(Emat(optialmal_row, 1));

  return result;
}

/**
 * Determine the optimal theta parameter based on the evaluation metrics.
 *
 * This function takes a matrix `Thetamat` with columns "theta", "rho", "mae", and "rmse".
 * It selects the optimal theta parameter by first maximizing "rho",
 * then minimizing "rmse", and finally minimizing "mae" if necessary.
 *
 * @param Thetamat A NumericMatrix with four columns: "theta", "rho", "mae", and "rmse".
 * @return The optimal theta parameter as a double.
 */
// [[Rcpp::export(rng = false)]]
double OptThetaParm(Rcpp::NumericMatrix Thetamat) {
  // Check if the input matrix has exactly 4 columns
  if (Thetamat.ncol() != 4) {
    Rcpp::stop("Input matrix must have exactly 4 columns: theta, rho, mae, and rmse.");
  }

  // Initialize variables to store the optialmal row index and its metrics
  int optialmal_row = 0;
  double optialmal_rho = Thetamat(0, 1);  // Initialize with the first row's rho
  double optialmal_rmse = Thetamat(0, 3); // Initialize with the first row's rmse
  double optialmal_mae = Thetamat(0, 2);  // Initialize with the first row's mae

  // Iterate through each row of the matrix
  for (int i = 1; i < Thetamat.nrow(); ++i) {
    double current_rho = Thetamat(i, 1);   // Current row's rho
    double current_rmse = Thetamat(i, 3);  // Current row's rmse
    double current_mae = Thetamat(i, 2);   // Current row's mae

    // Compare rho values first
    if (current_rho > optialmal_rho) {
      optialmal_row = i;
      optialmal_rho = current_rho;
      optialmal_rmse = current_rmse;
      optialmal_mae = current_mae;
    } else if (current_rho == optialmal_rho) {
      // If rho is equal, compare rmse values
      if (current_rmse < optialmal_rmse) {
        optialmal_row = i;
        optialmal_rho = current_rho;
        optialmal_rmse = current_rmse;
        optialmal_mae = current_mae;
      } else if (current_rmse == optialmal_rmse) {
        // If rmse is also equal, compare mae values
        if (current_mae < optialmal_mae) {
          optialmal_row = i;
          optialmal_rho = current_rho;
          optialmal_rmse = current_rmse;
          optialmal_mae = current_mae;
        }
      }
    }
  }

  // Return the optimal theta param from the optialmal row
  return Thetamat(optialmal_row, 0);
}

/**
 * This function takes a NumericMatrix as input and returns a matrix
 * containing the row and column indices of all non-NA elements in the input matrix.
 *
 * The processing order can be controlled using the `byrow` parameter:
 *   - If `byrow` is true, the matrix is processed row by row.
 *   - If `byrow` is false, the matrix is processed column by column.
 *
 * Parameters:
 *   - mat: A NumericMatrix object that is to be processed.
 *   - byrow: A boolean parameter to control the processing order.
 *     - If true, the matrix is processed row by row (default is true).
 *     - If false, the matrix is processed column by column.
 *
 * Returns:
 *   - A NumericMatrix with two columns:
 *     - The first column contains the row indices,
 *     - The second column contains the column indices of non-NA elements.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix MatNotNAIndice(Rcpp::NumericMatrix mat, bool byrow = true) {
  // Initialize vectors to store the row and column indices of non-NA elements
  std::vector<double> row_indices;
  std::vector<double> col_indices;

  // Get the number of rows and columns in the input matrix
  int nrow = mat.nrow();
  int ncol = mat.ncol();

  // Loop through the matrix depending on the value of 'byrow'
  if (byrow) {
    // Process by row (row-wise iteration)
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        // Check if the element is not NA
        if (!Rcpp::NumericMatrix::is_na(mat(i, j))) {
          // Record the row and column indices (1-based indexing for R compatibility)
          row_indices.push_back(i + 1);
          col_indices.push_back(j + 1);
        }
      }
    }
  } else {
    // Process by column (column-wise iteration)
    for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < nrow; i++) {
        // Check if the element is not NA
        if (!Rcpp::NumericMatrix::is_na(mat(i, j))) {
          // Record the row and column indices (1-based indexing for R compatibility)
          row_indices.push_back(i + 1);
          col_indices.push_back(j + 1);
        }
      }
    }
  }

  // Create a NumericMatrix to store the result
  int n = row_indices.size();
  Rcpp::NumericMatrix result(n, 2);

  // Fill the result matrix with the row and column indices
  for (int i = 0; i < n; i++) {
    result(i, 0) = row_indices[i];
    result(i, 1) = col_indices[i];
  }

  // Return the result matrix
  return result;
}
