#include <vector>
#include <limits>
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
 * This function selects the best (E, k) combination based on:
 *   1. Maximizing rho
 *   2. Minimizing rmse
 *   3. Minimizing mae
 *   4. If still tied, choosing smallest k, then smallest E
 * A warning is issued when tie-breaking by k and E is used.
 *
 * @param Emat A NumericMatrix with 5 columns: E, k, rho, mae, rmse.
 * @return IntegerVector with optimal E and k.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptEmbedDim(Rcpp::NumericMatrix Emat) {
  if (Emat.ncol() != 5) {
    Rcpp::stop("Input matrix must have exactly 5 columns: E, k, rho, mae, and rmse.");
  }

  const double tol = 1e-10;
  int n = Emat.nrow();
  int opt_row = 0;

  double best_rho = Emat(0, 2);
  double best_rmse = Emat(0, 4);
  double best_mae = Emat(0, 3);
  int best_k = static_cast<int>(Emat(0, 1));
  int best_E = static_cast<int>(Emat(0, 0));

  bool used_kE_tiebreak = false;

  for (int i = 1; i < n; ++i) {
    double rho = Emat(i, 2);
    double rmse = Emat(i, 4);
    double mae = Emat(i, 3);
    int k = static_cast<int>(Emat(i, 1));
    int E = static_cast<int>(Emat(i, 0));

    bool rho_better = (rho - best_rho) > tol;
    bool rho_equal = std::abs(rho - best_rho) <= tol;
    bool rmse_better = (best_rmse - rmse) > tol;
    bool rmse_equal = std::abs(rmse - best_rmse) <= tol;
    bool mae_better = (best_mae - mae) > tol;
    bool mae_equal = std::abs(mae - best_mae) <= tol;

    if (rho_better ||
        (rho_equal && rmse_better) ||
        (rho_equal && rmse_equal && mae_better)) {
      opt_row = i;
      best_rho = rho;
      best_rmse = rmse;
      best_mae = mae;
      best_k = k;
      best_E = E;
      used_kE_tiebreak = false;
    } else if (rho_equal && rmse_equal && mae_equal) {
      // Tie on all three metrics: resolve using k then E
      if (k < best_k || (k == best_k && E < best_E)) {
        opt_row = i;
        best_k = k;
        best_E = E;
        used_kE_tiebreak = true;
      }
    }
  }

  if (used_kE_tiebreak) {
    Rcpp::warning("Ties in evaluation metrics resolved by selecting the smallest k, then smallest E.");
  }

  Rcpp::IntegerVector result(2);
  result[0] = static_cast<int>(Emat(opt_row, 0)); // E
  result[1] = static_cast<int>(Emat(opt_row, 1)); // k
  return result;
}

/**
 * Determine the optimal theta parameter based on evaluation metrics.
 *
 * This function takes a NumericMatrix `Thetamat` with columns:
 * "theta", "rho", "mae", and "rmse".
 * The selection criteria are:
 *  - Maximize "rho"
 *  - Minimize "rmse" if "rho" ties
 *  - Minimize "mae" if "rho" and "rmse" tie
 * If multiple rows tie on these metrics (within a tolerance of 1e-10),
 * preference is given to theta == 1, or else to the theta closest to 1.
 * Warnings are issued when tie-breaking occurs or when all metrics are identical.
 *
 * @param Thetamat A NumericMatrix with four columns: theta, rho, mae, and rmse.
 * @return The optimal theta parameter as a double.
 */
// [[Rcpp::export(rng = false)]]
double OptThetaParm(Rcpp::NumericMatrix Thetamat) {
  if (Thetamat.ncol() != 4) {
    Rcpp::stop("Input matrix must have exactly 4 columns: theta, rho, mae, and rmse.");
  }

  int n = Thetamat.nrow();
  std::vector<int> best_rows;
  double best_rho = Thetamat(0, 1);
  double best_rmse = Thetamat(0, 3);
  double best_mae = Thetamat(0, 2);
  const double tol = 1e-10;

  best_rows.push_back(0);

  for (int i = 1; i < n; ++i) {
    double rho = Thetamat(i, 1);
    double rmse = Thetamat(i, 3);
    double mae = Thetamat(i, 2);

    bool rho_better = (rho - best_rho) > tol;
    bool rho_equal = std::abs(rho - best_rho) <= tol;
    bool rmse_better = (best_rmse - rmse) > tol;
    bool rmse_equal = std::abs(rmse - best_rmse) <= tol;
    bool mae_better = (best_mae - mae) > tol;
    bool mae_equal = std::abs(mae - best_mae) <= tol;

    if (rho_better ||
        (rho_equal && rmse_better) ||
        (rho_equal && rmse_equal && mae_better)) {
      best_rows.clear();
      best_rows.push_back(i);
      best_rho = rho;
      best_rmse = rmse;
      best_mae = mae;
    } else if (rho_equal && rmse_equal && mae_equal) {
      best_rows.push_back(i);
    }
  }

  // If only one best row, return its theta
  if (best_rows.size() == 1) {
    return Thetamat(best_rows[0], 0);
  }

  // Check if all metrics are identical (within tolerance)
  bool all_equal = true;
  for (size_t i = 1; i < best_rows.size(); ++i) {
    int r0 = best_rows[0];
    int r1 = best_rows[i];
    if (std::abs(Thetamat(r0, 1) - Thetamat(r1, 1)) > tol ||
        std::abs(Thetamat(r0, 2) - Thetamat(r1, 2)) > tol ||
        std::abs(Thetamat(r0, 3) - Thetamat(r1, 3)) > tol) {
      all_equal = false;
      break;
    }
  }

  // Select theta == 1 if exists, else theta closest to 1
  double selected_theta = std::numeric_limits<double>::quiet_NaN();
  double min_dist_to_1 = std::numeric_limits<double>::max();

  for (int i : best_rows) {
    double theta = Thetamat(i, 0);
    if (std::abs(theta - 1.0) <= tol) {
      selected_theta = theta;
      break;
    } else {
      double dist = std::abs(theta - 1.0);
      if (dist < min_dist_to_1) {
        min_dist_to_1 = dist;
        selected_theta = theta;
      }
    }
  }

  if (all_equal) {
    Rcpp::warning("All evaluation metrics are identical within tolerance; choosing theta == 1 if available, otherwise closest to 1.");
  } else if (best_rows.size() > 1) {
    Rcpp::warning("Tied best evaluation metrics within tolerance; choosing theta == 1 if available, otherwise closest to 1.");
  }

  return selected_theta;
}

/**
 * Select the optimal embedding dimension (E) and number of nearest neighbors (k)
 * from a 4-column matrix: E, k, performance metric, and p-value.
 *
 * Only rows with p-value <= 0.05 are considered.
 * Among them, select the row with:
 *   1. Highest metric (compared using relative tolerance for robustness),
 *   2. If tie, smallest k,
 *   3. If still tie, smallest E.
 *
 * If multiple rows tie on the best metric (within tolerance), a warning is issued
 * and the combination with the smallest k and E is chosen.
 *
 * If no valid rows (p <= 0.05) exist, the function stops with an error.
 *
 * @param Emat NumericMatrix with columns: E, k, metric, and p-value.
 * @return IntegerVector of length 2: optimal E and k.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector OptICparm(Rcpp::NumericMatrix Emat) {
  if (Emat.ncol() != 4) {
    Rcpp::stop("Input matrix must have exactly 4 columns: E, k, metric, and p-value.");
  }

  // Helper: check if two doubles are nearly equal (relative + absolute tolerance)
  auto nearlyEqual = [](double a, double b, double rel_tol = 1e-6, double abs_tol = 1e-12) -> bool {
    return std::abs(a - b) <= std::max(rel_tol * std::max(std::abs(a), std::abs(b)), abs_tol);
  };

  std::vector<int> valid_rows;
  for (int i = 0; i < Emat.nrow(); ++i) {
    if (Emat(i, 3) <= 0.05) {
      valid_rows.push_back(i);
    }
  }

  if (valid_rows.empty()) {
    Rcpp::stop("No valid rows with p-value <= 0.05. The chosen neighborhood parameter may be unreasonable or there's no causal relationship. Please consider resetting.");
  }

  int optimal_row = valid_rows[0];
  double best_metric = Emat(optimal_row, 2);
  int best_k = static_cast<int>(Emat(optimal_row, 1));
  int best_E = static_cast<int>(Emat(optimal_row, 0));
  int tie_count = 1;

  for (size_t i = 1; i < valid_rows.size(); ++i) {
    int row = valid_rows[i];
    double current_metric = Emat(row, 2);
    int current_k = static_cast<int>(Emat(row, 1));
    int current_E = static_cast<int>(Emat(row, 0));

    if (current_metric > best_metric && !nearlyEqual(current_metric, best_metric)) {
      optimal_row = row;
      best_metric = current_metric;
      best_k = current_k;
      best_E = current_E;
      tie_count = 1;
    } else if (nearlyEqual(current_metric, best_metric)) {
      ++tie_count;
      if (current_k < best_k || (current_k == best_k && current_E < best_E)) {
        optimal_row = row;
        best_k = current_k;
        best_E = current_E;
      }
    }
  }

  if (tie_count > 1) {
    Rcpp::warning("Multiple parameter combinations have identical optimal metric (within tolerance); selected one with smallest k and E.");
  }

  Rcpp::IntegerVector result(2);
  result[0] = best_E;
  result[1] = best_k;
  return result;
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
