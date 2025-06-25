#ifndef HelperFuns_H
#define HelperFuns_H

#include <vector>
#include <RcppThread.h>
#include <RcppArmadillo.h>

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
Rcpp::IntegerVector OptEmbedDim(Rcpp::NumericMatrix Emat);

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
double OptThetaParm(Rcpp::NumericMatrix Thetamat);

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
Rcpp::IntegerVector OptICparm(Rcpp::NumericMatrix Emat);

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
Rcpp::NumericMatrix MatNotNAIndice(Rcpp::NumericMatrix mat, bool byrow = true);

#endif // HelperFuns
