#' logistic map
#'
#' @param x value x.
#' @param y (optional) value y.
#' @param z (optional) value z.
#' @param step (optional) number of simulation time steps.
#' @param alpha_x (optional) growth parameter for x.
#' @param alpha_y (optional) growth parameter for y.
#' @param alpha_z (optional) growth parameter for y.
#' @param beta_xy (optional) cross-inhibition from x to y.
#' @param beta_xz (optional) cross-inhibition from x to z.
#' @param beta_yx (optional) cross-inhibition from y to x.
#' @param beta_yz (optional) cross-inhibition from y to z.
#' @param beta_zx (optional) cross-inhibition from z to x.
#' @param beta_zy (optional) cross-inhibition from z to y.
#' @param threshold (optional) set to `NaN` if the absolute value exceeds this threshold.
#' @param transient (optional) transients to be excluded from the results.
#'
#' @return A data.frame
#' @export
#'
#' @examples
#' logistic_map(x = 0.2)
#'
logistic_map = \(x, y = NULL, z = NULL, step = 15, alpha_x = 3.6, alpha_y = 3.72, alpha_z = 3.68,
                 beta_xy = 0.05, beta_xz = 0.05, beta_yx = 0.2, beta_yz = 0.2, beta_zx = 0.35, beta_zy = 0.35,
                 threshold = Inf, transient = 1){
  xl = x; yl = y; zl = z;
  if (is.null(x)) xl = 0;
  if (is.null(y)) yl = 0;
  if (is.null(z)) zl = 0;
  res = lapply(RcppLogisticMap(xl,yl,zl,step,alpha_x,alpha_y,alpha_z,beta_xy,beta_xz,beta_yx,beta_yz,beta_zx,beta_zy,threshold),
               \(.x) .x[-unique(abs(transient))])
  if (is.null(y)) res$y = NULL
  if (is.null(z)) res$z = NULL
  return(as.data.frame(res))
}
