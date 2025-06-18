.simplex_ts_method = \(data,column,target,lib = NULL,pred = NULL,
                       E = 1:10, tau = 0, k = E+1, threads = length(E)){
  vx = .uni_ts(data,column)
  vy = .uni_ts(data,target)
  if (is.null(lib)) lib = .internal_library(cbind(vx,vy))
  if (is.null(pred)) pred = lib
  res = RcppSimplex4TS(vx,vy,lib,pred,E,k,tau,threads)
  return(.bind_xmapself(res,target))
}

#' simplex forecast
#'
#' @inheritParams embedded
#' @param column name of library variable.
#' @param lib (optional) libraries indices.
#' @param pred (optional) predictions indices.
#' @param k (optional) number of nearest neighbors used in prediction.
#' @param threads (optional) number of threads to use.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{self mapping prediction results}
#' \item{\code{varname}}{name of target variable}
#' }
#' @export
#' @name simplex
#' @aliases simplex,data.frame-method
#' @references
#' Sugihara G. and May R. 1990. Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series. Nature, 344:734-741.
#'
#' @examples
#' sim = as.data.frame(logistic_map(x = 0.2, y = 0.4, step = 45, beta_xy = 0.5, beta_yx = 0))
#' simplex(sim,"x","y")
#' ccm(sim,"x","y",libsizes = seq(5,45,5))
#'
methods::setMethod("simplex", "data.frame", .simplex_ts_method)
