.smap_ts_method = \(data,column,target,lib = NULL,pred = NULL,E = 3,tau = 0,k = E+1,
                    theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03,
                              0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
                    threads = length(theta)){
  vx = .uni_ts(data,column)
  vy = .uni_ts(data,target)
  if (is.null(lib)) lib = .internal_library(cbind(vx,vy))
  if (is.null(pred)) pred = lib
  res = RcppSMap4TS(vx,vy,lib,pred,theta,E,tau,k,threads)
  return(.bind_xmapself(res,target,"smap"))
}

#' smap forecast
#'
#' @inheritParams simplex
#' @param theta (optional) weighting parameter for distances.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{forecast performance}
#' \item{\code{varname}}{name of target variable}
#' \item{\code{method}}{method of cross mapping}
#' }
#' @export
#' @name smap
#' @aliases smap,data.frame-method
#' @references
#' Sugihara G. 1994. Nonlinear forecasting for the classification of natural time series. Philosophical Transactions: Physical Sciences and Engineering, 348 (1688):477-495.
#'
#' @examples
#' sim = as.data.frame(logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0))
#' smap(sim,"x","y",E = 7,k = 7,threads = 1)
#'
methods::setMethod("smap", "data.frame", .smap_ts_method)
