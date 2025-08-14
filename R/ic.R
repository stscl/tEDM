.ic_ts_method = \(data, column, target, lib = NULL, pred = NULL, E = 2:10, tau = 1, k = E+2,
                  dist.metric = "L1", threads = length(pred), parallel.level = "low"){
  vx = .uni_ts(data,column)
  vy = .uni_ts(data,target)
  pl = .check_parallellevel(parallel.level)

  if (is.null(lib)) lib = .internal_library(cbind(vx,vy))
  if (is.null(pred)) pred = lib
  if (threads == 0) threads = length(pred)

  res = RcppIC4TS(vx,vy,lib,pred,E,k,tau,0,.check_distmetric(dist.metric),threads,pl)
  return(.bind_xmapself(res,target,"ic",tau))
}

#' intersection cardinality
#'
#' @inheritParams simplex
#' @param parallel.level (optional) level of parallelism, `low` or `high`.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping performance}
#' \item{\code{varname}}{name of target variable}
#' \item{\code{method}}{method of cross mapping}
#' \item{\code{tau}}{step of time lag}
#' }
#' @export
#' @name ic
#' @aliases ic,data.frame-method
#' @references
#' Tao, P., Wang, Q., Shi, J., Hao, X., Liu, X., Min, B., Zhang, Y., Li, C., Cui, H., Chen, L., 2023. Detecting dynamical causality by intersection cardinal concavity. Fundamental Research.
#'
#' @examples
#' sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
#' ic(sim,"x","y",E = 4,k = 15:30,threads = 1)
#'
methods::setMethod("ic", "data.frame", .ic_ts_method)
