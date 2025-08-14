.cmc_ts_method = \(data, cause, effect, libsizes = NULL, E = 3, tau = 1, k = pmin(E^2), lib = NULL, pred = NULL,
                   dist.metric = "L1", threads = length(pred), parallel.level = "low", bidirectional = TRUE, progressbar = TRUE){
  varname = .check_character(cause, effect)
  E = .check_inputelementnum(E,2)
  tau = .check_inputelementnum(tau,2)
  k = .check_inputelementnum(k,2)
  pl = .check_parallellevel(parallel.level)

  cause = .uni_ts(data,cause)
  effect = .uni_ts(data,effect)

  if (is.null(lib)) lib = .internal_library(data)
  if (is.null(pred)) pred = lib
  if (is.null(libsizes)) libsizes = length(lib)
  if (threads == 0) threads = length(pred)

  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppCMC(cause,effect,libsizes,lib,pred,E,tau,k[1],0,.check_distmetric(dist.metric),threads,pl,progressbar)
  }
  y_xmap_x = RcppCMC(effect,cause,libsizes,lib,pred,rev(E),rev(tau),k[2],0,.check_distmetric(dist.metric),threads,pl,progressbar)

  return(.bind_intersectdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' cross mapping cardinality
#'
#' @param data observation data.
#' @param cause name of causal variable.
#' @param effect name of effect variable.
#' @param libsizes (optional) number of time points used.
#' @param E (optional) embedding dimensions.
#' @param tau (optional) step of time lags.
#' @param k (optional) number of nearest neighbors.
#' @param dist.metric (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).
#' @param lib (optional) libraries indices.
#' @param pred (optional) predictions indices.
#' @param threads (optional) number of threads to use.
#' @param parallel.level (optional) level of parallelism, `low` or `high`.
#' @param bidirectional (optional) whether to examine bidirectional causality.
#' @param progressbar (optional) whether to show the progress bar.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping results}
#' \item{\code{cs}}{causal strength}
#' \item{\code{varname}}{names of causal and effect variable}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name cmc
#' @aliases cmc,data.frame-method
#' @references
#' Tao, P., Wang, Q., Shi, J., Hao, X., Liu, X., Min, B., Zhang, Y., Li, C., Cui, H., Chen, L., 2023. Detecting dynamical causality by intersection cardinal concavity. Fundamental Research.
#'
#' @examples
#' sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
#' cmc(sim,"x","y",E = 4,k = 15,threads = 1)
#'
methods::setMethod("cmc", "data.frame", .cmc_ts_method)
