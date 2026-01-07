.ccm_ts_method = \(data, cause, effect, libsizes = NULL, E = 3, tau = 0, k = E+1, theta = 1, algorithm = "simplex", lib = NULL, pred = NULL,
                   dist.metric = "L1",dist.average = TRUE,threads = length(pred),parallel.level = "low",bidirectional = TRUE,progressbar = TRUE){
  varname = .check_character(cause,effect)
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

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppCCM(cause,effect,libsizes,lib,pred,E[1],tau[1],k[1],simplex,theta,threads,pl,.check_distmetric(dist.metric),dist.average,progressbar)
  }
  y_xmap_x = RcppCCM(effect,cause,libsizes,lib,pred,E[2],tau[2],k[2],simplex,theta,threads,pl,.check_distmetric(dist.metric),dist.average,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' convergent cross mapping
#'
#' @inheritParams cmc
#' @param theta (optional) weighting parameter for distances, useful when `algorithm` is `smap`.
#' @param algorithm (optional) prediction algorithm.
#' @param dist.average (optional) whether to average distance.
#'
#' @return A list
#' \describe{
#' \item{\code{xmap}}{cross mapping results}
#' \item{\code{varname}}{names of causal and effect variables}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name ccm
#' @aliases ccm,data.frame-method
#' @references
#' Sugihara, G., May, R., Ye, H., Hsieh, C., Deyle, E., Fogarty, M., Munch, S., 2012. Detecting Causality in Complex Ecosystems. Science 338, 496–500.
#'
#' @examples
#' sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
#' ccm(sim,"x","y",libsizes = seq(5,45,5),E = 10,k = 7,threads = 1)
#'
methods::setMethod("ccm", "data.frame", .ccm_ts_method)
