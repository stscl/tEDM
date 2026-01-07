.pcm_ts_method = \(data, cause, effect, conds, libsizes = NULL, E = 3, tau = 0, k = E+1, theta = 1, algorithm = "simplex", lib = NULL, pred = NULL, dist.metric = "L1",
                   dist.average = TRUE, threads = length(pred), parallel.level = "low", bidirectional = TRUE, cumulate = FALSE, progressbar = TRUE){
  varname = .check_character(c(cause, effect, conds))
  E = .check_inputelementnum(E,length(varname),length(conds))
  tau = .check_inputelementnum(tau,length(varname))
  k = .check_inputelementnum(k,length(varname))
  pl = .check_parallellevel(parallel.level)
  cause = data[,cause,drop = TRUE]
  effect = data[,effect,drop = TRUE]
  condmat = as.matrix(data[,conds,drop = FALSE])

  if (is.null(lib)) lib = .internal_library(data)
  if (is.null(pred)) pred = lib
  if (is.null(libsizes)) libsizes = length(lib)
  if (threads == 0) threads = length(pred)

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppPCM(cause,effect,condmat,libsizes,lib,pred,E[-2],tau[-2],k[-2],simplex,theta,threads,pl,cumulate,.check_distmetric(dist.metric),dist.average,progressbar)
  }
  y_xmap_x = RcppPCM(effect,cause,condmat,libsizes,lib,pred,E[-1],tau[-1],k[-1],simplex,theta,threads,pl,cumulate,.check_distmetric(dist.metric),dist.average,progressbar)

  return(.bind_xmapdf2(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' partial cross mapping
#'
#' @inheritParams ccm
#' @param conds name of conditioning variables.
#' @param cumulate (optional) serial or cumulative computation of partial cross mapping.
#'
#' @return A list
#' \describe{
#' \item{\code{pxmap}}{partial cross mapping results}
#' \item{\code{xmap}}{cross mapping results}
#' \item{\code{varname}}{names of causal, effect and conditioning variables}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name pcm
#' @aliases pcm,data.frame-method
#' @references
#' Leng, S., Ma, H., Kurths, J. et al. Partial cross mapping eliminates indirect causal influences. Nat Commun 11, 2632 (2020).
#'
#' @examples
#' sim = logistic_map(x = 0.4,y = 0.4,z = 0.4,step = 45,
#'                    beta_xy = 0.5, beta_xz = 0,
#'                    beta_yx = 0, beta_yz = 0.5,
#'                    beta_zx = 0, beta_zy = 0)
#' pcm(sim,"x","z","y",libsizes = seq(5,45,5),E = 10,k = 7,threads = 1)
#'
methods::setMethod("pcm", "data.frame", .pcm_ts_method)
