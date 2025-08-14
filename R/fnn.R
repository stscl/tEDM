.fnn_ts_method = \(data, target, lib = NULL, pred = NULL, E = 2:10, tau = 1, 
                   dist.metric = "L1", rt = 10, eps = 2, threads = length(E)){
  vec = .uni_ts(data,target)
  rt = .check_inputelementnum(rt,max(E))
  eps = .check_inputelementnum(eps,max(E))
  if (is.null(lib)) lib = which(!is.na(vec))
  if (is.null(pred)) pred = lib
  return(RcppFNN4TS(vec,rt,eps,lib,pred,E,tau,.check_distmetric(dist.metric),threads))
}

#' false nearest neighbours
#'
#' @inheritParams embedded
#' @param lib (optional) libraries indices.
#' @param pred (optional) predictions indices.
#' @param dist.metric (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).
#' @param rt (optional) escape factor.
#' @param eps (optional) neighborhood diameter.
#' @param threads (optional) number of threads to use.
#'
#' @return A vector
#' @export
#' @name fnn
#' @aliases fnn,data.frame-method
#' @references
#' Kennel M. B., Brown R. and Abarbanel H. D. I., Determining embedding dimension for phase-space reconstruction using a geometrical construction, Phys. Rev. A, Volume 45, 3403 (1992).
#'
#' @examples
#' sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
#' fnn(sim,"x",threads = 1)
#'
methods::setMethod("fnn", "data.frame", .fnn_ts_method)
