.ccm_ts_method = \(data, cause, effect, libsizes, E = 3, tau = 1, k = E+2, theta = 1, algorithm = "simplex", lib = NULL,
                   pred = NULL, threads = 1, parallel.level = "low", bidirectional = TRUE, progressbar = TRUE){
  E = .check_inputelementnum(E,2)
  tau = .check_inputelementnum(tau,2)
  k = .check_inputelementnum(k,2)
  pl = .check_parallellevel(parallel.level)
  if (is.null(lib)) lib = .internal_library(data)
  if (is.null(pred)) pred = lib
  cause = .uni_ts(data,cause)
  effect = .uni_ts(data,effect)

  simplex = ifelse(algorithm == "simplex", TRUE, FALSE)
  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppCCM(cause,effect,libsizes,lib,pred,E[1],tau[1],k[1],simplex,theta,threads,pl,progressbar)
  }
  y_xmap_x = RcppCCM(effect,cause,libsizes,lib,pred,E[2],tau[2],k[2],simplex,theta,threads,pl,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' convergent cross mapping
#'
#' @param data observation data.
#' @param cause name of causal variable.
#' @param effect name of effect variable.
#' @param libsizes number of time points used in prediction.
#' @param E (optional) embedding dimensions.
#' @param tau (optional) step of time lags.
#' @param k (optional) number of nearest neighbors used in prediction.
#' @param theta (optional) weighting parameter for distances, useful when `algorithm` is `smap`.
#' @param algorithm (optional) prediction algorithm.
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
#' \item{\code{varname}}{names of causal and effect variable}
#' \item{\code{bidirectional}}{whether to examine bidirectional causality}
#' }
#' @export
#' @name ccm
#' @aliases ccm,ts-method
#' @references
#' Sugihara, G., May, R., Ye, H., Hsieh, C., Deyle, E., Fogarty, M., Munch, S., 2012. Detecting Causality in Complex Ecosystems. Science 338, 496–500.
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' g = gccm(columbus,"hoval","crime",libsizes = seq(5,45,5),E = 6)
#' g
#' plot(g, ylimits = c(0,0.85))
#' }

# methods::setMethod("ccm", "data.frame", .ccm_ts_method)
