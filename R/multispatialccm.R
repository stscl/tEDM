.multispatialccm_ts_method = \(data, cause, effect, libsizes, E = 3, tau = 0, k = E+1, boot = 99, seed = 42,
                               threads = length(libsizes), parallel.level = "low", bidirectional = TRUE, progressbar = TRUE){
  varname = .check_character(cause,effect)
  E = .check_inputelementnum(E,2)
  tau = .check_inputelementnum(tau,2)
  k = .check_inputelementnum(k,2)
  pl = .check_parallellevel(parallel.level)
  cause = as.matrix(data[[cause]])
  effect = as.matrix(data[[effect]])

  x_xmap_y = NULL
  if (bidirectional){
    x_xmap_y = RcppMultispatialCCM(cause,effect,libsizes,E[1],tau[1],k[1],boot,seed,threads,pl,progressbar)
  }
  y_xmap_x = RcppMultispatialCCM(effect,cause,libsizes,E[2],tau[2],k[2],boot,seed,threads,pl,progressbar)

  return(.bind_xmapdf(varname,x_xmap_y,y_xmap_x,bidirectional))
}

#' multispatial convergent cross mapping
#'
#' @param data observation data.
#' @param cause name of causal variable.
#' @param effect name of effect variable.
#' @param libsizes number of time points used in prediction.
#' @param E (optional) embedding dimensions.
#' @param tau (optional) step of time lags.
#' @param k (optional) number of nearest neighbors used in prediction.
#' @param boot (optional) number of bootstraps to perform.
#' @param seed (optional) random seed.
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
#' @name multispatialccm
#' @aliases multispatialccm,list-method
#' @references
#' Clark, A.T., Ye, H., Isbell, F., Deyle, E.R., Cowles, J., Tilman, G.D., Sugihara, G., 2015. Spatial convergent cross mapping to detect causal relationships from short time series. Ecology 96, 1174–1181.
#'
#' @examples
#' set.seed(42)
#' obs = runif(15,0,0.1)
#' sim = vector("list",15)
#' for (i in seq_along(obs)){
#'   sim[[i]] = logistic_map(x = obs[i],y = obs[i],step = 15,beta_xy = 0.5,beta_yx = 0)
#' }
#'
#' lst = list(x = do.call(rbind, lapply(sim, function(df) df$x)),
#'            y = do.call(rbind, lapply(sim, function(df) df$y)))
#' multispatialccm(lst,"x","y",libsizes = seq(5,15,1),E = c(2,4),k = 5,threads = 1)
#'
methods::setMethod("multispatialccm", "list", .multispatialccm_ts_method)
