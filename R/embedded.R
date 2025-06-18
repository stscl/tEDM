.embedded_ts_method = \(data,target,E = 3,tau = 0){
  vec = .uni_ts(data,target)
  return(RcppEmbed(vec,E,tau))
}

#' embedding time series data
#'
#' @param data observation data.
#' @param target name of target variable.
#' @param E (optional) embedding dimensions.
#' @param tau (optional) step of time lags.
#'
#' @return A matrix
#' @export
#' @name embedded
#' @aliases embedded,ts-method
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' v = embedded(columbus,"crime")
#' v[1:5,]
#'
#' cu = terra::rast(system.file("case/cu.tif", package="spEDM"))
#' r = embedded(cu,"cu")
#' r[1:5,]
#'
methods::setMethod("embedded", "data.frame", .embedded_ts_method)
