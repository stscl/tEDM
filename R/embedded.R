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
#' @aliases embedded,data.frame-method
#'
#' @examples
#' embedded(data.frame(t = 1:5),"t",3)
#'
methods::setMethod("embedded", "data.frame", .embedded_ts_method)
