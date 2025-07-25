% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logostic_map.R
\name{logistic_map}
\alias{logistic_map}
\title{logistic map}
\usage{
logistic_map(
  x,
  y = NULL,
  z = NULL,
  step = 15,
  alpha_x = 3.6,
  alpha_y = 3.72,
  alpha_z = 3.68,
  beta_xy = 0.05,
  beta_xz = 0.05,
  beta_yx = 0.2,
  beta_yz = 0.2,
  beta_zx = 0.35,
  beta_zy = 0.35,
  threshold = Inf,
  transient = 1
)
}
\arguments{
\item{x}{value x.}

\item{y}{(optional) value y.}

\item{z}{(optional) value z.}

\item{step}{(optional) number of simulation time steps.}

\item{alpha_x}{(optional) growth parameter for x.}

\item{alpha_y}{(optional) growth parameter for y.}

\item{alpha_z}{(optional) growth parameter for y.}

\item{beta_xy}{(optional) cross-inhibition from x to y.}

\item{beta_xz}{(optional) cross-inhibition from x to z.}

\item{beta_yx}{(optional) cross-inhibition from y to x.}

\item{beta_yz}{(optional) cross-inhibition from y to z.}

\item{beta_zx}{(optional) cross-inhibition from z to x.}

\item{beta_zy}{(optional) cross-inhibition from z to y.}

\item{threshold}{(optional) set to \code{NaN} if the absolute value exceeds this threshold.}

\item{transient}{(optional) transients to be excluded from the results.}
}
\value{
A data.frame
}
\description{
logistic map
}
\examples{
logistic_map(x = 0.2)

}
