# simplex forecast

simplex forecast

## Usage

``` r
# S4 method for class 'data.frame'
simplex(
  data,
  column,
  target,
  lib = NULL,
  pred = NULL,
  E = 2:10,
  tau = 1,
  k = E + 1,
  dist.metric = "L1",
  dist.average = TRUE,
  threads = length(E)
)

# S4 method for class 'list'
simplex(
  data,
  column,
  target,
  lib = NULL,
  pred = NULL,
  E = 2:10,
  tau = 1,
  k = E + 1,
  dist.metric = "L1",
  dist.average = TRUE,
  threads = length(E)
)
```

## Arguments

- data:

  observation data.

- column:

  name of library variable.

- target:

  name of target variable.

- lib:

  (optional) libraries indices.

- pred:

  (optional) predictions indices.

- E:

  (optional) embedding dimensions.

- tau:

  (optional) step of time lags.

- k:

  (optional) number of nearest neighbors used in prediction.

- dist.metric:

  (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).

- dist.average:

  (optional) whether to average distance.

- threads:

  (optional) number of threads to use.

## Value

A list

- `xmap`:

  forecast performance

- `varname`:

  name of target variable

- `method`:

  method of cross mapping

## References

Sugihara G. and May R. 1990. Nonlinear forecasting as a way of
distinguishing chaos from measurement error in time series. Nature,
344:734-741.

## Examples

``` r
sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
simplex(sim,"x","y",E = 4:10,k = 7,threads = 1)
#> The suggested E,k,tau for variable y is 10, 7 and 1 
```
