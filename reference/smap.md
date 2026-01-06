# smap forecast

smap forecast

## Usage

``` r
# S4 method for class 'data.frame'
smap(
  data,
  column,
  target,
  lib = NULL,
  pred = NULL,
  E = 3,
  tau = 1,
  k = E + 1,
  dist.metric = "L1",
  dist.average = TRUE,
  theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3,
    4, 6, 8),
  threads = length(theta)
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

- theta:

  (optional) weighting parameter for distances.

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

Sugihara G. 1994. Nonlinear forecasting for the classification of
natural time series. Philosophical Transactions: Physical Sciences and
Engineering, 348 (1688):477-495.

## Examples

``` r
sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
smap(sim,"x","y",E = 10,k = 7,threads = 1)
#> The suggested theta for variable y is 8 
```
