# convergent cross mapping

convergent cross mapping

## Usage

``` r
# S4 method for class 'data.frame'
ccm(
  data,
  cause,
  effect,
  libsizes = NULL,
  E = 3,
  tau = 0,
  k = E + 1,
  theta = 1,
  algorithm = "simplex",
  lib = NULL,
  pred = NULL,
  dist.metric = "L1",
  dist.average = TRUE,
  threads = length(pred),
  parallel.level = "low",
  bidirectional = TRUE,
  progressbar = TRUE
)
```

## Arguments

- data:

  observation data.

- cause:

  name of causal variable.

- effect:

  name of effect variable.

- libsizes:

  (optional) number of time points used.

- E:

  (optional) embedding dimensions.

- tau:

  (optional) step of time lags.

- k:

  (optional) number of nearest neighbors.

- theta:

  (optional) weighting parameter for distances, useful when `algorithm`
  is `smap`.

- algorithm:

  (optional) prediction algorithm.

- lib:

  (optional) libraries indices.

- pred:

  (optional) predictions indices.

- dist.metric:

  (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).

- dist.average:

  (optional) whether to average distance.

- threads:

  (optional) number of threads to use.

- parallel.level:

  (optional) level of parallelism, `low` or `high`.

- bidirectional:

  (optional) whether to examine bidirectional causality.

- progressbar:

  (optional) whether to show the progress bar.

## Value

A list

- `xmap`:

  cross mapping results

- `varname`:

  names of causal and effect variable

- `bidirectional`:

  whether to examine bidirectional causality

## References

Sugihara, G., May, R., Ye, H., Hsieh, C., Deyle, E., Fogarty, M., Munch,
S., 2012. Detecting Causality in Complex Ecosystems. Science 338,
496–500.

## Examples

``` r
sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
ccm(sim,"x","y",libsizes = seq(5,45,5),E = 10,k = 7,threads = 1)
#> Computing: [========================================] 100% (done)                         
#> Computing: [========================================] 100% (done)                         
#>   libsizes      x->y      y->x
#> 1        5 0.7573109 0.5602115
#> 2       10 0.8528876 0.6196130
#> 3       15 0.8662672 0.6233221
#> 4       20 0.8631069 0.6221401
#> 5       25 0.8604142 0.6195158
#> 6       30 0.8574624 0.6180640
#> 7       35 0.8564398 0.6189599
#> 8       36 0.9069745 0.7032761
```
