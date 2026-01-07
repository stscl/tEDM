# cross mapping cardinality

cross mapping cardinality

## Usage

``` r
# S4 method for class 'data.frame'
cmc(
  data,
  cause,
  effect,
  libsizes = NULL,
  E = 3,
  tau = 1,
  k = pmin(E^2),
  lib = NULL,
  pred = NULL,
  dist.metric = "L1",
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

- lib:

  (optional) libraries indices.

- pred:

  (optional) predictions indices.

- dist.metric:

  (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).

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

- `cs`:

  causal strength

- `varname`:

  names of causal and effect variables

- `bidirectional`:

  whether to examine bidirectional causality

## References

Tao, P., Wang, Q., Shi, J., Hao, X., Liu, X., Min, B., Zhang, Y., Li,
C., Cui, H., Chen, L., 2023. Detecting dynamical causality by
intersection cardinal concavity. Fundamental Research.

## Examples

``` r
sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
cmc(sim,"x","y",E = 4,k = 15,threads = 1)
#> Computing: [========================================] 100% (done)                         
#> Computing: [========================================] 100% (done)                         
#>   neighbors      x->y      y->x
#> 1        15 0.2666667 0.1333333
```
