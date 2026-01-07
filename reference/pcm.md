# partial cross mapping

partial cross mapping

## Usage

``` r
# S4 method for class 'data.frame'
pcm(
  data,
  cause,
  effect,
  conds,
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
  cumulate = FALSE,
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

- conds:

  name of conditioning variables.

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

- cumulate:

  (optional) serial or cumulative computation of partial cross mapping.

- progressbar:

  (optional) whether to show the progress bar.

## Value

A list

- `pxmap`:

  partial cross mapping results

- `xmap`:

  cross mapping results

- `varname`:

  names of causal and effect variable

- `bidirectional`:

  whether to examine bidirectional causality

## References

Leng, S., Ma, H., Kurths, J. et al. Partial cross mapping eliminates
indirect causal influences. Nat Commun 11, 2632 (2020).

## Examples

``` r
sim = logistic_map(x = 0.4,y = 0.4,z = 0.4,step = 45,
                   beta_xy = 0.5, beta_xz = 0,
                   beta_yx = 0, beta_yz = 0.5,
                   beta_zx = 0, beta_zy = 0)
pcm(sim,"x","z","y",libsizes = seq(5,45,5),E = 10,k = 7,threads = 1)
#> Computing: [========================================] 100% (done)                         
#> Computing: [========================================] 100% (done)                         
#> -------------------------------------- 
#> ***partial cross mapping prediction*** 
#> -------------------------------------- 
#>   libsizes        x->z          z->x
#> 1        5  0.35911303  1.260417e-01
#> 2       10  0.19561345 -4.068307e-05
#> 3       15  0.08032524  4.125853e-03
#> 4       20  0.01890635  5.573376e-02
#> 5       25 -0.01338321  1.036854e-01
#> 6       27  0.05144222  2.845048e-01
#> 
#> ------------------------------ 
#> ***cross mapping prediction*** 
#> ------------------------------ 
#>   libsizes      x->z      z->x
#> 1        5 0.7714550 0.5182871
#> 2       10 0.8545976 0.5570433
#> 3       15 0.8718242 0.5648568
#> 4       20 0.8733124 0.5715475
#> 5       25 0.8732259 0.5722011
#> 6       27 0.9304557 0.5953293
```
