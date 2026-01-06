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
#>   libsizes        x->z       z->x
#> 1        5 0.410387689 0.11654973
#> 2       10 0.316685463 0.04464733
#> 3       15 0.256396689 0.02586971
#> 4       20 0.250855810 0.07990501
#> 5       25 0.250303599 0.09410492
#> 6       27 0.007456209 0.15394446
#> 
#> ------------------------------ 
#> ***cross mapping prediction*** 
#> ------------------------------ 
#>   libsizes      x->z      z->x
#> 1        5 0.8104853 0.6090928
#> 2       10 0.8695076 0.6226428
#> 3       15 0.8802282 0.6245086
#> 4       20 0.8808134 0.6285323
#> 5       25 0.8801786 0.6288845
#> 6       27 0.9314565 0.6174386
```
