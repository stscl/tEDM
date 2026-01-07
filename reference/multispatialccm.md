# multispatial convergent cross mapping

multispatial convergent cross mapping

## Usage

``` r
# S4 method for class 'list'
multispatialccm(
  data,
  cause,
  effect,
  libsizes,
  E = 3,
  tau = 0,
  k = E + 1,
  boot = 99,
  seed = 42,
  dist.metric = "L1",
  dist.average = TRUE,
  threads = length(libsizes),
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

  number of time points used in prediction.

- E:

  (optional) embedding dimensions.

- tau:

  (optional) step of time lags.

- k:

  (optional) number of nearest neighbors used in prediction.

- boot:

  (optional) number of bootstraps to perform.

- seed:

  (optional) random seed.

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

  names of causal and effect variables

- `bidirectional`:

  whether to examine bidirectional causality

## References

Clark, A.T., Ye, H., Isbell, F., Deyle, E.R., Cowles, J., Tilman, G.D.,
Sugihara, G., 2015. Spatial convergent cross mapping to detect causal
relationships from short time series. Ecology 96, 1174–1181.

## Examples

``` r
set.seed(42)
obs = runif(15,0,0.1)
sim = vector("list",15)
for (i in seq_along(obs)){
  sim[[i]] = logistic_map(x = obs[i],y = obs[i],step = 15,beta_xy = 0.5,beta_yx = 0)
}
lst = list(x = do.call(cbind, lapply(sim, function(df) df$x)),
           y = do.call(cbind, lapply(sim, function(df) df$y)))
multispatialccm(lst,"x","y",libsizes = 5:15,E = 2,k = 3,threads = 1)
#> Computing: [========================================] 100% (done)                         
#> Computing: [========================================] 100% (done)                         
#>    libsizes      x->y      y->x
#> 1         5 0.7667238 0.5809166
#> 2         6 0.8105505 0.6238409
#> 3         7 0.8353561 0.6530469
#> 4         8 0.8518396 0.6786115
#> 5         9 0.8641090 0.6994746
#> 6        10 0.8754431 0.7134046
#> 7        11 0.8893507 0.7334001
#> 8        12 0.9021861 0.7570492
#> 9        13 0.9118122 0.7803526
#> 10       14 0.9218980 0.7985672
#> 11       15 0.9300771 0.8139998
```
