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
  tau = 1,
  k = E + 1,
  lib = NULL,
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

- lib:

  (optional) libraries indices.

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
obs = runif(10,0,0.1)
sim = vector("list",5)
for (i in seq_along(obs)){
  sim[[i]] = logistic_map(x = obs[i],y = obs[i],step = 15,beta_xy = 0.5,beta_yx = 0)
}
lst = list(x = do.call(cbind, lapply(sim, function(df) df$x)),
           y = do.call(cbind, lapply(sim, function(df) df$y)))
multispatialccm(lst,"x","y",libsizes = 1:5,E = 2,k = 3,threads = 1)
#> Computing: [========================================] 100% (done)                         
#> Computing: [========================================] 100% (done)                         
#>   libsizes      x->y      y->x
#> 1        1 0.3131501 0.3161107
#> 2        2 0.6804184 0.3629746
#> 3        3 0.7935703 0.4772453
#> 4        4 0.8390027 0.5316734
#> 5        5 0.8691234 0.5819187
```
