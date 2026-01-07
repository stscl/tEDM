# intersection cardinality

intersection cardinality

## Usage

``` r
# S4 method for class 'data.frame'
ic(
  data,
  column,
  target,
  lib = NULL,
  pred = NULL,
  E = 2:10,
  tau = 1,
  k = E + 2,
  dist.metric = "L1",
  threads = length(pred),
  parallel.level = "low"
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

- threads:

  (optional) number of threads to use.

- parallel.level:

  (optional) level of parallelism, `low` or `high`.

## Value

A list

- `xmap`:

  cross mapping performance

- `varname`:

  name of target variable

- `method`:

  method of cross mapping

## References

Tao, P., Wang, Q., Shi, J., Hao, X., Liu, X., Min, B., Zhang, Y., Li,
C., Cui, H., Chen, L., 2023. Detecting dynamical causality by
intersection cardinal concavity. Fundamental Research.

## Examples

``` r
sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
ic(sim,"x","y",E = 4,k = 15:30,threads = 1)
#> The suggested E,k,tau for variable y is 4, 15 and 1 
```
