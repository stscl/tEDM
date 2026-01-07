# false nearest neighbours

false nearest neighbours

## Usage

``` r
# S4 method for class 'data.frame'
fnn(
  data,
  target,
  lib = NULL,
  pred = NULL,
  E = 2:10,
  tau = 1,
  dist.metric = "L1",
  rt = 10,
  eps = 2,
  threads = length(E)
)
```

## Arguments

- data:

  observation data.

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

- dist.metric:

  (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).

- rt:

  (optional) escape factor.

- eps:

  (optional) neighborhood diameter.

- threads:

  (optional) number of threads to use.

## Value

A vector

## References

Kennel M. B., Brown R. and Abarbanel H. D. I., Determining embedding
dimension for phase-space reconstruction using a geometrical
construction, Phys. Rev. A, Volume 45, 3403 (1992).

## Examples

``` r
sim = logistic_map(x = 0.4,y = 0.4,step = 45,beta_xy = 0.5,beta_yx = 0)
fnn(sim,"x",threads = 1)
#>       E:1       E:2       E:3       E:4       E:5       E:6       E:7       E:8 
#> 0.1944444 0.1388889 0.1111111 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 
#>       E:9 
#> 0.0000000 
```
