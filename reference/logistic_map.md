# logistic map

logistic map

## Usage

``` r
logistic_map(
  x,
  y = NULL,
  z = NULL,
  step = 15,
  alpha_x = 3.6,
  alpha_y = 3.72,
  alpha_z = 3.68,
  beta_xy = 0.05,
  beta_xz = 0.05,
  beta_yx = 0.2,
  beta_yz = 0.2,
  beta_zx = 0.35,
  beta_zy = 0.35,
  threshold = Inf,
  transient = 1
)
```

## Arguments

- x:

  value x.

- y:

  (optional) value y.

- z:

  (optional) value z.

- step:

  (optional) number of simulation time steps.

- alpha_x:

  (optional) growth parameter for x.

- alpha_y:

  (optional) growth parameter for y.

- alpha_z:

  (optional) growth parameter for y.

- beta_xy:

  (optional) cross-inhibition from x to y.

- beta_xz:

  (optional) cross-inhibition from x to z.

- beta_yx:

  (optional) cross-inhibition from y to x.

- beta_yz:

  (optional) cross-inhibition from y to z.

- beta_zx:

  (optional) cross-inhibition from z to x.

- beta_zy:

  (optional) cross-inhibition from z to y.

- threshold:

  (optional) set to `NaN` if the absolute value exceeds this threshold.

- transient:

  (optional) transients to be excluded from the results.

## Value

A data.frame

## Examples

``` r
logistic_map(x = 0.2)
#>            x
#> 1  0.5760000
#> 2  0.8792064
#> 3  0.3823290
#> 4  0.8501527
#> 5  0.4586150
#> 6  0.8938342
#> 7  0.3416206
#> 8  0.8096975
#> 9  0.5547149
#> 10 0.8892226
#> 11 0.3546207
#> 12 0.8239135
#> 13 0.5222881
#> 14 0.8982117
#> 15 0.3291389
```
