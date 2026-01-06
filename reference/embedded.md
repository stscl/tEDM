# embedding time series data

embedding time series data

## Usage

``` r
# S4 method for class 'data.frame'
embedded(data, target, E = 3, tau = 1)
```

## Arguments

- data:

  observation data.

- target:

  name of target variable.

- E:

  (optional) embedding dimensions.

- tau:

  (optional) step of time lags.

## Value

A matrix

## Examples

``` r
embedded(data.frame(t = 1:5),"t",3)
#>      [,1] [,2] [,3]
#> [1,]    1  NaN  NaN
#> [2,]    2    1  NaN
#> [3,]    3    2    1
#> [4,]    4    3    2
#> [5,]    5    4    3
```
