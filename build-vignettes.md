# NA

Due to the time-consuming computations involved in the vignette of the
*tEDM* package, it is necessary to pre-build the vignette prior to
package submission.

``` r
.prebuild_vignettes = \(name){
  out = paste0("vignettes/",name,".Rmd")
  inp = paste0(out,".orig")
  knitr::knit(inp,out)
}

.prebuild_vignettes("tEDM")
```
