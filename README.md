# PCEV
[![Build Status](https://travis-ci.org/turgeonmaxime/multiKernel.svg?branch=master)](https://travis-ci.org/turgeonmaxime/multiKernel)

This ```R``` package implements multivariate prediction using kernel-machine regression.

## Installation

This package can be directly installed from GitHub using the [devtools](http://cran.r-project.org/web/packages/devtools/index.html) package:

``` r
library(devtools)
devtools::install_github('turgeonmaxime/multiKernel')
```

However, you should be aware that, even though the fitting functions are most probably correct, my implementation of cross-validation is most probably wrong.
