
# n2oflux

<!-- badges: start -->
<!-- badges: end -->

The goal of n2oflux is to help users of the Licor 8200 Smart Chamber to create a table based on the json files for a more flexible analysis of the data collected. The package also allows to calculate fluxes using a linear and non-linear regression following the same methodology used by the smart chamber. The deadband can also be optimized within a certain range.

## Installation

You can install the development version of n2oflux from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("antodi/n2oflux")
```

or 

``` r
# install.packages("devtools")
install_github("antodi/n2oflux")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(n2oflux)
process_json_files(PATH)
calculate_n2o_flux(data,deadband=30,deadband_c=0,stop_time_ag=120,offset_k="json",opt_db="no",use_json_parameters=0)
```

