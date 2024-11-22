
# n2oflux

<!-- badges: start -->
<!-- badges: end -->

The goal of n2oflux is to help users of the Licor 8200 Smart Chamber to create a table based on the json files for a more flexible analysis of the data collected. The package also allows to calculate fluxes using a linear and non-linear regression following the same methodology used by the smart chamber. The deadband can also be optimized within a certain range.

## Installation

You can install the development version of n2oflux from [GitHub](https://github.com/) with:


``` r
# install.packages("devtools")
install_github("antodi/n2oflux")
```

## Example

``` r
#load package
library(n2oflux)

#store path
PATH<- c("C:/Users/YOURNAME/Documents/file.json")

#process json file
n2o_obs <- process_json_files(PATH)

#calculate fluxes
n2o_flux <- calculate_n2o_flux(n2o_obs,deadband=30,deadband_c=0,stop_time_ag=120,offset_k="json",opt_db="no")

#show results
head(n2o_flux$data_n2o)
```

