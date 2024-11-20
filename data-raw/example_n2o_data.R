## code to prepare `example_n2o_data` dataset goes here

set.seed(1)

example_n2o_data <- data.frame(
  "ETIME" = c(1:120),
  "N2O_DRY" = c(seq(370,375,c((375-370)/119))+runif(120,-0.2,0.2))
  )

usethis::use_data(example_n2o_data)
