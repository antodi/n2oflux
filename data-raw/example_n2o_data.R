## code to prepare `example_n2o_data` dataset goes here

set.seed(1)

example_n2o_data <- rbind (
  data.frame(
    "TIME"="12:00:00",
    "date"="20241120",
    "LABEL"="1",
    "DATE_TIME"= "11/20/2024  12:00:00 PM" ,
    "DeadBand"=30,
    "Area"=318,
    "Offset"=3,
    "ChamVolume"=4244.1,
    "IrgaVolume"=59.6,
    "Remark"="",
    "F_0"=7.66	,
    "t_0"=4.808,
    "C_0"=365,
    "alpha_v"=0.00036,
    "C_x"=4484,
    "ETIME" = c(1:120),
    "PA"=runif(120,98.1,98.2),
    "T_PA"=32,
    "TA"=runif(120,26.4,26.65),
    "TS_1"=23,
    "EC_2"=0.05,
    "SWC_2"=0.282,
    "TS_2"=23.7,
    "H2O"= c(seq(27.5,30,c((30-27.5)/119))+runif(120,-0.2,0.2)),
    "N2O_DRY" =c(c(seq(375.382 ,376.1983  ,c((376.1983 -375.382)/119))+runif(120,-0.1,0.1))) ,
    "DIAGNOSIS"=0
  ),
  data.frame(
    "TIME"="11:00:00",
    "date"="20241120",
    "LABEL"="2",
    "DATE_TIME"= "11/20/2024  11:00:00 PM" ,
    "DeadBand"=30,
    "Area"=318,
    "Offset"=3,
    "ChamVolume"=4244.1,
    "IrgaVolume"=59.6,
    "Remark"="",
    "F_0"=7.66	,
    "t_0"=4.808,
    "C_0"=365,
    "alpha_v"=0.00036,
    "C_x"=4484,
    "ETIME" = c(1:120),
    "PA"=runif(120,98.1,98.2),
    "T_PA"=32,
    "TA"=runif(120,26.4,26.65),
    "TS_1"=23,
    "EC_2"=0.05,
    "SWC_2"=0.282,
    "TS_2"=23.7,
    "H2O"= c(seq(27.5,30,c((30-27.5)/119))+runif(120,-0.2,0.2)),
    "N2O_DRY" =379-4*exp(-0.01*(1:120-20))+runif(120,-0.1,0.1) ,
    "DIAGNOSIS"=0
  ), data.frame(
    "TIME"="11:00:00",
    "date"="20241120",
    "LABEL"="3",
    "DATE_TIME"= "11/20/2024  11:00:00 PM" ,
    "DeadBand"=30,
    "Area"=318,
    "Offset"=3,
    "ChamVolume"=4244.1,
    "IrgaVolume"=59.6,
    "Remark"="",
    "F_0"=7.66	,
    "t_0"=4.808,
    "C_0"=365,
    "alpha_v"=0.00036,
    "C_x"=4484,
    "ETIME" = c(1:120),
    "PA"=runif(120,98.1,98.2),
    "T_PA"=32,
    "TA"=runif(120,26.4,26.65),
    "TS_1"=23,
    "EC_2"=0.05,
    "SWC_2"=0.282,
    "TS_2"=23.7,
    "H2O"= c(seq(27.5,30,c((30-27.5)/119))+runif(120,-0.2,0.2)),
    "N2O_DRY" =c(c(seq(376.382 ,378.1983  ,c((378.1983 -376.382)/39))+runif(40,-0.3,0.3)),
                 c(seq(378.1983 ,379.6983  ,c((379.6983 -378.1983)/39))+runif(40,-0.3,0.3)),
                 c(seq(379.6983  ,381 ,c((381 -379.6983)/39))+runif(40,-0.3,0.3))) ,
    "DIAGNOSIS"=0
  ), data.frame(
    "TIME"="11:00:00",
    "date"="20241120",
    "LABEL"="4",
    "DATE_TIME"= "11/20/2024  11:00:00 PM" ,
    "DeadBand"=30,
    "Area"=318,
    "Offset"=3,
    "ChamVolume"=4244.1,
    "IrgaVolume"=59.6,
    "Remark"="",
    "F_0"=7.66	,
    "t_0"=4.808,
    "C_0"=375.512,
    "alpha_v"=0.00458992,
    "C_x"=379.582,
    "ETIME" = c(1:120),
    "PA"=runif(120,98.1,98.2),
    "T_PA"=34.6354 ,
    "TA"=runif(120,21.4,21.65),
    "TS_1"=17,
    "EC_2"=0.05,
    "SWC_2"=0.282,
    "TS_2"=23.7,
    "H2O"= c(seq(27.5,30,c((30-27.5)/119))+runif(120,-0.2,0.2)),
    "N2O_DRY" =379-5*exp(-0.01*(1:120-40))+runif(120,-0.1,0.1) ,
    "DIAGNOSIS"=0
  )
)

usethis::use_data(example_n2o_data)
