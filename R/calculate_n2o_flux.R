#' Calculate N2O fluxes from processed json files
#'
#' This function calculate N2O fluxes (nmol m-2 s-1) using the table created by process_json_files function using a linear and nonlinear regression following the method described by Licor in the 8200 smart chamber manual. When the nonlinear function cannot be fit, the flux value for the linear regression is attributed by default.
#' It also allows to optimize the deadband within a given range. Before steady mixing is attained, N2O concentration can decreases. The optimum deadband is the time (opt_db) at which N2O concentration ceases to decrease and starts increasing. The optimum deadband is found by fitting two linear regression at start-opt_db and opt_db-end with opt_db taking all values between the given range. The opt_db value maximizing the R2 for both regressions is the optimum deadband. In case both regressions are positive, opt_db defaults to the given deadband.
#'
#' @param data The output from the process_json_files function.
#' @param deadband The desired value (sec) for the deadband. N2O observations before that time will be discarded.
#' @param deadband_c Use only if the json file contains CO2 data. The desired value (sec) for the deadband. CO2 observations before that time will be discarded. A value of 0 (default) will ignore this argument.
#' @param stop_time_ag The duration of observation (sec).
#' @param offset_k The value of the offset (cm) for the collar. If set to "json", the offset will be retrieved from the json file.
#' @param opt_db The range at which the deadband should be optimized. For intance, if opt_db is set to '20-50', the optimum deadband will tested for all values between 20 and 50sec. The deadband won't be optimzed if set to 'no' (default).
#' @param clean Value (clean) used in the following formula to detect outliers and discard them. Outliers are defined as observations (obs_i) for which: obs_i > upper_quantile_0.75 + Inter_Quantile x (clean) OR obs_i < lower_quantile_0.25 - Inter_Quantile x (clean). Those observations are discarded. Outliers are not removed by default.
#' @param show_bar Allow to display a progress bar. If set to "yes", a bar will display progress. Is set to "no" by default.
#'
#' @return A list containing; (1) data_n2o, the calculated fluxes. The column 'deadband' will show the optimized deadband if requested. The dNdt, flux, R squared, and RMSE values for the linear and nonlinear regressions are reported in the column containing _LIN (for linear) or _nLIN (for nonlinear regression) in their header. The pvalue for the linear model and for the alpha term of the non-linear model are reported. The column F_N2O reports the flux for the best model based on the RMSE. If the alpha term of the non-linear model is not significant, the flux estimated based on the linear model is given. (2) linear_model_n2o, the results for the linear regression. (3) nonlinear_model_n2o, the results for the nonlinear regression. (4) When outliers were asked to be removed, a list of the measurements where outliers were detected is returned. All observations for those measurements are saved and the discarded outliers are flagged.
#' @export
#'
#' @examples
#' head(example_n2o_data)
#' n2o_flux<-calculate_n2o_flux(example_n2o_data,opt_db="20-50",clean=2)
#' n2o_flux$data_n2o #table with the results
#' names(n2o_flux$outliers) #observations containing outliers, empty in this example
#' names(n2o_flux$linear_model_n2o) #list of linear models for all observations
calculate_n2o_flux <- function(data,deadband=30,deadband_c=0,stop_time_ag=120,offset_k="json",opt_db="no",clean="no", show_bar="no"){

  #set objects
  groups <- unique(interaction(data$date,data$LABEL))

  db_default <- deadband

  res_n2o_flux <- list()
  data_n2o<-c()
  linear_model_n2o <-list()
  nonlinear_model_n2o <-list()
  outliers<-list()

  if(show_bar == "yes"){
    # Initializes the progress bar
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = length(groups), # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
    }

  for (i in 1:length(groups) ){

    if(show_bar=="yes"){
      #progress bar
      setTxtProgressBar(pb, i)
    }

    #subset plot data
    date <- strsplit(as.character(groups[i]),"\\.")[[1]][1]
    LABEL <- strsplit(as.character(groups[i]),"\\.")[[1]][2]
    sub_sam <- data[which(data$date == date & data$LABEL == LABEL ),]

    #if bug in data collection, then return 99999 and proceed to next measurements
    if(length(unique(sub_sam$ETIME)) < 10 ){
      if(deadband_c>0){
        #compile table
        plot_data <- data.frame(date, "time"=99999, LABEL,  "stop_time"=99999, "deadband"=99999, "offset"=99999,
                                "DIAGNOSIS"=99999,"N2O_log_repeats"=99999, "Remark"=99999,
                                "TA_m"=99999, "TS1_m"=99999, "EC2_m"=99999, "SWC2_m"=99999, "TS2_m"=99999,
                                "FN2O_DRY_LIN_dNdt"=99999, "FN2O_DRY_LIN"=99999, "FN2O_DRY_LIN_R2"=99999, "FN2O_DRY_LIN_RMSE"=99999,
                                "FN2O_LIN_pval"=99999,
                                "FN2O_DRY_nLIN_dNdt0"=99999, "FN2O_DRY_nLIN"=99999,"FN2O_DRY_nLIN_R2"=99999, "FN2O_DRY_nLIN_RMSE"=99999,
                                "Cx"=99999,"alpha_v"=99999,"FN2O_alpha_pval"=99999 ,"ETIME0"=99999,"N2O_CV"=99999,
                                "FCO2_DRY_LIN_dNdt"=99999, "FCO2_DRY_LIN"=99999, "FCO2_DRY_LIN_R2"=99999, "FCO2_DRY_LIN_RMSE"=99999,
                                "FCO2_LIN_pval"=99999,
                                "FCO2_DRY_nLIN_dNdt0"=99999, "FCO2_DRY_nLIN"=99999, "FCO2_DRY_nLIN_R2"=99999, "FCO2_DRY_nLIN_RMSE"=99999,
                                "Cx_co2"=99999,"alpha_co2"=99999,"FCO2_alpha_pval"=99999 ,"ETIME0_co2"=99999,"CO2_CV"=99999)

      }else{

        #compile table
        plot_data <- data.frame(date, "time"=99999, LABEL,  "stop_time"=99999, "deadband"=99999, "offset"=99999,
                                "DIAGNOSIS"=99999,"N2O_log_repeats"=99999,"Remark"=99999,
                                "TA_m"=99999, "TS1_m"=99999, "EC2_m"=99999, "SWC2_m"=99999, "TS2_m"=99999,
                                "FN2O_DRY_LIN_dNdt"=99999, "FN2O_DRY_LIN"=99999, "FN2O_DRY_LIN_R2"=99999, "FN2O_DRY_LIN_RMSE"=99999,
                                "FN2O_LIN_pval"=99999,
                                "FN2O_DRY_nLIN_dNdt0"=99999, "FN2O_DRY_nLIN"=99999,"FN2O_DRY_nLIN_R2"=99999, "FN2O_DRY_nLIN_RMSE"=99999,
                                "Cx"=99999,"alpha_v"=99999,"ETIME0"=99999,"N2O_CV"=99999)
      }

      data_n2o <- rbind(data_n2o, plot_data)

      next }

    #check for technical issues
    n2o_repeats<-1
    for(rep_dia in 2:length(sub_sam$N2O_DRY)){

      tryCatch(  #if parameters cannot be estimated, then use dC/dt from linear regression
        { if(sub_sam$N2O_DRY[rep_dia] == sub_sam$N2O_DRY[c(rep_dia-1)] ){
          n2o_repeats[c(rep_dia)] <- n2o_repeats[c(rep_dia-1)]+1  }else{
            n2o_repeats[c(rep_dia)] <- 1 } } ,
        error=function(e){ n2o_repeats[c(rep_dia)] <- 1 } )

    }
    N2O_log_repeats <- max(n2o_repeats,na.rm = TRUE)

    #adjust ETIME
    if(0 %in% sub_sam$ETIME){sub_sam$ETIME <- sub_sam$ETIME +1 } ###NEW

    #cleaning by removing outliers
    if(clean != "no"){
      tresh<- as.numeric(clean)
      upper <- rollapply(sub_sam$N2O_DRY, width = 30, FUN = quantile, p = 0.75, na.rm = TRUE, align = "center", fill = NA)
      lower <- rollapply(sub_sam$N2O_DRY, width = 30, FUN = quantile, p = 0.25, na.rm = TRUE, align = "center", fill = NA)
      iqr <-rollapply(sub_sam$N2O_DRY, width = 30, FUN = IQR, na.rm = TRUE, align = "center", fill = NA)
      sub_sam$outliers <- 0
      sub_sam$outliers[sub_sam$N2O_DRY < lower - iqr * tresh |
                         sub_sam$N2O_DRY > upper + iqr * tresh] <- 1

      if(sum(sub_sam$outliers,na.rm = TRUE)>0){
        outliers[[as.character(groups[i])]] <- sub_sam
      }

      sub_sam <- sub_sam[which(sub_sam$outliers == 0), ]

    }


    #if number of observation below stop_time, then set stop_time to the number of obs
    if(max(sub_sam$ETIME,na.rm = TRUE)<stop_time_ag){stop_time<-max(sub_sam$ETIME,na.rm = TRUE)}else{stop_time<-stop_time_ag}

    #calculate total volume chamber
    #li870 volume = 33.5 cm3
    #li7820 volume = 28 cm3
    #li7810 volume = 28 cm3

    if(offset_k == "json"){
      offset <- unique(sub_sam$Offset) } else{ offset <- offset} #cm

    #optimize deadband
    if(length(sub_sam$N2O_DRY>90) & table(is.na(sub_sam$N2O_DRY))[["FALSE"]] > 89){
      if(str_detect(opt_db,"[0-9]")==TRUE){
        min_db <- strsplit(opt_db,"-")[[1]][1]
        max_db <- strsplit(opt_db,"-")[[1]][2]
        res_tab <-c()
        coef_reg <-c()
        for(j in min_db:max_db){

          sub_sam_1 <- sub_sam[c(1:j),]
          sub_sam_2 <- sub_sam[c(j+1:nrow(sub_sam)),]

          lm_beg <- lm(N2O_DRY~ETIME,data=sub_sam_1)
          lm_end <- lm(N2O_DRY~ETIME,data=sub_sam_2)

          beg <- round(summary(lm_beg)$r.squared,3)
          end <- round(summary(lm_end)$r.squared,3)

          res <- data.frame(beg,end,"deadband"=j)
          res_tab <- rbind(res_tab,res)

        }

        res_tab$opt <- res_tab$beg+res_tab$end

        #pick value with max R2 for both regression
        deadband <- res_tab$deadband[which.max(res_tab$opt)] } #sec

      #check if first regression is positive or negative
      sub_sam_1 <- sub_sam[c(1:deadband),]
      lm_beg <- lm(N2O_DRY~ETIME,data=sub_sam_1)

      if( coef(lm_beg)[[2]] > 0){ deadband <- db_default }
    }else{deadband<-deadband}

    Scham <- unique(sub_sam$Area) #cm2
    ChamVolume <- unique(sub_sam$ChamVolume) #cm3
    IrgaVolume <- unique(sub_sam$IrgaVolume) # = analyzer volume (cm3) + tubing (cm) * tubing area (0.158 cm2)
    Vcham<- Scham*offset + IrgaVolume + ChamVolume #cm3

    if("CO2_DRY" %in% colnames(sub_sam)){Vcham <- Vcham+unique(sub_sam$IrgaVolume_li870) }

    #set constant
    R <- 8.314 #gas constant (Pa m3 K-1 mol-1)

    #get time measurement begins
    time <- strsplit(sub_sam$DATE_TIME[1]," ")[[1]][2]

    #get remark
    Remark <- unique(sub_sam$Remark)

    #check if first 10 sec are usable
    test<-c()
    if(length(unique(is.na(sub_sam$ETIME[1:10]))) ==1 ){
      if(unique(is.na(sub_sam$ETIME[1:10]))!=TRUE) {test <- "pass"}else{
        test <- "fail" }
    }


    if(length(unique(is.na(sub_sam$ETIME[1:10]))) ==2 ){
      if(table(is.na(sub_sam$ETIME[1:10]))["FALSE"]>6 ){test <- "pass"}else{
        test <- "fail" }
    }

    if(test=="pass"){

      #get initial values
      P  <- summary(lm(data=sub_sam[1:10,], PA~ ETIME))$coef[1,1]
      W0 <- summary(lm(data=sub_sam[1:10,], H2O~ ETIME))$coef[1,1]
      T0 <- summary(lm(data=sub_sam[1:10,], TA~ ETIME))$coef[1,1]
      C0 <- summary(lm(data=sub_sam[1:10,], N2O_DRY~ ETIME))$coef[1,1]



      #trim for n2o
      sub_sam <- sub_sam[which(sub_sam$ETIME %in% c(c(1:stop_time)) ),] # select observation length
      sub_sam <- sub_sam[-which(sub_sam$ETIME %in% c(c(1:deadband))  ) ,] # remove deadband

      #extract chamber condition
      TA_m <- mean(sub_sam$TA[which(sub_sam$TA<100)])
      TS1_m <- mean(sub_sam$TS_1[which(sub_sam$TS_1<100)])
      EC2_m <- mean(sub_sam$EC_2[which(sub_sam$EC_2<1)])
      SWC2_m <- mean(sub_sam$SWC_2[which(sub_sam$SWC_2<1)])
      TS2_m <- mean(sub_sam$TS_2[which(sub_sam$TS_2<100)])

      #extract diagnosis
      DIAGNOSIS <- mean(sub_sam$DIAGNOSIS)

      #extract CV
      N2O_CV <- mean(rollapply(sub_sam$N2O_DRY, width = 30, FUN = function(x)
      {sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)}, align = "center", fill = NA)*100, na.rm = TRUE)

      #calculate linear N2O flux
      l_model <- lm(data=sub_sam, N2O_DRY~ ETIME)
      linear_model_n2o[[i]]<-l_model
      names(linear_model_n2o)[i]<-as.character(groups[i])

      FN2O_DRY_LIN_R2 <- summary(l_model)$r.squared
      FN2O_DRY_LIN_RMSE <- sqrt(mean((sub_sam$N2O_DRY - predict(l_model,x=sub_sam$ETIME))^2))

      FN2O_DRY_LIN_dNdt <- summary(l_model)$coef[2,1]

      FN2O_DRY_LIN <- (10*Vcham*P*(1-W0/1000))/(R*Scham*(T0+273.15))* FN2O_DRY_LIN_dNdt

      lin_f <- summary(l_model)$fstatistic
      FN2O_LIN_pval <- pf(lin_f[1],lin_f[2],lin_f[3],lower.tail=F)

      # #calculate non-linear N2O flux
      est_nlin<-list()
      est_nlin_mod<-list()

      tryCatch(  #if parameters cannot be estimated, then use dC/dt from linear regression
        {
          all_par<-c()
          res_fm2 <- nls2(N2O_DRY ~ Cx + (C0-Cx)*exp( -alpha_v*(ETIME-ETIME0)), #get better starting values for nls
                          start = list(Cx=c(1,1000), alpha_v=c(-0.1,0.1),ETIME0=c(5,50)),
                          alg = "brute",data=sub_sam)


          nl_model <- nlsLM(N2O_DRY ~ Cx + (C0-Cx)*exp( -alpha_v*(ETIME-ETIME0)),
                            start=coef(res_fm2),
                            control = nls.lm.control(maxiter=150),
                            data = sub_sam)
          FN2O_DRY_nLIN_R2 <- 1 - (deviance(nl_model)/sum((sub_sam$N2O_DRY-mean(sub_sam$N2O_DRY))^2))
          FN2O_DRY_nLIN_RMSE <- sqrt(mean((sub_sam$N2O_DRY - predict(nl_model,x=sub_sam$ETIME))^2))

          Cx <- coef(nl_model)[[1]]
          alpha_v <- coef(nl_model)[[2]]
          ETIME0 <- coef(nl_model)[[3]]

          dN_dtp <- alpha_v*(Cx-C0)*exp(-alpha_v*(sub_sam$ETIME[nrow(sub_sam)]-ETIME0 ) )
          FN2O_DRY_nLIN_dNdt0 <- alpha_v*(Cx-C0)  #slope at t=t0

          FN2O_alpha_pval <- summary(nl_model)$parameters['alpha_v','Pr(>|t|)']

          all_par <- data.frame(FN2O_DRY_nLIN_R2,FN2O_DRY_nLIN_RMSE,Cx,alpha_v,ETIME0,
                                FN2O_DRY_nLIN_dNdt0,FN2O_alpha_pval)

          est_nlin_mod[[1]] <- nl_model
          est_nlin[[1]] <-  all_par


        },
        error=function(e){
          FN2O_DRY_nLIN_dNdt0 <<- FN2O_DRY_LIN_dNdt
          FN2O_DRY_nLIN_R2 <<- 0
          FN2O_DRY_nLIN_RMSE <<-99999
          Cx <<- 99999
          alpha_v <<- 99999
          FN2O_alpha_pval <<- 99999
          ETIME0 <<-99999

          all_par <<- data.frame(FN2O_DRY_nLIN_R2,FN2O_DRY_nLIN_RMSE,Cx,alpha_v,ETIME0,
                                 FN2O_DRY_nLIN_dNdt0,FN2O_alpha_pval)

          est_nlin_mod[[1]] <<- NA
          est_nlin[[1]] <<-  all_par

        }
      )

      tryCatch(  #if parameters cannot be estimated, then use dC/dt from linear regression
        {
          all_par<-c()

          nl_model <- nlsLM(N2O_DRY ~ Cx + (C0-Cx)*exp( -alpha_v*(ETIME-ETIME0)),
                            start=list(Cx=c(300), alpha_v=c(-.1),ETIME0=c(100)),
                            control = nls.lm.control(maxiter=150),
                            data = sub_sam)
          FN2O_DRY_nLIN_R2 <- 1 - (deviance(nl_model)/sum((sub_sam$N2O_DRY-mean(sub_sam$N2O_DRY))^2))
          FN2O_DRY_nLIN_RMSE <- sqrt(mean((sub_sam$N2O_DRY - predict(nl_model,x=sub_sam$ETIME))^2))

          Cx <- coef(nl_model)[[1]]
          alpha_v <- coef(nl_model)[[2]]
          ETIME0 <- coef(nl_model)[[3]]

          dN_dtp <- alpha_v*(Cx-C0)*exp(-alpha_v*(sub_sam$ETIME[nrow(sub_sam)]-ETIME0 ) )
          FN2O_DRY_nLIN_dNdt0 <- alpha_v*(Cx-C0)  #slope at t=t0

          FN2O_alpha_pval <- summary(nl_model)$parameters['alpha_v','Pr(>|t|)']

          all_par <- data.frame(FN2O_DRY_nLIN_R2,FN2O_DRY_nLIN_RMSE,Cx,alpha_v,ETIME0,
                                FN2O_DRY_nLIN_dNdt0,FN2O_alpha_pval)
          est_nlin_mod[[2]] <- nl_model
          est_nlin[[2]] <-  all_par

        },
        error=function(e){
          FN2O_DRY_nLIN_dNdt0 <<- FN2O_DRY_LIN_dNdt
          FN2O_DRY_nLIN_R2 <<- 0
          FN2O_DRY_nLIN_RMSE <<-99999
          Cx <<- 99999
          alpha_v <<- 99999
          FN2O_alpha_pval <<- 99999
          ETIME0 <<-99999

          all_par <<- data.frame(FN2O_DRY_nLIN_R2,FN2O_DRY_nLIN_RMSE,Cx,alpha_v,ETIME0,
                                 FN2O_DRY_nLIN_dNdt0,FN2O_alpha_pval)

          est_nlin_mod[[2]] <<- NA
          est_nlin[[2]] <<-  all_par
        }
      )

      tryCatch(  #if parameters cannot be estimated, then use dC/dt from linear regression
        {
          all_par<-c()

          nl_model <- nlsLM(N2O_DRY ~ Cx + (C0-Cx)*exp( -alpha_v*(ETIME-ETIME0)),
                            start=list(Cx=c(300), alpha_v=c(.1),ETIME0=c(100)),
                            control = nls.lm.control(maxiter=150),
                            data = sub_sam)
          FN2O_DRY_nLIN_R2 <- 1 - (deviance(nl_model)/sum((sub_sam$N2O_DRY-mean(sub_sam$N2O_DRY))^2))
          FN2O_DRY_nLIN_RMSE <- sqrt(mean((sub_sam$N2O_DRY - predict(nl_model,x=sub_sam$ETIME))^2))

          Cx <- coef(nl_model)[[1]]
          alpha_v <- coef(nl_model)[[2]]
          ETIME0 <- coef(nl_model)[[3]]

          dN_dtp <- alpha_v*(Cx-C0)*exp(-alpha_v*(sub_sam$ETIME[nrow(sub_sam)]-ETIME0 ) )
          FN2O_DRY_nLIN_dNdt0 <- alpha_v*(Cx-C0)  #slope at t=t0

          FN2O_alpha_pval <- summary(nl_model)$parameters['alpha_v','Pr(>|t|)']

          all_par <- data.frame(FN2O_DRY_nLIN_R2,FN2O_DRY_nLIN_RMSE,Cx,alpha_v,ETIME0,
                                FN2O_DRY_nLIN_dNdt0,FN2O_alpha_pval)
          est_nlin_mod[[3]] <- nl_model
          est_nlin[[3]] <-  all_par

        },
        error=function(e){
          FN2O_DRY_nLIN_dNdt0 <<- FN2O_DRY_LIN_dNdt
          FN2O_DRY_nLIN_R2 <<- 0
          FN2O_DRY_nLIN_RMSE <<-99999
          Cx <<- 99999
          alpha_v <<- 99999
          FN2O_alpha_pval <<- 99999
          ETIME0 <<-99999

          all_par <<- data.frame(FN2O_DRY_nLIN_R2,FN2O_DRY_nLIN_RMSE,Cx,alpha_v,ETIME0,
                                 FN2O_DRY_nLIN_dNdt0,FN2O_alpha_pval)

          est_nlin_mod[[3]] <<- NA
          est_nlin[[3]] <<-  all_par
        }
      )

      best_nlin_fit <- which.min(c(est_nlin[[1]]$FN2O_DRY_nLIN_RMSE,
                                   est_nlin[[2]]$FN2O_DRY_nLIN_RMSE,
                                   est_nlin[[3]]$FN2O_DRY_nLIN_RMSE))

      nl_model <- est_nlin_mod[[best_nlin_fit]]
      FN2O_DRY_nLIN_dNdt0 <- est_nlin[[best_nlin_fit]]$FN2O_DRY_nLIN_dNdt0
      FN2O_DRY_nLIN_R2 <- est_nlin[[best_nlin_fit]]$FN2O_DRY_nLIN_R2
      FN2O_DRY_nLIN_RMSE <- est_nlin[[best_nlin_fit]]$FN2O_DRY_nLIN_RMSE
      Cx <- est_nlin[[best_nlin_fit]]$Cx
      alpha_v <- est_nlin[[best_nlin_fit]]$alpha_v
      FN2O_alpha_pval <- est_nlin[[best_nlin_fit]]$FN2O_alpha_pval
      ETIME0 <- est_nlin[[best_nlin_fit]]$ETIME0

      nonlinear_model_n2o[[i]]<-nl_model
      names(nonlinear_model_n2o)[i]<-as.character(groups[i])

      FN2O_DRY_nLIN <- (10*Vcham*P*(1-W0/1000))/(R*Scham*(T0+273.15))* FN2O_DRY_nLIN_dNdt0

    }else{

      var_n <- c("TA_m", "TS1_m", "EC2_m", "SWC2_m", "TS2_m",
                 "FN2O_DRY_LIN_dNdt", "FN2O_DRY_LIN", "FN2O_DRY_LIN_R2", "FN2O_DRY_LIN_RMSE",
                 "FN2O_LIN_pval",
                 "FN2O_DRY_nLIN_dNdt0", "FN2O_DRY_nLIN","FN2O_DRY_nLIN_R2", "FN2O_DRY_nLIN_RMSE",
                 "Cx","alpha_v","FN2O_alpha_pval","ETIME0")

      for (i in seq(var_n)) assign(var_n[i],9999)

    }


    # calculate CO2 fluxes
    if( deadband_c > 0 ){

      #adjust ETIME
      if(0 %in% sub_sam$ETIME_co2){sub_sam$ETIME_co2 <- sub_sam$ETIME_co2 +1 } ###NEW

      #trim for co2
      sub_sam <- data[which(data$date == date & data$LABEL == LABEL ),]

      sub_sam <- sub_sam[which(sub_sam$ETIME_co2 %in% c(c(1:stop_time)) ),] # select observation length
      sub_sam <- sub_sam[-which(sub_sam$ETIME_co2 %in% c(c(1:deadband_c))  ) ,] # remove deadband

      W0_co2 <- summary(lm(data=sub_sam[1:10,], H2O_co2~ ETIME_co2))$coef[1,1]
      C0_co2 <- summary(lm(data=sub_sam[1:10,], CO2_DRY~ ETIME_co2))$coef[1,1]

      #extract CV
      CO2_CV <- mean(rollapply(sub_sam$CO2_DRY, width = 30, FUN = function(x)
      {sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)}, align = "center", fill = NA)*100, na.rm = TRUE)


      #calculate linear CO2 flux
      l_model <- lm(data=sub_sam, CO2_DRY~ ETIME_co2)
      FCO2_DRY_LIN_R2 <- summary(l_model)$r.squared
      FCO2_DRY_LIN_RMSE <- sqrt(mean((sub_sam$CO2_DRY - predict(l_model,x=sub_sam$ETIME_co2))^2))

      FCO2_DRY_LIN_dNdt <- summary(l_model)$coef[2,1]

      FCO2_DRY_LIN <- (10*Vcham*P*(1-W0_co2/1000))/(R*Scham*(T0+273.15))* FCO2_DRY_LIN_dNdt

      lin_f <- summary(l_model)$fstatistic
      FCO2_LIN_pval <- pf(lin_f[1],lin_f[2],lin_f[3],lower.tail=F)

      # #calculate non-linear CO2 flux
      res_fm3 <- nls2(CO2_DRY ~ Cx_co2 + (C0_co2-Cx_co2)*exp( -alpha_co2*(ETIME_co2-ETIME0_co2)), #get better starting values for nls
                      start = list(Cx_co2=c(1000,5000), alpha_co2=c(-0.1,0.1),ETIME0_co2=c(5,50)), alg = "brute",data=sub_sam)

      tryCatch(  #if parameters cannot be estimated, then use dC/dt from linear regression
        { nl_model = nlsLM(CO2_DRY ~ Cx_co2 + (C0_co2-Cx_co2)*exp( -alpha_co2*(ETIME_co2-ETIME0_co2)),
                           start=coef(res_fm3),
                           control = nls.lm.control(maxiter=100),
                           data = sub_sam)


        FCO2_DRY_nLIN_R2 <- 1 - (deviance(nl_model)/sum((sub_sam$CO2_DRY-mean(sub_sam$CO2_DRY))^2))
        FCO2_DRY_nLIN_RMSE <- sqrt(mean((sub_sam$CO2_DRY - predict(nl_model,x=sub_sam$ETIME_co2))^2))

        Cx_co2 <- coef(nl_model)[[1]]
        alpha_co2 <- coef(nl_model)[[2]]
        ETIME0_co2 <- coef(nl_model)[[3]]

        dN_dtp <- alpha_co2*(Cx_co2-C0_co2)*exp(-alpha_co2*(sub_sam$ETIME_co2[nrow(sub_sam)]-ETIME0_co2 ) )
        FCO2_DRY_nLIN_dNdt0 <- alpha_co2*(Cx_co2-C0_co2)  #slope at t=t0

        FCO2_alpha_pval <- summary(nl_model)$parameters['alpha_co2','Pr(>|t|)']

        },
        error=function(e){
          FCO2_DRY_nLIN_dNdt0 <<- FCO2_DRY_LIN_dNdt
          FCO2_DRY_nLIN_R2 <<- 0
          FCO2_DRY_nLIN_RMSE <<-99999
          Cx_co2 <<- 99999
          FCO2_alpha_pval<<-99999
          alpha_co2 <<- 99999
          ETIME0_co2 <<-99999
        }
      )

      FCO2_DRY_nLIN <- (10*Vcham*P*(1-W0_co2/1000))/(R*Scham*(T0+273.15))* FCO2_DRY_nLIN_dNdt0

      #compile table
      plot_data <- data.frame(date, time, LABEL,  stop_time, deadband, offset,
                              DIAGNOSIS, N2O_log_repeats, Remark,
                              TA_m, TS1_m, EC2_m, SWC2_m, TS2_m,
                              FN2O_DRY_LIN_dNdt, FN2O_DRY_LIN, FN2O_DRY_LIN_R2, FN2O_DRY_LIN_RMSE,
                              FN2O_LIN_pval,
                              FN2O_DRY_nLIN_dNdt0, FN2O_DRY_nLIN,FN2O_DRY_nLIN_R2, FN2O_DRY_nLIN_RMSE,
                              Cx,alpha_v,FN2O_alpha_pval,ETIME0,N2O_CV,
                              FCO2_DRY_LIN_dNdt, FCO2_DRY_LIN, FCO2_DRY_LIN_R2, FCO2_DRY_LIN_RMSE,
                              FCO2_LIN_pval,
                              FCO2_DRY_nLIN_dNdt0, FCO2_DRY_nLIN, FCO2_DRY_nLIN_R2, FCO2_DRY_nLIN_RMSE,
                              Cx_co2,alpha_co2,FCO2_alpha_pval,ETIME0_co2,CO2_CV)

    }else{

      #compile table
      plot_data <- data.frame(date, time, LABEL, stop_time, deadband, offset,
                              DIAGNOSIS, N2O_log_repeats, Remark,
                              TA_m, TS1_m, EC2_m, SWC2_m, TS2_m,
                              FN2O_DRY_LIN_dNdt, FN2O_DRY_LIN, FN2O_DRY_LIN_R2, FN2O_DRY_LIN_RMSE,
                              FN2O_LIN_pval,
                              FN2O_DRY_nLIN_dNdt0, FN2O_DRY_nLIN,FN2O_DRY_nLIN_R2, FN2O_DRY_nLIN_RMSE,
                              Cx,alpha_v,FN2O_alpha_pval,ETIME0,N2O_CV )
    }




    data_n2o <- rbind(data_n2o, plot_data)

  }

  #select best fit based on RMSE
  data_n2o$F_N2O<-c()
  data_n2o$best_model<-c()
  for( i in 1:nrow(data_n2o)){

    if(data_n2o$FN2O_DRY_LIN_RMSE[i] <= data_n2o$FN2O_DRY_nLIN_RMSE[i] |
       data_n2o$FN2O_alpha_pval[i] > 0.05){
      data_n2o$F_N2O[i] <- data_n2o$FN2O_DRY_LIN[i]
      data_n2o$best_model[i] <- "LIN"
    }else{
      data_n2o$F_N2O[i] <- data_n2o$FN2O_DRY_nLIN[i]
      data_n2o$best_model[i] <- "nLIN"
    }
  }
  data_n2o$best_model[which(data_n2o$FN2O_DRY_nLIN_RMSE==99999)]<-"LIN"


  if( deadband_c > 0 ){

    #select best fit based on RMSE
    data_n2o$F_CO2<-c()
    data_n2o$best_model_co2<-c()
    for( i in 1:nrow(data_n2o)){

      if(data_n2o$FCO2_DRY_LIN_RMSE[i] <= data_n2o$FCO2_DRY_nLIN_RMSE[i] |
         data_n2o$FCO2_alpha_pval[i]  > 0.05){
        data_n2o$F_CO2[i] <- data_n2o$FCO2_DRY_LIN[i]
        data_n2o$best_model_co2[i] <- "LIN"
      }else{
        data_n2o$F_CO2[i] <- data_n2o$FCO2_DRY_nLIN[i]
        data_n2o$best_model_co2[i] <- "nLIN"
      }
    }
    data_n2o$best_model[which(data_n2o$FCO2_DRY_nLIN_RMSE==99999)]<-"LIN"

  }

  res_n2o_flux$data_n2o <- data_n2o
  res_n2o_flux$linear_model_n2o<-linear_model_n2o
  res_n2o_flux$nonlinear_model_n2o<-nonlinear_model_n2o
  res_n2o_flux$outliers <- outliers

  if(show_bar=="yes"){close(pb) }#close progress bar

  return(res_n2o_flux)
}
