#' Show results of linear and non-linear fits for N2O data, along with deadband and outliers when used
#'
#' @param id The ID (written as date.LABEL) of measurements to plot. Several measurements can be processed at the same time.
#' @param data The data containing the measurements seconds by seconds, given by the process_json_files function.
#' @param deadband The desired value (sec) for the deadband. N2O observations before that time will be discarded.
#' @param stop_time_ag The duration of observation (sec).
#' @param offset_k The value of the offset (cm) for the collar. If set to "json", the offset will be retrieved from the json file.
#' @param opt_db The range at which the deadband should be optimized. For intance, if opt_db is set to '20-50', the optimum deadband will tested for all values between 20 and 50sec. The deadband won't be optimzed if set to 'no' (default).
#' @param clean Value (clean) used in the following formula to detect outliers and discard them. Outliers are defined as observations (obs_i) for which: obs_i > upper_quantile_0.75 + Inter_Quantile x (clean) OR obs_i < lower_quantile_0.25 - Inter_Quantile x (clean). Those observations are discarded. Outliers are not removed by default.
#'
#' @returns A list containing the plots for all IDs
#' @export
#'
#' @examples
#' plot_res <- plot_n2o_flux<-calculate_n2o_flux(id="20250310.plot10",data=n2o_obs,opt_db="20-50",clean=2)
plot_fit <- function(id,data,deadband=30,stop_time_ag=120,offset_k="json",opt_db="no",clean="no"){

  all_fig <- list()
  all_plot <-list()
  db_default <- deadband

  for( mesu in 1:length(id)){

    measurement_id<-id[mesu]

    #subset plot data
    date <- strsplit(as.character(measurement_id),"\\.")[[1]][1]
    LABEL <- strsplit(as.character(measurement_id),"\\.")[[1]][2]
    sub_sam <- data[which(data$date == date & data$LABEL == LABEL ),]

    #get deadband
    ymin <- min(sub_sam$N2O_DRY)
    ymax <- max(sub_sam$N2O_DRY)

    #get offset
    if(offset_k == "json"){
      offset <- unique(sub_sam$Offset) } else{ offset <- offset} #cm

    #adjust ETIME
    if(0 %in% sub_sam$ETIME){sub_sam$ETIME <- sub_sam$ETIME +1 } ###NEW

    sub_sam_full <- sub_sam

    #cleaning by removing outliers

    are_there_outliers <- FALSE

    if(clean != "no"){

      tresh<- as.numeric(clean)
      upper <- rollapply(sub_sam$N2O_DRY, width = 30, FUN = quantile, p = 0.75, na.rm = TRUE, align = "center", fill = NA)
      lower <- rollapply(sub_sam$N2O_DRY, width = 30, FUN = quantile, p = 0.25, na.rm = TRUE, align = "center", fill = NA)
      iqr <-rollapply(sub_sam$N2O_DRY, width = 30, FUN = IQR, na.rm = TRUE, align = "center", fill = NA)
      sub_sam$outliers <- 0
      sub_sam$outliers[sub_sam$N2O_DRY < lower - iqr * tresh |
                         sub_sam$N2O_DRY > upper + iqr * tresh] <- 1

      if( 1 %in% sub_sam$outliers){
        all_outliers <- sub_sam
        are_there_outliers <- TRUE } ## save outliers

      sub_sam <- sub_sam[which(sub_sam$outliers == 0), ]

    }



    #if number of observation below stop_time, then set stop_time to the number of obs
    if(max(sub_sam$ETIME,na.rm = TRUE)<stop_time_ag){stop_time<-max(sub_sam$ETIME,na.rm = TRUE)}else{stop_time<-stop_time_ag}

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



    #get initial values
    P  <- summary(lm(data=sub_sam[1:10,], PA~ ETIME))$coef[1,1]
    W0 <- summary(lm(data=sub_sam[1:10,], H2O~ ETIME))$coef[1,1]
    T0 <- summary(lm(data=sub_sam[1:10,], TA~ ETIME))$coef[1,1]
    C0 <- summary(lm(data=sub_sam[1:10,], N2O_DRY~ ETIME))$coef[1,1]

    #get factors
    R <- 8.314 #gas constant (Pa m3 K-1 mol-1)
    Scham <- unique(sub_sam$Area) #cm2
    ChamVolume <- unique(sub_sam$ChamVolume) #cm3
    IrgaVolume <- unique(sub_sam$IrgaVolume) # = analyzer volume (cm3) + tubing (cm) * tubing area (0.158 cm2)
    Vcham<- Scham*offset + IrgaVolume + ChamVolume #cm3

    # if(deadband_c>0){Vcham <- Vcham+unique(sub_sam$IrgaVolume_li870) }
    if("CO2_DRY" %in% colnames(sub_sam)){Vcham <- Vcham+unique(sub_sam$IrgaVolume_li870) } #BETTER!!!!


    #trim for n2o
    sub_sam <- sub_sam[which(sub_sam$ETIME %in% c(c(1:stop_time)) ),] # select observation length ###NEW
    sub_sam <- sub_sam[-which(sub_sam$ETIME %in% c(c(1:deadband))  ) ,] # remove deadband ###NEW

    #calculate linear N2O flux
    l_model <- lm(data=sub_sam, N2O_DRY~ ETIME)

    FN2O_DRY_LIN_R2 <- summary(l_model)$r.squared
    FN2O_DRY_LIN_RMSE <- sqrt(mean((sub_sam$N2O_DRY - predict(l_model,x=sub_sam$ETIME))^2))

    FN2O_DRY_LIN_dNdt <- summary(l_model)$coef[2,1]

    FN2O_DRY_LIN <- (10*Vcham*P*(1-W0/1000))/(R*Scham*(T0+273.15))* FN2O_DRY_LIN_dNdt


    # #calculate non-linear N2O flux
    res_fm <- nls2(N2O_DRY ~ Cx + (C0-Cx)*exp( -alpha_v*(ETIME-ETIME0)), #get better starting values for nls
                   start = list(Cx=c(350,500), alpha_v=c(-0.1,0.1),ETIME0=c(5,50)),
                   alg = "brute",data=sub_sam)


    tryCatch(  #if parameters cannot be estimated, then use dC/dt from linear regression
      {nl_model <- nlsLM(N2O_DRY ~ Cx + (C0-Cx)*exp( -alpha_v*(ETIME-ETIME0)),
                         start=coef(res_fm),
                         control = nls.lm.control(maxiter=100),
                         data = sub_sam)

      FN2O_DRY_nLIN_R2 <- 1 - (deviance(nl_model)/sum((sub_sam$N2O_DRY-mean(sub_sam$N2O_DRY))^2))
      FN2O_DRY_nLIN_RMSE <- sqrt(mean((sub_sam$N2O_DRY - predict(nl_model,x=sub_sam$ETIME))^2))

      Cx <- coef(nl_model)[[1]]
      alpha_v <- coef(nl_model)[[2]]
      ETIME0 <- coef(nl_model)[[3]]

      dN_dtp <- alpha_v*(Cx-C0)*exp(-alpha_v*(sub_sam$ETIME[nrow(sub_sam)]-ETIME0 ) )
      FN2O_DRY_nLIN_dNdt0 <- alpha_v*(Cx-C0)  #slope at t=t0
      },
      error=function(e){
        FN2O_DRY_nLIN_dNdt0 <<- FN2O_DRY_LIN_dNdt
        FN2O_DRY_nLIN_R2 <<- 0
        FN2O_DRY_nLIN_RMSE <<-0
        Cx <<- 99999
        alpha_v <<- 99999
        ETIME0 <<-99999

      }
    )

    FN2O_DRY_nLIN <- (10*Vcham*P*(1-W0/1000))/(R*Scham*(T0+273.15))* FN2O_DRY_nLIN_dNdt0

    nlin_val <- paste0("alpha=",round(alpha_v,6),", etime0=",round(ETIME0,2),", Cx", "=",round(Cx,1))


    #linear regression
    fig<-ggplot(sub_sam_full,aes(x=ETIME,y=N2O_DRY))

    if(are_there_outliers == TRUE){
      fig <- fig + geom_point(data=all_outliers, size=2,shape=21,stroke=1, aes(color = as.factor(outliers)) )+
        scale_colour_manual(values=c("black","red"))
    }else{ fig <- fig+ geom_point(size=2,shape=21,stroke=1)+ scale_colour_manual(values=c("black"))
    }

    fig <- fig +

      annotate('rect', xmin=0,xmax=deadband, ymin=ymin, ymax=ymin+(ymax-ymin)*1.2, alpha=.2, fill='black')+ #deadband
      annotate("text", label = paste0("deadband=",deadband),
               x = 2, y = ymin, size = 4, colour = "black",hjust=0  ) + #deadband value

      annotate("text", label = measurement_id,
               x = 2, y = ymin+(ymax-ymin)*1.2, size = 6, colour = "black",hjust=0  ) + #ID

      annotate("text", label = paste0("lin: flux=",round(FN2O_DRY_LIN,3),", R2=",round(FN2O_DRY_LIN_R2,3),
                                      " , RMSE=",round(FN2O_DRY_LIN_RMSE,3)),
               x = 2, y = ymin+(ymax-ymin)*1.1, size = 4.5, colour = "blue",hjust=0  ) + #LIN

      annotate("text", label = paste0("nlin: flux=",round(FN2O_DRY_nLIN,3), ", R2=",round(FN2O_DRY_nLIN_R2,3),
                                      " , RMSE=",round(FN2O_DRY_nLIN_RMSE,3)),
               x = 2, y = ymin+(ymax-ymin)*1.0, size = 4.5, colour = "red" ,hjust=0 ) +#nLIN

      annotate("text", label = nlin_val,
               x = 2, y = ymin+(ymax-ymin)*.9, size = 4.5, colour = "red",hjust=0  ) +#nLIN parameters

      geom_vline(xintercept = deadband)+
      geom_smooth(data= sub_sam,method = "lm", formula = y~x)+
      guides(color="none")+
      theme_bw()+
      theme(axis.text=element_text(size=12,color="black"),
            axis.title=element_text(size=14,color="black"),
            panel.border = element_rect(size =  1.5),
            plot.margin = margin(.5, .5, .5, .5, "cm") )

    if(FN2O_DRY_nLIN_RMSE>0){
      fig <- fig +  geom_line(data=sub_sam,aes(ETIME,predict(nl_model) ),
                              inherit.aes = FALSE,color="red",linewidth=1)}

    all_fig[[mesu]]<- ggplotGrob(fig)

    grid::grid.draw(all_fig[[mesu]])
    all_plot[[mesu]] <- recordPlot()
    names(all_plot)[length(all_plot)]<-measurement_id


  }

  return(all_plot)

}
