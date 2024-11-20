process_json_files = function(filepath) {

  data_n2o_obs<-c()

  for (fi in 1:length(filepath)){

    #set variables
    plot<-c()
    Date<-c()
    DeadBand<-c()
    Area <- c()
    Offset<-c()
    IrgaVolume<-c()
    ChamVolume<-c()

    F_0<-c()
    t_0<-c()
    C_0<-c()
    alpha_v<-c()
    C_x<-c()


    tims_dat<-c()
    tims_f<-c()
    tims_line<-"false"

    champ_dat<-c()
    champ_f<-c()
    champ_line<-"false"

    champt_dat<-c()
    champt_f<-c()
    champt_line<-"false"

    chamt_dat<-c()
    chamt_f<-c()
    chamt_line<-"false"

    soilt_dat<-c()
    soilt_f<-c()
    soilt_line<-"false"

    soilpc_dat<-c()
    soilpc_f<-c()
    soilpc_line<-"false"

    soilpm_dat<-c()
    soilpm_f<-c()
    soilpm_line<-"false"

    soilpt_dat<-c()
    soilpt_f<-c()
    soilpt_line<-"false"

    n2o_dat<-c()
    n2o_f<-c()
    n2o_line<-"false"

    h2o_dat<-c()
    h2o_f<-c()
    h2o_line<-"false"

    err_dat<-c()
    err_f<-c()
    err_line<-"false"

    Remark<-c()

    InstrumentModel<-"unknown"
    two_instruments <- "unknown"

    IrgaVolume_li870<-c()

    tims_co2_dat<-c()
    tims_co2_f<-c()
    tims_co2_line<-"false"

    F_0_co2<-c()
    t_0_co2<-c()
    C_0_co2<-c()
    alpha_co2<-c()
    C_x_co2<-c()

    co2_dat<-c()
    co2_f<-c()
    co2_line<-"false"

    co2_wet_dat<-c()
    co2_wet_f<-c()
    co2_wet_line<-"false"

    h2o_co2_dat<-c()
    h2o_co2_f<-c()
    h2o_co2_line<-"false"


    data<-"false"
    header<-"false"

    #read files
    con = file(filepath[fi], "r")

    while ( TRUE ) {
      line = readLines(con, n = 1)

      #stop when no more lines
      if ( length(line) == 0 ) {
        break
      }


      #store metadata
      if(str_detect(line, "header\":")==TRUE ){
        header<- "reading"
      }

      if(str_detect(line, "^\\{.")==TRUE & str_detect(line, ":$")==TRUE){

        #if log with no fluxes
        if(length(plot)>length(F_0) & length(plot) == length(Date) ){

          plot <- plot[-length(plot)]
          Date<-Date[-length(Date)]
          DeadBand<-DeadBand[-length(DeadBand)]
          Area<-Area[-length(Area)]
          Offset<-Offset[-length(Offset)]
          ChamVolume<-ChamVolume[-length(ChamVolume)]
          IrgaVolume<-IrgaVolume[-length(IrgaVolume)]
          Remark<-Remark[-length(Remark)]

          tims_dat<-tims_dat[-length(tims_dat)]
          champ_dat<-champ_dat[-length(champ_dat)]
          champt_dat<-champt_dat[-length(champt_dat)]
          chamt_dat<-chamt_dat[-length(chamt_dat)]
          soilt_dat<-soilt_dat[-length(soilt_dat)]
          soilpc_dat<-soilpc_dat[-length(soilpc_dat)]
          soilpm_dat<-soilpm_dat[-length(soilpm_dat)]
          soilpt_dat<-soilpt_dat[-length(soilpt_dat)]

        }

        #if log with no data at all
        if(length(plot)>length(F_0) & length(plot) > length(Date) ){
          plot <- plot[-length(plot)]
          Remark<-Remark[-length(Remark)]

        }

        plot <- c(plot, str_extract(line, "[a-zA-Z0-9_\\-.]+"))
      }


      if(str_detect(line, "remark")==TRUE & (length(plot)>0 ) ){
        rem_ext <- str_extract(line, "(k\":)(.+)(,)",group=2)
        rem_ext <- gsub("\\\"","",rem_ext)
        Remark <- c(Remark, rem_ext )
      }


      if(str_detect(line, "InstrumentModel")==TRUE & (length(plot)>0 ) ){
        InstrumentModel <- str_extract(line, "(Model\":\")(.+)(\")",group=2)

        if(InstrumentModel == "li870"){
          two_instruments <- "TRUE"
        }

      }

      ## compile data from li7820
      if(InstrumentModel == "LI-7820" ){

        if(str_detect(line, "Date\":")==TRUE & header=="reading" ){
          line<-gsub("\"","",line)
          line<-gsub(",","",line)
          Date <- c(Date, str_extract(line, "[0-9].+"))
        }

        if(str_detect(line, "DeadBand\":")==TRUE ){
          DeadBand <- c(DeadBand, str_extract(line, "[0-9.]+"))
        }

        if(str_detect(line, "Area\":")==TRUE ){
          Area<- c(Area, str_extract(line, "[0-9.]+"))
        }

        if(str_detect(line, "Offset\":")==TRUE ){
          Offset<- c(Offset, str_extract(line, "[0-9.]+"))
        }

        if(str_detect(line, "ChamVolume\":")==TRUE ){
          ChamVolume<- c(ChamVolume, str_extract(line, "[0-9.]+"))
        }

        if(str_detect(line, "IrgaVolume\":")==TRUE ){
          IrgaVolume<- c(IrgaVolume, str_extract(line, "[0-9.]+"))
        }


        #store nonlinear reg parameters
        if(str_detect(line, "F_o\":")==TRUE ){
          F_0<- c(F_0, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }

        if(str_detect(line, "t_o\":")==TRUE ){
          t_0<- c(t_0, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }

        if(str_detect(line, "C_o\":")==TRUE ){
          C_0<- c(C_0, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }

        if(str_detect(line, "\"a\":")==TRUE ){
          alpha_v<- c(alpha_v, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }

        if(str_detect(line, "C_x\":")==TRUE ){
          C_x<- c(C_x, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }



        #check position in the document
        if(str_detect(line, "data\":")==TRUE ){
          data<- "reading"
          header<- "false"

        }

        if(str_detect(line, "summary\":")==TRUE ){
          data<- "false"
        }


        #store timestamp data
        if( tims_line == "reading" & str_detect(line, "\\]")==TRUE ){
          tims_dat <- c(tims_dat,paste0(tims_f,str_extract(line, "[0-9,.]+") ) )
          tims_line <- "false"
        }

        if( tims_line == "reading"){
          tims_f <- paste0(tims_f,str_extract(line, "[0-9,.]+") )
        }

        if(str_detect(line, "timestamp\":")==TRUE & data=="reading" ){
          tims_line <- "reading"
          tims_f <- str_extract(line, "[0-9,.]+")
        }

        #extract chamber_p data
        if( champ_line == "reading" & str_detect(line, "\\]")==TRUE ){
          champ_dat <- c(champ_dat,paste0(champ_f,str_extract(line, "[0-9,.-]+") ) )
          champ_line <- "false"
        }

        if( champ_line == "reading"){
          champ_f <- paste0(champ_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "chamber_p\":")==TRUE & data=="reading" ){
          champ_line <- "reading"
          champ_f <- str_extract(line, "[0-9,.-]+")
        }


        #extract chamber_p_t data
        if( champt_line == "reading" & str_detect(line, "\\]")==TRUE ){
          champt_dat <- c(champt_dat,paste0(champt_f,str_extract(line, "[0-9,.-]+") ) )
          champt_line <- "false"
        }

        if( champt_line == "reading"){
          champt_f <- paste0(champt_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "chamber_p_t\":")==TRUE & data=="reading" ){
          champt_line <- "reading"
          champt_f <- str_extract(line, "[0-9,.-]+")
        }


        #extract chamber_t data
        if( chamt_line == "reading" & str_detect(line, "\\]")==TRUE ){
          chamt_dat <- c(chamt_dat,paste0(chamt_f,str_extract(line, "[0-9,.-]+") ) )
          chamt_line <- "false"
        }

        if( chamt_line == "reading"){
          chamt_f <- paste0(chamt_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "chamber_t\":")==TRUE & data=="reading" ){
          chamt_line <- "reading"
          chamt_f <- str_extract(line, "[0-9,.-]+")
        }


        #extract soil_t data
        if( soilt_line == "reading" & str_detect(line, "\\]")==TRUE ){
          soilt_dat <- c(soilt_dat,paste0(soilt_f,str_extract(line, "[0-9,.-]+") ) )
          soilt_line <- "false"
        }

        if( soilt_line == "reading"){
          soilt_f <- paste0(soilt_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "soil_t\":")==TRUE & data=="reading" ){
          soilt_line <- "reading"
          soilt_f <- str_extract(line, "[0-9,.-]+")
        }


        #extract soilp_c data
        if( soilpc_line == "reading" & str_detect(line, "\\]")==TRUE ){
          soilpc_dat <- c(soilpc_dat,paste0(soilpc_f,str_extract(line, "[0-9,.-]+") ) )
          soilpc_line <- "false"
        }

        if( soilpc_line == "reading"){
          soilpc_f <- paste0(soilpc_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "soilp_c\":")==TRUE & data=="reading" ){
          soilpc_line <- "reading"
          soilpc_f <- str_extract(line, "[0-9,.-]+")
        }


        #extract soilp_m data
        if( soilpm_line == "reading" & str_detect(line, "\\]")==TRUE ){
          soilpm_dat <- c(soilpm_dat,paste0(soilpm_f,str_extract(line, "[0-9,.-]+") ) )
          soilpm_line <- "false"
        }

        if( soilpm_line == "reading"){
          soilpm_f <- paste0(soilpm_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "soilp_m\":")==TRUE & data=="reading" ){
          soilpm_line <- "reading"
          soilpm_f <- str_extract(line, "[0-9,.-]+")
        }


        #extract soilp_t data
        if( soilpt_line == "reading" & str_detect(line, "\\]")==TRUE ){
          soilpt_dat <- c(soilpt_dat,paste0(soilpt_f,str_extract(line, "[0-9,.-]+") ) )
          soilpt_line <- "false"
        }

        if( soilpt_line == "reading"){
          soilpt_f <- paste0(soilpt_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "soilp_t\":")==TRUE & data=="reading" ){
          soilpt_line <- "reading"
          soilpt_f <- str_extract(line, "[0-9,.-]+")
        }



        #extract n2o data
        if( n2o_line == "reading" & str_detect(line, "\\]")==TRUE ){
          n2o_dat <- c(n2o_dat,paste0(n2o_f,str_extract(line, "[0-9,.-]+") ) )
          n2o_line <- "false"
        }

        if( n2o_line == "reading"){
          n2o_f <- paste0(n2o_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "n2o\":")==TRUE & data=="reading" ){
          n2o_line <- "reading"
          line<-gsub("n2o","",line)
          n2o_f <- str_extract(line, "[0-9,.-]+")
        }


        #extract h2o data
        if( h2o_line == "reading" & str_detect(line, "\\]")==TRUE ){
          h2o_dat <- c(h2o_dat,paste0(h2o_f,str_extract(line, "[0-9,.-]+") ) )
          h2o_line <- "false"
        }

        if( h2o_line == "reading"){
          h2o_f <- paste0(h2o_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "h2o\":")==TRUE & data=="reading" ){
          line<-gsub("h2o","",line)
          h2o_line <- "reading"
          h2o_f <- str_extract(line, "[0-9,.-]+")
        }

        #extract err data (diagnosis)
        if( err_line == "reading" & str_detect(line, "\\]")==TRUE ){
          err_dat <- c(err_dat,paste0(err_f,str_extract(line, "[0-9,.]+") ) )
          err_line <- "false"
        }

        if( err_line == "reading"){
          err_f <- paste0(err_f,str_extract(line, "[0-9,.]+") )
        }

        if(str_detect(line, "err\":")==TRUE & data=="reading" ){
          line<-gsub("err","",line)
          err_line <- "reading"
          err_f <- str_extract(line, "[0-9,.]+")
        }

      }

      ### compile data from li870
      if(InstrumentModel == "li870"){

        if(str_detect(line, "IrgaVolume\":")==TRUE ){
          IrgaVolume_li870<- c(IrgaVolume_li870, str_extract(line, "[0-9.]+"))
        }

        #store nonlinear reg parameters
        if(str_detect(line, "F_o\":")==TRUE ){
          F_0_co2<- c(F_0_co2, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }

        if(str_detect(line, "t_o\":")==TRUE ){
          t_0_co2<- c(t_0_co2, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }

        if(str_detect(line, "C_o\":")==TRUE ){
          C_0_co2<- c(C_0_co2, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }

        if(str_detect(line, "\"a\":")==TRUE ){
          alpha_co2<- c(alpha_co2, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }

        if(str_detect(line, "C_x\":")==TRUE ){
          C_x_co2<- c(C_x_co2, as.numeric(str_extract(line, "[0-9.e+-]+")))
        }


        #check position in the document
        if(str_detect(line, "data\":")==TRUE ){
          data<- "reading"
          header<- "false"

        }

        if(str_detect(line, "summary\":")==TRUE ){
          data<- "false"
        }


        #extract co2 data
        if( co2_line == "reading" & str_detect(line, "\\]")==TRUE ){
          co2_dat <- c(co2_dat,paste0(co2_f,str_extract(line, "[0-9,.-]+") ) )
          co2_line <- "false"
        }

        if( co2_line == "reading"){
          co2_f <- paste0(co2_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "co2\":")==TRUE & data=="reading" ){
          co2_line <- "reading"
          line<-gsub("co2","",line)
          co2_f <- str_extract(line, "[0-9,.-]+")
        }

        #store timestamp data for co2
        if( tims_co2_line == "reading" & str_detect(line, "\\]")==TRUE ){
          tims_co2_dat <- c(tims_co2_dat,paste0(tims_co2_f,str_extract(line, "[0-9,.]+") ) )
          tims_co2_line <- "false"
        }

        if( tims_co2_line == "reading"){
          tims_co2_f <- paste0(tims_co2_f,str_extract(line, "[0-9,.]+") )
        }

        if(str_detect(line, "timestamp\":")==TRUE & data=="reading" ){
          tims_co2_line <- "reading"
          tims_co2_f <- str_extract(line, "[0-9,.]+")
        }

        #extract co2_wet data
        if( co2_wet_line == "reading" & str_detect(line, "\\]")==TRUE ){
          co2_wet_dat <- c(co2_wet_dat,paste0(co2_wet_f,str_extract(line, "[0-9,.-]+") ) )
          co2_wet_line <- "false"
        }

        if( co2_wet_line == "reading"){
          co2_wet_f <- paste0(co2_wet_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "co2_wet\":")==TRUE & data=="reading" ){
          co2_wet_line <- "reading"
          line<-gsub("co2_wet","",line)
          co2_wet_f <- str_extract(line, "[0-9,.-]+")
        }


        #extract h2o data
        if( h2o_co2_line == "reading" & str_detect(line, "\\]")==TRUE ){
          h2o_co2_dat <- c(h2o_co2_dat,paste0(h2o_co2_f,str_extract(line, "[0-9,.-]+") ) )
          h2o_co2_line <- "false"
        }

        if( h2o_co2_line == "reading"){
          h2o_co2_f <- paste0(h2o_co2_f,str_extract(line, "[0-9,.-]+") )
        }

        if(str_detect(line, "h2o\":")==TRUE & data=="reading" ){
          line<-gsub("h2o","",line)
          h2o_co2_line <- "reading"
          h2o_co2_f <- str_extract(line, "[0-9,.-]+")
        }

      }
    }

    close(con)

    #check if last log is bogus
    if(length(plot) == length(Remark) & length(plot) > length(F_0)){
      plot <- plot[-length(plot)]
      Remark<-Remark[-length(Remark)]
    }


    if(two_instruments == "TRUE" ){

      #compile data
      metadata<-rbind(plot,Date,DeadBand,Area,Offset,ChamVolume,IrgaVolume,IrgaVolume_li870,Remark)
      nonlinear_para<-rbind(F_0,t_0,C_0,alpha_v,C_x,F_0_co2,t_0_co2,C_0_co2,alpha_co2,C_x_co2)
      all_n2o <- list("metadata"=metadata,
                      "nonlinear_para"=nonlinear_para,
                      "ETIME"=tims_dat,
                      "PA"=champ_dat,"T_PA"=champt_dat,"TA"=chamt_dat,
                      "TS_1"=soilt_dat,"EC_2"=soilpc_dat,"SWC_2"=soilpm_dat,"TS_2"=soilpt_dat,
                      "H2O"=h2o_dat,"N2O_DRY"=n2o_dat, "DIAGNOSIS"=err_dat,
                      "ETIME_co2"=tims_co2_dat,
                      "CO2_DRY"=co2_dat,"CO2_WET"=co2_wet_dat,"H2O_co2"=h2o_co2_dat )

      #create a table with row.names matching etime#

      #process compiled data
      data_n2o_json<-c()
      for ( i in 1:ncol(all_n2o$metadata)){


        n2otab<- data.frame("ETIME"=strsplit(all_n2o$ETIME[i],",")[[1]],
                            "PA"=strsplit(all_n2o$PA[i],",")[[1]],
                            "T_PA"=strsplit(all_n2o$T_PA[i],",")[[1]],
                            "TA"=strsplit(all_n2o$TA[i],",")[[1]],
                            "TS_1"=strsplit(all_n2o$TS_1[i],",")[[1]],
                            "EC_2"=strsplit(all_n2o$EC_2[i],",")[[1]],
                            "SWC_2"=strsplit(all_n2o$SWC_2[i],",")[[1]],
                            "TS_2"=strsplit(all_n2o$TS_2[i],",")[[1]],
                            "H2O"=strsplit(all_n2o$H2O[i],",")[[1]],
                            "N2O_DRY"=strsplit(all_n2o$N2O_DRY[i],",")[[1]],
                            "DIAGNOSIS"=strsplit(all_n2o$DIAGNOSIS[i],",")[[1]],
                            row.names = NULL)

        n2otab$ETIME<-as.numeric(n2otab$ETIME)

        #if two etime has the same value, set one half a second earlier
        if(length(names(which(table(n2otab$ETIME)>1)) )>0 ){

          dupl_time <- names(which(table(n2otab$ETIME)>1))

          for(t in 1:length(dupl_time)){
            n2otab$ETIME[which(n2otab$ETIME %in% dupl_time[t])][1] <- n2otab$ETIME[which(n2otab$ETIME %in% dupl_time[t])][1]-0.5
          }
        }

        if(TRUE %in% is.na(n2otab$ETIME) ){

        }else{
          row.names(n2otab)<-n2otab$ETIME
        }


        co2tab<-data.frame("ETIME_co2"=strsplit(all_n2o$ETIME_co2[i],",")[[1]],
                           "F_0_co2"= all_n2o$nonlinear_para["F_0_co2",i],
                           "t_0_co2"= all_n2o$nonlinear_para["t_0_co2",i],
                           "C_0_co2"= all_n2o$nonlinear_para["C_0_co2",i],
                           "alpha_co2"= all_n2o$nonlinear_para["alpha_co2",i],
                           "C_x_co2"= all_n2o$nonlinear_para["C_x_co2",i],
                           "CO2_DRY"=strsplit(all_n2o$CO2_DRY[i],",")[[1]],
                           "CO2_WET"=strsplit(all_n2o$CO2_WET[i],",")[[1]],
                           "H2O_co2"=strsplit(all_n2o$H2O_co2[i],",")[[1]],
                           row.names = NULL)
        co2tab$ETIME<-as.numeric(co2tab$ETIME)

        #if two etime has the same value, set one half a second earlier
        if(length(names(which(table(co2tab$ETIME)>1)) )>0 ){

          dupl_time <- names(which(table(co2tab$ETIME)>1))

          for(t in 1:length(dupl_time)){
            co2tab$ETIME_co2[which(co2tab$ETIME_co2 %in% dupl_time[t])][1] <- co2tab$ETIME[which(co2tab$ETIME_co2 %in% dupl_time[t])][1]-0.5
          }
        }


        if(TRUE %in% is.na(n2otab$ETIME) ){

        }else{
          row.names(co2tab)<-co2tab$ETIME_co2
        }


        #combine based on etime
        tab_anal <- merge(n2otab,co2tab, by=0, all=TRUE)
        colnames(tab_anal)[which(colnames(tab_anal)=="ETIME.x")]<-"ETIME"
        tab_anal <- tab_anal[order(as.numeric(tab_anal$Row.names)),]
        row.names(tab_anal)<-NULL

        tab1 <- data.frame("LABEL"= all_n2o$metadata["plot",i],
                           "DATE_TIME"= all_n2o$metadata["Date",i],
                           "DeadBand"= all_n2o$metadata["DeadBand",i],
                           "Area"= all_n2o$metadata["Area",i],
                           "Offset"= all_n2o$metadata["Offset",i],
                           "ChamVolume"= all_n2o$metadata["ChamVolume",i],
                           "IrgaVolume"= all_n2o$metadata["IrgaVolume",i],
                           "IrgaVolume_li870"= all_n2o$metadata["IrgaVolume_li870",i],
                           "Remark"= all_n2o$metadata["Remark",i],

                           "F_0"= all_n2o$nonlinear_para["F_0",i],
                           "t_0"= all_n2o$nonlinear_para["t_0",i],
                           "C_0"= all_n2o$nonlinear_para["C_0",i],
                           "alpha_v"= all_n2o$nonlinear_para["alpha_v",i],
                           "C_x"= all_n2o$nonlinear_para["C_x",i],
                           tab_anal,
                           row.names = NULL)


        data_n2o_json<-rbind(data_n2o_json,tab1,row.names = NULL)


      }

    }else{

      #compile data
      metadata<-rbind(plot,Date,DeadBand,Area,Offset,ChamVolume,IrgaVolume,Remark)

      nonlinear_para<-rbind(F_0,t_0,C_0,alpha_v,C_x)
      all_n2o <- list("metadata"=metadata,
                      "nonlinear_para"=nonlinear_para,
                      "ETIME"=tims_dat,
                      "PA"=champ_dat,"T_PA"=champt_dat,"TA"=chamt_dat,
                      "TS_1"=soilt_dat,"EC_2"=soilpc_dat,"SWC_2"=soilpm_dat,"TS_2"=soilpt_dat,
                      "H2O"=h2o_dat,"N2O_DRY"=n2o_dat, "DIAGNOSIS"=err_dat)

      #process compiled data
      data_n2o_json<-c()
      for ( i in 1:ncol(all_n2o$metadata)){

        tab1<-data.frame("LABEL"= all_n2o$metadata["plot",i],
                         "DATE_TIME"= all_n2o$metadata["Date",i],
                         "DeadBand"= all_n2o$metadata["DeadBand",i],
                         "Area"= all_n2o$metadata["Area",i],
                         "Offset"= all_n2o$metadata["Offset",i],
                         "ChamVolume"= all_n2o$metadata["ChamVolume",i],
                         "IrgaVolume"= all_n2o$metadata["IrgaVolume",i],
                         "Remark"= all_n2o$metadata["Remark",i],

                         "F_0"= all_n2o$nonlinear_para["F_0",i],
                         "t_0"= all_n2o$nonlinear_para["t_0",i],
                         "C_0"= all_n2o$nonlinear_para["C_0",i],
                         "alpha_v"= all_n2o$nonlinear_para["alpha_v",i],
                         "C_x"= all_n2o$nonlinear_para["C_x",i],

                         "ETIME"=strsplit(all_n2o$ETIME[i],",")[[1]],
                         "PA"=strsplit(all_n2o$PA[i],",")[[1]],
                         "T_PA"=strsplit(all_n2o$T_PA[i],",")[[1]],
                         "TA"=strsplit(all_n2o$TA[i],",")[[1]],
                         "TS_1"=strsplit(all_n2o$TS_1[i],",")[[1]],
                         "EC_2"=strsplit(all_n2o$EC_2[i],",")[[1]],
                         "SWC_2"=strsplit(all_n2o$SWC_2[i],",")[[1]],
                         "TS_2"=strsplit(all_n2o$TS_2[i],",")[[1]],
                         "H2O"=strsplit(all_n2o$H2O[i],",")[[1]],
                         "N2O_DRY"=strsplit(all_n2o$N2O_DRY[i],",")[[1]],
                         "DIAGNOSIS"=strsplit(all_n2o$DIAGNOSIS[i],",")[[1]],
                         row.names = NULL )
        data_n2o_json<-rbind(data_n2o_json,tab1,row.names = NULL)

      }

    }


    #set to numeric
    cols.num <- colnames(data_n2o_json)[-which(colnames(data_n2o_json) %in% c("LABEL","DATE_TIME", "Remark") )]
    data_n2o_json[cols.num] <- sapply(data_n2o_json[cols.num],as.numeric)

    #add date column
    data_n2o_json$date <- gsub("\\D", "",  do.call(rbind,str_split(data_n2o_json$DATE_TIME," "))[,1] )

    #append final table
    data_n2o_obs<-rbind(data_n2o_obs,data_n2o_json)

  }

  return(data_n2o_obs)


}
