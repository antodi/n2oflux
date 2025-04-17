#' Show results of linear and non-linear fits for CO2 data, along with deadband and outliers when used.
#'
#' @param id The ID (written as date.LABEL) of measurements to plot. Several measurements can be processed at the same time.
#' @param data The data containing the measurements seconds by seconds, given by the process_json_files function.
#' @param deadband_c The desired value (sec) for the deadband. CO2 observations before that time will be discarded.
#' @param stop_time_ag The duration of observation (sec).
#' @param offset_k The value of the offset (cm) for the collar. If set to "json", the offset will be retrieved from the json file.
#'
#' @returns A list containing the plots for all IDs
#' @export
#'
#' @examples
plot_fit_co2 <- function (id, data, deadband_c = 30, stop_time_ag = 120, offset_k = "json"){

  all_fig <- list()
  all_plot <- list()

  for (mesu in 1:length(id)) {

    measurement_id <- id[mesu]
    date <- strsplit(as.character(measurement_id), "\\.")[[1]][1]
    LABEL <- strsplit(as.character(measurement_id), "\\.")[[1]][2]
    sub_sam <- data[which(data$date == date & data$LABEL == LABEL), ]
    ymin <- min(sub_sam$CO2_DRY)
    ymax <- max(sub_sam$CO2_DRY)

    if (0 %in% sub_sam$ETIME_co2) {
      sub_sam$ETIME_co2 <- sub_sam$ETIME_co2 + 1
    }

    sub_sam_full <- sub_sam

    if (max(sub_sam$ETIME_co2, na.rm = TRUE) < stop_time_ag) {
      stop_time <- max(sub_sam$ETIME_co2, na.rm = TRUE)    } else {
        stop_time <- stop_time_ag
      }

    if (offset_k == "json") {
      offset <- unique(sub_sam$Offset) }else {
        offset <- offset
      }



    sub_sam <- sub_sam[which(sub_sam$ETIME_co2 %in% c(c(1:stop_time))), ]
    sub_sam <- sub_sam[-which(sub_sam$ETIME_co2 %in% c(c(1:deadband_c))), ]

    P <- summary(lm(data = sub_sam[1:10, ], PA ~ ETIME_co2))$coef[1, 1]
    T0 <- summary(lm(data = sub_sam[1:10, ], TA ~ ETIME_co2))$coef[1, 1]

    W0_co2 <- summary(lm(data = sub_sam[1:10, ], H2O_co2 ~ ETIME_co2))$coef[1, 1]
    C0_co2 <- summary(lm(data = sub_sam[1:10, ], CO2_DRY ~ ETIME_co2))$coef[1, 1]
    CO2_CV <- mean(rollapply(sub_sam$CO2_DRY, width = 30,
                             FUN = function(x) {
                               sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
                             }, align = "center", fill = NA) * 100, na.rm = TRUE)

    R <- 8.314
    Scham <- unique(sub_sam$Area)
    ChamVolume <- unique(sub_sam$ChamVolume)
    IrgaVolume <- unique(sub_sam$IrgaVolume)
    Vcham <- Scham * offset + IrgaVolume + ChamVolume

    if ("CO2_DRY" %in% colnames(sub_sam)) {
      Vcham <- Vcham + unique(sub_sam$IrgaVolume_li870)
    }


    l_model <- lm(data = sub_sam, CO2_DRY ~ ETIME_co2)
    FCO2_DRY_LIN_R2 <- summary(l_model)$r.squared
    FCO2_DRY_LIN_RMSE <- sqrt(mean((sub_sam$CO2_DRY - predict(l_model, x = sub_sam$ETIME_co2))^2))
    FCO2_DRY_LIN_dNdt <- summary(l_model)$coef[2, 1]
    FCO2_DRY_LIN <- (10 * Vcham * P * (1 - W0_co2/1000))/(R * Scham * (T0 + 273.15)) * FCO2_DRY_LIN_dNdt

    lin_f <- summary(l_model)$fstatistic

    FCO2_LIN_pval <- pf(lin_f[1], lin_f[2], lin_f[3], lower.tail = F)

    if (FCO2_LIN_pval<0.001){FCO2_LIN_pval<-"<0.001"}else{
      FCO2_LIN_pval <- paste0("=",round(as.numeric(FCO2_LIN_pval),4) ) }

    res_fm3 <- nls2(CO2_DRY ~ Cx_co2 + (C0_co2 - Cx_co2) * exp(-alpha_co2 * (ETIME_co2 - ETIME0_co2)),
                    start = list(Cx_co2 = c(1000, 5000), alpha_co2 = c(-0.1, 0.1), ETIME0_co2 = c(5, 50)),
                    alg = "brute", data = sub_sam)
    tryCatch({
      nl_model = nlsLM(CO2_DRY ~ Cx_co2 + (C0_co2 -
                                             Cx_co2) * exp(-alpha_co2 * (ETIME_co2 - ETIME0_co2)),
                       start = coef(res_fm3), control = nls.lm.control(maxiter = 100),
                       data = sub_sam)
      FCO2_DRY_nLIN_R2 <- 1 - (deviance(nl_model)/sum((sub_sam$CO2_DRY -
                                                         mean(sub_sam$CO2_DRY))^2))
      FCO2_DRY_nLIN_RMSE <- sqrt(mean((sub_sam$CO2_DRY -
                                         predict(nl_model, x = sub_sam$ETIME_co2))^2))
      Cx_co2 <- coef(nl_model)[[1]]
      alpha_co2 <- coef(nl_model)[[2]]
      ETIME0_co2 <- coef(nl_model)[[3]]
      dN_dtp <- alpha_co2 * (Cx_co2 - C0_co2) * exp(-alpha_co2 *
                                                      (sub_sam$ETIME_co2[nrow(sub_sam)] - ETIME0_co2))
      FCO2_DRY_nLIN_dNdt0 <- alpha_co2 * (Cx_co2 -
                                            C0_co2)
      FCO2_alpha_pval <- summary(nl_model)$parameters["alpha_co2",
                                                      "Pr(>|t|)"]
    }, error = function(e) {
      FCO2_DRY_nLIN_dNdt0 <<- FCO2_DRY_LIN_dNdt
      FCO2_DRY_nLIN_R2 <<- 0
      FCO2_DRY_nLIN_RMSE <<- 99999
      Cx_co2 <<- 99999
      FCO2_alpha_pval <<- 99999
      alpha_co2 <<- 99999
      ETIME0_co2 <<- 99999
    })

    FCO2_DRY_nLIN <- (10 * Vcham * P * (1 - W0_co2/1000))/(R * Scham * (T0 + 273.15)) * FCO2_DRY_nLIN_dNdt0

    if (FCO2_alpha_pval < 0.001) {
      FCO2_alpha_pval <- "<0.001"
    }
    else {
      FCO2_alpha_pval <- paste0("=", round(as.numeric(FCO2_alpha_pval),4))
    }

    nlin_val <- paste0("alpha=", round(alpha_co2, 6),
                       ", etime0=",round(ETIME0_co2, 2),
                       ", Cx", "=", round(Cx_co2, 1))

    fig <- ggplot(sub_sam_full, aes(x = ETIME, y = CO2_DRY))
    fig <- fig + geom_point(size = 2, shape = 21, stroke = 1) +
        scale_colour_manual(values = c("black")) +
      annotate('rect', xmin=0,xmax=deadband_c, ymin=ymin, ymax=ymin+(ymax-ymin)*1.2, alpha=.2, fill='black')+ #deadband
      annotate("text", label = paste0("deadband=",deadband_c),
               x = 2, y = ymin, size = 4, colour = "black",hjust=0  )  #deadband value
    fig <- fig + annotate("text",
                          label = measurement_id, x = 2, y = ymin + (ymax - ymin) * 1.2, size = 6, colour = "black", hjust = 0) +
      annotate("text", label = paste0("lin: flux=", round(FCO2_DRY_LIN, 3),
                                      ", R2=", round(FCO2_DRY_LIN_R2, 3), " , RMSE=",
                                      round(FCO2_DRY_LIN_RMSE, 3), ", pval", FCO2_LIN_pval),
               x = 2, y = ymin + (ymax - ymin) * 1.1, size = 4.5, colour = "blue", hjust = 0) +
      annotate("text", label = paste0("nlin: flux=", round(FCO2_DRY_nLIN, 3), ", R2=", round(FCO2_DRY_nLIN_R2, 3),
                                      " , RMSE=", round(FCO2_DRY_nLIN_RMSE, 3), ", pval", FCO2_alpha_pval),
                                                      x = 2, y = ymin + (ymax - ymin) * 1, size = 4.5,
                                                      colour = "red", hjust = 0) +
      annotate("text", label = nlin_val, x = 2, y = ymin + (ymax - ymin) * 0.9, size = 4.5, colour = "red", hjust = 0) +
      geom_vline(xintercept = deadband_c) +
      geom_smooth(data = sub_sam, method = "lm", formula = y ~ x) + guides(color = "none") +
      theme_bw() + theme(axis.text = element_text(size = 12, color = "black"),
                         axis.title = element_text(size = 14, color = "black"),
                         panel.border = element_rect(linewidth = 1.5) )

    if (FCO2_DRY_nLIN_RMSE != 99999) {
      fig <- fig + geom_line(data = sub_sam, aes(ETIME,
                                                 predict(nl_model)), inherit.aes = FALSE, color = "red",
                             linewidth = 1)
    }
    all_fig[[mesu]] <- ggplotGrob(fig)
    grid::grid.draw(all_fig[[mesu]])
    all_plot[[mesu]] <- recordPlot()
    names(all_plot)[length(all_plot)] <- measurement_id
  }
  return(all_plot)
}
