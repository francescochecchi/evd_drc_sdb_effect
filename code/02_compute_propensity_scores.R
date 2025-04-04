#...............................................................................
### +++ EFFECT OF SAFE AND DIGNIFIED BURIALS ON EBOV TRANSMISSION IN DRC +++ ###
#...............................................................................

#...............................................................................
## ---- R SCRIPT TO COMPUTE (GENERALISED OR BINOMIAL) PROPENSITY SCORES ----- ##
#...............................................................................

                            # Written by Francesco Checchi, LSHTM (March 2023)
                            # francesco.checchi@lshtm.ac.uk 



#...............................................................................
### Exploring univariate associations and candidate lags / time windows
#................................................................................

  #...................................
  ## Read pre-processed dataset
  tsw <- read.csv(paste(dir_output, "/out_tsw.csv", sep = "") )  


  #...................................
  ## Observe graphically the correlation between different indicators and 
    # transmissibility
    
    # Continuous variables 
    x <- c("p_success_prev", "p_timely_prev", "r_sdb_to_pop_prev")
    for (i in x) {
      plot <- ggplot(data = tsw, mapping = aes(tsw[, i], p_curr)) +
        geom_point(colour = palette_cb[4], fill = palette_cb[4], alpha = 0.5) +
        scale_y_continuous("proportion of cases in current period", 
          limits = c(0, 1)) +
        scale_x_continuous(i) +
        theme_bw() +
        geom_smooth(method = "loess")
      print(plot)
    }

    # # Categorical variables
    # x <- paste0(c("p_success_prev", "p_timely_prev","r_sdb_to_pop_prev"),"_cat")
    # for (i in x) {
    #   plot <- ggplot(data = tsw, mapping = aes(tsw[, i], p_curr)) +
    #     geom_boxplot(colour = palette_cb[4], fill = palette_cb[4], alpha = 0.5)+
    #     scale_y_continuous("proportion of cases in current period", 
    #       limits = c(0, 1)) +
    #     scale_x_discrete(i) +
    #     theme_bw()
    #   print(plot)
    # }
    # 
    
  #...................................
  ## Compare mixed and simple generalised linear models for SDB performance on 
    # relative period incidence
  x <- c("p_success_prev", "p_timely_prev", "r_sdb_to_pop_prev")
  
  for (i in x) {    
      print("//////////////////////")

      # fixed effects only model
      print("---Fixed effects ---------------------")
      fit_glm <- f_glm(i, tsw, F, window_transmission)
      print( model_parameters(fit_glm, exponentiate = T) )
      print(paste("AIC:", AIC(fit_glm)) )
            
      print("////")
      
      # mixed effects model
      print("----Random effects ---------------------")
      fit_glmm <- f_glmm(i, tsw, F, window_transmission)
      print( model_parameters(fit_glmm, exponentiate = T) )
      print(paste("AIC:", AIC(fit_glmm)) )
  }
    
  # retain mixed model (seems more appropriate)

  
  #...................................
  ## Select confounder time windows and lags with best fit for p success prev
  for (i in c("vcp_", "n_against_evd_", "n_suspensions_", "rate_events_", 
    "rate_fatalities_", "etc_open_") ) { 
    for (j in windows_other) {
      for (k in lags) {
        print("--------------------------------")
        print(" ")
        form <- as.formula(paste0("p_success_prev ~ ", 
          paste0(i, j, "w_lag", k)))
        fit <- lm(form, tsw)
        print(model_parameters(fit, exponentiate = F) )
        print(paste("AIC:", AIC(fit)) )
      }  
    }
  }
  
    confounders_inc_p_success <- c("p_nosocomial_prev", "p_epi_link_prev", 
      "vcp_4w_lag2", "etc_open_4w_lag2", "n_against_evd_4w_lag2", 
      "rate_events_4w_lag2", "n_suspensions_4w_lag2", 
      "n_networks", "n_frequencies", "rate_hf", "rate_road_length")    
  
  
  #...................................
  ## Select confounder time windows and lags with best fit for p success curr
  for (i in c("vcp_", "n_against_evd_", "n_suspensions_", "rate_events_", 
    "rate_fatalities_", "etc_open_") ) { 
    for (j in windows_other) {
      for (k in lags) {
        print("--------------------------------")
        print(" ")
        form <- as.formula(paste0("p_success_curr ~ ", 
          paste0(i, j, "w_lag", k)))
        fit <- lm(form, tsw)
        print(model_parameters(fit, exponentiate = F) )
        print(paste("AIC:", AIC(fit)) )
      }  
    }
  }
  
    confounders_rn_p_success <- c("p_nosocomial_curr", "p_epi_link_curr", 
      "vcp_1w_lag2", "etc_open_4w_lag2", "n_against_evd_4w_lag2", 
      "rate_events_4w_lag2", "n_suspensions_4w_lag2", 
      "n_networks", "n_frequencies", "rate_hf", "rate_road_length") 
  
  
  #...................................
  ## Select confounder time windows and lags with best fit for p timely prev
  for (i in c("vcp_", "n_against_evd_", "n_suspensions_", "rate_events_", 
    "rate_fatalities_", "etc_open_") ) { 
    for (j in windows_other) {
      for (k in lags) {
        print("--------------------------------")
        print(" ")
        form <- as.formula(paste0("p_timely_prev ~ ", 
          paste0(i, j, "w_lag", k)))
        fit <- lm(form, tsw)
        print(model_parameters(fit, exponentiate = F) )
        print(paste("AIC:", AIC(fit)) )
      }  
    }
  }
  
    confounders_inc_p_timely <- c("p_nosocomial_prev", "p_epi_link_prev", 
      "vcp_1w_lag2", "etc_open_4w_lag2", "n_against_evd_4w_lag2", 
      "rate_events_4w_lag2", "n_suspensions_4w_lag2", 
      "n_networks", "n_frequencies", "rate_hf", "rate_road_length")

  #...................................
  ## Select confounder time windows and lags with best fit for p timely curr
  for (i in c("vcp_", "n_against_evd_", "n_suspensions_", "rate_events_", 
    "rate_fatalities_", "etc_open_") ) { 
    for (j in windows_other) {
      for (k in lags) {
        print("--------------------------------")
        print(" ")
        form <- as.formula(paste0("p_timely_prev ~ ", 
          paste0(i, j, "w_lag", k)))
        fit <- lm(form, tsw)
        print(model_parameters(fit, exponentiate = F) )
        print(paste("AIC:", AIC(fit)) )
      }  
    }
  }
  
    confounders_rn_p_timely <- c("p_nosocomial_curr", "p_epi_link_curr", 
      "vcp_1w_lag2", "etc_open_4w_lag2", "n_against_evd_4w_lag2", 
      "rate_events_4w_lag2", "n_suspensions_4w_lag2", 
      "n_networks", "n_frequencies", "rate_hf", "rate_road_length")


  #...................................
  ## Select confounder time windows and lags with best fit for r sdb to pop prev
  for (i in c("vcp_", "n_against_evd_", "n_suspensions_", "rate_events_", 
    "rate_fatalities_", "etc_open_") ) { 
    for (j in windows_other) {
      for (k in lags) {
        print("--------------------------------")
        print(" ")
        form <- as.formula(paste0("r_sdb_to_pop_prev ~ ", 
          paste0(i, j, "w_lag", k)))
        fit <- lm(form, tsw)
        print(model_parameters(fit, exponentiate = F) )
        print(paste("AIC:", AIC(fit)) )
      }  
    }
  }
  
    confounders_inc_r_sdb_to_pop <- c("p_nosocomial_prev", "p_epi_link_prev", 
      "vcp_1w_lag0", "etc_open_4w_lag2", "n_against_evd_4w_lag2", 
      "rate_events_4w_lag2", "n_suspensions_4w_lag2", 
      "n_networks", "n_frequencies", "rate_hf", "rate_road_length")    
    
  #...................................
  ## Select confounder time windows and lags with best fit for r sdb to pop curr
  for (i in c("vcp_", "n_against_evd_", "n_suspensions_", "rate_events_", 
    "rate_fatalities_", "etc_open_") ) { 
    for (j in windows_other) {
      for (k in lags) {
        print("--------------------------------")
        print(" ")
        form <- as.formula(paste0("r_sdb_to_pop_curr ~ ", 
          paste0(i, j, "w_lag", k)))
        fit <- lm(form, tsw)
        print(model_parameters(fit, exponentiate = F) )
        print(paste("AIC:", AIC(fit)) )
      }  
    }
  }
  
    confounders_rn_r_sdb_to_pop <- c("p_nosocomial_curr", "p_epi_link_curr", 
      "vcp_1w_lag0", "etc_open_4w_lag2", "n_against_evd_4w_lag2", 
      "rate_events_4w_lag2", "n_suspensions_4w_lag2", 
      "n_networks", "n_frequencies", "rate_hf", "rate_road_length")    
    
    
        
#...............................................................................
### Preparing datasets for analysis
#...............................................................................

  #...................................
  ## Specify outcomes and exposures
    # Outcomes
    outcomes <- c("cases_curr", "cases_prev", "ln_r_mean", "r_mean")
  
    # Exposures 
    exposures <- c("p_success_prev", "p_timely_prev_cat", 
      "r_sdb_to_pop_prev", "p_success_curr", "p_timely_curr_cat",
      "r_sdb_to_pop_curr")
    
  #...................................
  ## Retain only observations with non-missing outcome, exposure and confounders
    # For incidence change
    x <- c(outcomes[1:2], exposures[1], confounders_inc_p_success, 
      "hz", "epi_year", "epi_week")
    tsw_inc_p_success <- tsw[complete.cases(tsw[, x]), ]
    nrow(tsw_inc_p_success)
    
    x <- c(outcomes[1:2], exposures[2], confounders_inc_p_timely, 
      "hz", "epi_year", "epi_week")
    tsw_inc_p_timely <- tsw[complete.cases(tsw[, x]), ]
    nrow(tsw_inc_p_timely)
    
    x <- c(outcomes[1:2], exposures[3], confounders_inc_r_sdb_to_pop, 
      "hz", "epi_year", "epi_week")
    tsw_inc_r_sdb_to_pop <- tsw[complete.cases(tsw[, x]), ]
    nrow(tsw_inc_r_sdb_to_pop)

    # For reproduction number
    x <- c(outcomes[3], exposures[4], confounders_rn_p_success, 
      "hz", "epi_year", "epi_week")
    tsw_rn_p_success <- tsw[complete.cases(tsw[, x]), ]
    nrow(tsw_rn_p_success)     

    x <- c(outcomes[3], exposures[5], confounders_rn_p_timely, 
      "hz", "epi_year", "epi_week")
    tsw_rn_p_timely <- tsw[complete.cases(tsw[, x]), ]
    nrow(tsw_rn_p_timely)     
    
    x <- c(outcomes[3], exposures[6], confounders_rn_r_sdb_to_pop, 
      "hz", "epi_year", "epi_week")
    tsw_rn_r_sdb_to_pop <- tsw[complete.cases(tsw[, x]), ]
    nrow(tsw_rn_r_sdb_to_pop)     
    
    
  #...................................
  ## Transform exposures and confounders
    # Exposures
      # proportion of successful SDB
      f_hist(exposures[1], tsw_inc_p_success, c(NA, NA))
      f_hist(exposures[4], tsw_rn_p_success, c(NA, NA))
        # no need to transform

      # proportion of timely SDB (binary)
      tsw_inc_p_timely[, exposures[2]] <- 
        as.character(tsw_inc_p_timely[, exposures[2]])
      tsw_inc_p_timely[which(tsw_inc_p_timely[, exposures[2]] == "< 1.00"),
        exposures[2]] <- 0
      tsw_inc_p_timely[which(tsw_inc_p_timely[, exposures[2]] == "1.00"),
        exposures[2]] <- 1
      tsw_inc_p_timely[, exposures[2]] <- as.integer(tsw_inc_p_timely[,
        exposures[2]])
      print(table(tsw_inc_p_timely[, exposures[2]]))

      tsw_rn_p_timely[, exposures[5]] <- 
        as.character(tsw_rn_p_timely[, exposures[5]])      
      tsw_rn_p_timely[which(tsw_rn_p_timely[, exposures[5]] == "< 1.00"),
        exposures[5]] <- 0
      tsw_rn_p_timely[which(tsw_rn_p_timely[, exposures[5]] == "1.00"),
        exposures[5]] <- 1
      tsw_rn_p_timely[, exposures[5]] <- as.integer(tsw_rn_p_timely[,
        exposures[5]])
      print(table(tsw_rn_p_timely[, exposures[5]]))
        # no need to transform
      
      # rate of SDBs to population
      f_hist(exposures[3], tsw_inc_r_sdb_to_pop, c(NA, NA))
      f_hist(exposures[6], tsw_rn_r_sdb_to_pop, c(NA, NA))
        # no need to transform

    # Scale and center all continuous confounders
    tsw_inc_p_success[, confounders_inc_p_success] <-
      scale(tsw_inc_p_success[, confounders_inc_p_success])
    tsw_rn_p_success[, confounders_rn_p_success] <-
      scale(tsw_rn_p_success[, confounders_rn_p_success])
    tsw_inc_p_timely[, confounders_inc_p_timely] <-
      scale(tsw_inc_p_timely[, confounders_inc_p_timely])
    tsw_rn_p_timely[, confounders_rn_p_timely] <-
      scale(tsw_rn_p_timely[, confounders_rn_p_timely])
    tsw_inc_r_sdb_to_pop[, confounders_inc_r_sdb_to_pop] <-
      scale(tsw_inc_r_sdb_to_pop[, confounders_inc_r_sdb_to_pop])
    tsw_rn_r_sdb_to_pop[, confounders_rn_r_sdb_to_pop] <-
      scale(tsw_rn_r_sdb_to_pop[, confounders_rn_r_sdb_to_pop])


#...............................................................................
### Computing propensity scores and resulting confounder-exposure balance
#................................................................................
   
  #...................................
  ## Proportion of successful SDBs (continuous)
    
    # Incidence change
    out <- f_gps(tsw_inc_p_success, exposures[1], confounders_inc_p_success, 
      palette_cb[5])
    x <- out$dataset
    colnames(x) <- gsub("gps", "gps_p_success", colnames(x))
    tsw_inc_p_success <- merge(
      tsw_inc_p_success, x[, c("hz", "epi_year", "epi_week", 
      grep("gps_p_success", colnames(x), value = T))],
      by = c("hz", "epi_year", "epi_week"))
    write.csv(out$balance_stats, 
      paste(dir_output, "/out_ps_bal_inc_p_success.csv", sep = ""),
      row.names = F)
    out$balance_plot
    ggsave(paste(dir_output, "/out_ps_bal_inc_p_success.png", sep = ""),
      dpi = "print", width = 15, height = 10, units = "cm")
      
      # plot correlation between GPS and exposure
      df <- tsw_inc_p_success
      df$exposure_cat <- cut(df[, exposures[1]], breaks = seq(0, 1, 0.2),
        include.lowest = T, right = F)
      ggplot(df, aes(y = gps_p_success, x = exposure_cat)) +
        geom_boxplot(alpha = 0.5, colour = palette_cb[5], 
          fill = palette_cb[5]) +
        theme_bw() +
        scale_x_discrete("percentage of successful SDBs (categorised)") +
        scale_y_continuous("distribution of generalised propensity scores")
      ggsave(paste(dir_output, "/out_ps_dist_inc_p_success.png", sep = ""),
        dpi = "print", width = 22, height = 13, units = "cm")
        
    # Reproduction number
    out <- f_gps(tsw_rn_p_success, exposures[4], confounders_rn_p_success, 
      palette_cb[11])
    x <- out$dataset    
    colnames(x) <- gsub("gps", "gps_p_success", colnames(x))
    tsw_rn_p_success <- merge(
      tsw_rn_p_success, x[, c("hz", "epi_year", "epi_week", 
      grep("gps_p_success", colnames(x), value = T))],
      by = c("hz", "epi_year", "epi_week"))
    write.csv(out$balance_stats, 
      paste(dir_output, "/out_ps_bal_rn_p_success.csv", sep = ""),
      row.names = F)
    out$balance_plot
    ggsave(paste(dir_output, "/out_ps_bal_rn_p_success.png", sep = ""),
      dpi = "print", width = 15, height = 10, units = "cm")

      # plot correlation between GPS and exposure
      df <- tsw_rn_p_success
      df$exposure_cat <- cut(df[, exposures[4]], breaks = seq(0, 1, 0.2),
        include.lowest = T, right = F)
      ggplot(df, aes(y = gps_p_success, x = exposure_cat)) +
        geom_boxplot(alpha = 0.5, colour = palette_cb[11], 
          fill = palette_cb[11]) +
        theme_bw() +
        scale_x_discrete("percentage of successful SDBs (categorised)") +
        scale_y_continuous("distribution of generalised propensity scores")
      ggsave(paste(dir_output, "/out_ps_dist_rn_p_success.png", sep = ""),
        dpi = "print", width = 22, height = 13, units = "cm")

    
  #...................................
  ## Proportion of timely SDBs (binary)
    
    # Incidence change
      # define formula for exposure as a function of confounders
      formula_exp <- formula(paste(exposures[2]," ~ ", 
        paste(confounders_inc_p_timely, collapse = " + ")))
    
      # random forest model of exposure as a function of confounders
      fit_exp <- MatchIt::matchit(formula_exp, data = tsw_inc_p_timely, 
        method = "nearest", ratio = 2, caliper = 0.2, distance = "glm")
      
      # extract PS values
      tsw_inc_p_timely[, "ps_p_timely"] <- fit_exp$distance
        # distribution of PS
        f_hist("ps_p_timely", tsw_inc_p_timely, c(NA, NA))

      # inspect covariate balance
      out <- cobalt::bal.tab(fit_exp, m.threshold = 0.10, un = T)
      out
      write.csv(out$Balance, paste(dir_output, "/out_ps_bal_inc_p_timely.csv", 
        sep = ""), row.names = T)
      cobalt::love.plot(fit_exp, binary = "std", continuous = "std", 
        m.threshold = 0.10, colors = palette_cb[c(12, 4)])
      ggsave(paste(dir_output, "/out_ps_bal_inc_p_timely.png", sep = ""), 
        dpi = "print", width = 20, height = 10, units = "cm")

      # extract matched dataset
      tsw_inc_p_timely <- MatchIt::match.data(fit_exp)
      
      # plot distributions of PS by exposure
      df <- tsw_inc_p_timely
      df$exposure <- factor(df[, exposures[2]],
        levels = c("0", "1"), labels = c("<100%", "100%"))
      
      ggplot(df, aes(x = subclass, y = ps_p_timely, colour = exposure,
        fill = exposure)) +
        geom_point(alpha = 0.5, size = 3) +
        theme_bw() +
        scale_y_continuous("propensity score") +
        scale_x_discrete("matched set") +
        scale_colour_manual("proportion of timely SDBs", 
          values = palette_cb[c(4,12)]) +
        scale_fill_manual("proportion of timely SDBs", 
          values = palette_cb[c(4,12)]) +
        theme(legend.position = "top")
      ggsave(paste(dir_output, "/out_ps_paired_inc_p_timely.png", sep = ""), 
        dpi = "print", width = 22, height = 13, units = "cm")
      
      ggplot(df, aes(x = ps_p_timely, colour = exposure,
        fill = exposure)) +
        geom_histogram(alpha = 0.5) +
        theme_bw() +
        scale_y_continuous("number of matched observations", 
          breaks = seq(0, 10, 2)) +
        scale_x_continuous("propensity score") +
        scale_colour_manual("proportion of timely SDBs", 
          values = palette_cb[c(4,12)]) +
        scale_fill_manual("proportion of timely SDBs", 
          values = palette_cb[c(4,12)]) +
        facet_wrap(exposure~., ncol = 1, nrow = 2) +
        theme(legend.position = "top")
      ggsave(paste(dir_output, "/out_ps_dist_inc_p_timely.png", sep = ""), 
        dpi = "print", width = 22, height = 13, units = "cm")
      
    # Reproduction number
      # define formula for exposure as a function of confounders
      formula_exp <- formula(paste(exposures[5]," ~ ", 
        paste(confounders_rn_p_timely, collapse = " + ")))
    
      # random forest model of exposure as a function of confounders
      fit_exp <- MatchIt::matchit(formula_exp, data = tsw_rn_p_timely, 
        method = "nearest", ratio = 2, caliper = 0.2, distance = "randomforest")
      
      # extract PS values
      tsw_rn_p_timely$ps_p_timely <- fit_exp$distance
        # distribution of PS
        f_hist("ps_p_timely", tsw_rn_p_timely, c(NA, NA))

      # inspect covariate balance
      out <- cobalt::bal.tab(fit_exp, m.threshold = 0.10, un = T)
      out
      write.csv(out$Balance, paste(dir_output, "/out_ps_bal_rn_p_timely.csv", 
        sep = ""), row.names = T)
      cobalt::love.plot(fit_exp, binary = "std", continuous = "std", 
        m.threshold = 0.10, colors = palette_cb[c(12, 4)])
      ggsave(paste(dir_output, "/out_ps_bal_rn_p_timely.png", sep = ""), 
        dpi = "print", width = 20, height = 10, units = "cm")
    
      # extract matched dataset
      tsw_rn_p_timely <- MatchIt::match.data(fit_exp)     

      # plot distributions of PS by exposure
      df <- tsw_rn_p_timely
      df$exposure <- factor(df[, exposures[5]],
        levels = c("0", "1"), labels = c("<100%", "100%"))
      
      ggplot(df, aes(x = subclass, y = ps_p_timely, colour = exposure,
        fill = exposure)) +
        geom_point(alpha = 0.5, size = 3) +
        theme_bw() +
        scale_y_continuous("propensity score") +
        scale_x_discrete("matched set") +
        scale_colour_manual("proportion of timely SDBs", 
          values = palette_cb[c(4,12)]) +
        scale_fill_manual("proportion of timely SDBs", 
          values = palette_cb[c(4,12)]) +
        theme(legend.position = "top")
      ggsave(paste(dir_output, "/out_ps_paired_rn_p_timely.png", sep = ""), 
        dpi = "print", width = 22, height = 13, units = "cm")
      
      ggplot(df, aes(x = ps_p_timely, colour = exposure,
        fill = exposure)) +
        geom_histogram(alpha = 0.5) +
        theme_bw() +
        scale_y_continuous("number of matched observations", 
          breaks = seq(0, 10, 2)) +
        scale_x_continuous("propensity score") +
        scale_colour_manual("proportion of timely SDBs", 
          values = palette_cb[c(4,12)]) +
        scale_fill_manual("proportion of timely SDBs", 
          values = palette_cb[c(4,12)]) +
        facet_wrap(exposure~., ncol = 1, nrow = 2) +
        theme(legend.position = "top")
      ggsave(paste(dir_output, "/out_ps_dist_rn_p_timely.png", sep = ""), 
        dpi = "print", width = 22, height = 13, units = "cm")
       

  #...................................
  ## Rate of SDBs per population (continuous)

    # Incidence change
    out <- f_gps(tsw_inc_r_sdb_to_pop,exposures[3],confounders_inc_r_sdb_to_pop, 
      palette_cb[13])
    x <- out$dataset
    colnames(x) <- gsub("gps", "gps_r_sdb_to_pop", colnames(x))
    tsw_inc_r_sdb_to_pop <- merge(
      tsw_inc_r_sdb_to_pop, x[, c("hz", "epi_year", "epi_week", 
      grep("gps_r_sdb_to_pop", colnames(x), value = T))],
      by = c("hz", "epi_year", "epi_week"))
    write.csv(out$balance_stats, 
      paste(dir_output, "/out_ps_bal_inc_r_sdb_to_pop.csv", sep = ""),
      row.names = F)
    out$balance_plot
    ggsave(paste(dir_output, "/out_ps_bal_inc_r_sdb_to_pop.png", sep = ""),
      dpi = "print", width = 15, height = 10, units = "cm")
      
      # plot correlation between GPS and exposure
      df <- tsw_inc_r_sdb_to_pop
      df$exposure_cat <- cut(df[, exposures[3]], breaks = seq(0, 18, 3),
        include.lowest = T, right = F)
      ggplot(df, aes(y = gps_r_sdb_to_pop, x = exposure_cat)) +
        geom_boxplot(alpha = 0.5, colour = palette_cb[5], 
          fill = palette_cb[5]) +
        theme_bw() +
        scale_x_discrete("SDBs per 100,000 person-weeks (categorised)") +
        scale_y_continuous("distribution of generalised propensity scores")
      ggsave(paste(dir_output, "/out_ps_dist_inc_r_sdb_to_pop.png", sep = ""),
        dpi = "print", width = 22, height = 13, units = "cm")
        
    # Reproduction number
    out <- f_gps(tsw_rn_r_sdb_to_pop, exposures[6], confounders_rn_r_sdb_to_pop, 
      palette_cb[13])
    x <- out$dataset    
    colnames(x) <- gsub("gps", "gps_r_sdb_to_pop", colnames(x))
    tsw_rn_r_sdb_to_pop <- merge(
      tsw_rn_r_sdb_to_pop, x[, c("hz", "epi_year", "epi_week", 
      grep("gps_r_sdb_to_pop", colnames(x), value = T))],
      by = c("hz", "epi_year", "epi_week"))
    write.csv(out$balance_stats, 
      paste(dir_output, "/out_ps_bal_rn_r_sdb_to_pop.csv", sep = ""),
      row.names = F)
    out$balance_plot
    ggsave(paste(dir_output, "/out_ps_bal_rn_r_sdb_to_pop.png", sep = ""),
      dpi = "print", width = 15, height = 10, units = "cm")

      # plot correlation between GPS and exposure
      df <- tsw_rn_r_sdb_to_pop
      df$exposure_cat <- cut(df[, exposures[6]], breaks = seq(0, 18, 3),
        include.lowest = T, right = F)
      ggplot(df, aes(y = gps_r_sdb_to_pop, x = exposure_cat)) +
        geom_boxplot(alpha = 0.5, colour = palette_cb[11], 
          fill = palette_cb[11]) +
        theme_bw() +
        scale_x_discrete("SDBs per 100,000 person-weeks (categorised)") +
        scale_y_continuous("distribution of generalised propensity scores")
      ggsave(paste(dir_output, "/out_ps_dist_rn_r_sdb_to_pop.png", sep = ""),
        dpi = "print", width = 22, height = 13, units = "cm")



#.........................................................................................
### ENDS
#.........................................................................................
  
