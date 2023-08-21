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
    # transmissibility / continuous variables 

    # 1-week time windows
    x <- c("p_success_prev_1w", "p_timely_prev_1w", "p_both_ok_prev_1w", 
      "r_sdb_to_pop_prev_1w")
    for (i in x) {
      plot <- ggplot(data = tsw, mapping = aes(tsw[, i], p_curr_1w)) +
        geom_point(colour = palette_cb[4], fill = palette_cb[4], alpha = 0.5) +
        scale_y_continuous("proportion of cases in current period", 
          limits = c(0, 1)) +
        scale_x_continuous(i) +
        theme_bw() +
        geom_smooth(method = "loess")
      print(plot)
    }
    
    # 3-week time windows  
    x <- c("p_success_prev_3w", "p_timely_prev_3w", "p_both_ok_prev_3w", 
      "r_sdb_to_pop_prev_3w")
    for (i in x) {
      plot <- ggplot(data = tsw, mapping = aes(tsw[, i], p_curr_3w)) +
        geom_point(colour = palette_cb[4], fill = palette_cb[4], alpha = 0.5) +
        scale_y_continuous("proportion of cases in current period", 
          limits = c(0, 1)) +
        scale_x_continuous(i) +
        theme_bw() +
        geom_smooth(method = "loess")
      print(plot)
    }

  #...................................
  ## Observe graphically the correlation between different indicators and 
    # transmissibility / categorical variables 
      
    # 1-week time windows
    x <- paste(c("p_success_prev_1w", "p_timely_prev_1w", "p_both_ok_prev_1w", 
      "r_sdb_to_pop_prev_1w"), "_cat", sep="")
    for (i in x) {
      plot <- ggplot(data = tsw, mapping = aes(tsw[, i], p_curr_1w)) +
        geom_boxplot(colour = palette_cb[4], fill = palette_cb[4], alpha = 0.5)+
        scale_y_continuous("proportion of cases in current period", 
          limits = c(0, 1)) +
        scale_x_discrete(i) +
        theme_bw()
      print(plot)
    }

    # 3-week time windows
    x <- paste(c("p_success_prev_3w", "p_timely_prev_3w", "p_both_ok_prev_3w", 
      "r_sdb_to_pop_prev_3w"), "_cat", sep="")
    for (i in x) {
      plot <- ggplot(data = tsw, mapping = aes(tsw[, i], p_curr_3w)) +
        geom_boxplot(colour = palette_cb[4], fill = palette_cb[4], alpha = 0.5)+
        scale_y_continuous("proportion of cases in current period", 
          limits = c(0, 1)) +
        scale_x_discrete(i) +
        theme_bw()
      print(plot)
    }
    
    
  #...................................
  ## Compare mixed and simple generalised linear models for SDB performance on 
    # relative period incidence
  x <- c("p_success_prev_", "p_timely_prev_", "p_both_ok_prev_", 
    "r_sdb_to_pop_prev_")
  
  for (j in x) {    
    for (i in windows_transmission) {
      print("//////////////////////")

      # fixed effects only model
      print("---Fixed effects ---------------------")
      fit_glm <- f_glm(paste(j, i, "w_cat", sep=""), tsw, FALSE, i)
      print( model_parameters(fit_glm, exponentiate = TRUE) )
      print(paste("AIC:", AIC(fit_glm)) )
            
      print("////")
      
      # mixed effects model
      print("----Random effects ---------------------")
      fit_glmm <- f_glmm(paste(j, i, "w_cat", sep=""), tsw, FALSE, i)
      print( model_parameters(fit_glmm, exponentiate = TRUE) )
      print(paste("AIC:", AIC(fit_glmm)) )
    }
  }
    
  # retain mixed model (better AIC, seems more appropriate), 
  # choose 3 weeks window (clearer effects),
  # keep all SDBs (not just community or success), 
  # keep ratio of SDBs (all) to population (best measure of coverage 
    # among those available)

  
  #...................................
  ## Select confounder time windows and lags with best fit / plausibility   
  for (i in windows_other) {
    for (j in lags) {
      for (k in c("vc_", "n_against_evd_", "n_suspensions_", "rate_events_", 
        "rate_fatalities_", "etc_open_") ) { 
        print("/")
        print("--------------------------------")
        fit_glmm <- f_glmm(paste(k, i, "w_lag", j, "_cat", sep=""), tsw,FALSE,3)
        print( model_parameters(fit_glmm, exponentiate = TRUE) )
        print(paste("AIC:", AIC(fit_glmm)) )
      }  
    }
  }
      
    # retain 3w instead of 5w (generally much better fit) and lag of 2 weeks 
      # (better fit and more plausible: length of 1 serial interval)
    # retain rate of insecurity events instead of insecurity deaths (better fit)

  
#...............................................................................
### Preparing datasets for analysis
#...............................................................................

  #...................................
  ## Specify outcomes, exposures and confounders
    # Outcomes
    outcomes <- c("cases_curr_3w", "cases_prev_3w", "ln_r_mean", "r_mean")
  
    # Exposures 
    exposures <- c("p_success_prev_3w", "p_timely_prev_3w_cat", 
      "p_both_ok_prev_3w", "p_success_curr_3w", "p_timely_curr_3w_cat", 
      "p_both_ok_curr_3w")
    
    # List confounders (static and time-varying - all in their continuous form)
      # For incidence change outcome
      confounders_inc <- c("p_nosocomial_prev_3w", "p_epi_link_prev_3w", 
        "vc_3w_lag2", "etc_open_3w_lag2", "n_against_evd_3w_lag2", 
        "rate_events_3w_lag2", "n_networks", "n_frequencies", "rate_hf", 
        "rate_road_length", "r_sdb_to_pop_prev_3w")    
        # "n_suspensions_3w_lag2" excluded because it appears to be 
        # highly collinear (rank deficient models)

      # For reproduction number outcome
      confounders_rn <- c("p_nosocomial_curr_3w", "p_epi_link_curr_3w", 
        "vc_3w_lag2", "etc_open_3w_lag2", "n_against_evd_3w_lag2", 
        "rate_events_3w_lag2", "n_networks", "n_frequencies", "rate_hf", 
        "rate_road_length", "r_sdb_to_pop_curr_3w")    
        # "n_suspensions_3w_lag2" excluded because it appears to be 
        # highly collinear (rank deficient models)
        # stick with lag2 because model fits better than with lag0
       
  #...................................
  ## Retain only observations with non-missing outcome, exposure and confounders
    # For incidence change
    tsw_inc <- tsw[complete.cases(tsw[, c(outcomes[1:2], exposures[1:3], 
      confounders_inc)]), ]
    nrow(tsw_inc)     

    # For reproduction number
    tsw_rn <- tsw[complete.cases(tsw[, c(outcomes[4], exposures[4:6], 
      confounders_rn)]), ]
    nrow(tsw_rn)     

    
  #...................................
  ## Transform exposures and confounders
    # Exposures
      # proportion of successful SDB
      f_hist(exposures[1], tsw_inc, c(NA, NA))
      f_hist(exposures[4], tsw_inc, c(NA, NA))
      f_hist(exposures[1], tsw_rn, c(NA, NA))
      f_hist(exposures[4], tsw_rn, c(NA, NA))
        # no need to transform

      # proportion of timely SDB (binary)
      for (i in c(2, 5)) {
        tsw_inc[which(tsw_inc[, exposures[i]] == "< 1.00"), exposures[i]] <- 0
        tsw_inc[which(tsw_inc[, exposures[i]] == "1.00"), exposures[i]] <- 1
        tsw_inc[, exposures[i]] <- as.integer(tsw_inc[, exposures[i]])
        print(table(tsw_inc[, exposures[i]]))
        
        tsw_rn[which(tsw_rn[, exposures[i]] == "< 1.00"), exposures[i]] <- 0
        tsw_rn[which(tsw_rn[, exposures[i]] == "1.00"), exposures[i]] <- 1
        tsw_rn[, exposures[i]] <- as.integer(tsw_rn[, exposures[i]])
        print(table(tsw_rn[, exposures[i]]))
      }
  
      # proportion of successful and timely SDB
      f_hist(exposures[3], tsw_inc, c(NA, NA))
      f_hist(exposures[6], tsw_inc, c(NA, NA))
      f_hist(exposures[3], tsw_rn, c(NA, NA))
      f_hist(exposures[6], tsw_rn, c(NA, NA))
        # no need to transform

    # Scale and center all confounders
    tsw_inc[, confounders_inc] <- scale(tsw_inc[, confounders_inc])
    tsw_rn[, confounders_rn] <- scale(tsw_rn[, confounders_rn])
    

#...............................................................................
### Computing propensity scores and resulting confounder-exposure balance
#................................................................................
   
  #...................................
  ## Proportion of successful SDBs (continuous)
    
    # Incidence change
    out <- f_gps(tsw_inc, exposures[1], confounders_inc, palette_cb[6])
    x <- out$dataset
    colnames(x) <- gsub("gps", "gps_p_success", colnames(x))
    tsw_inc <- merge(tsw_inc, x[, c("hz", "epi_year", "epi_week", 
      grep("gps_p_success", colnames(x), value = TRUE))],
      by = c("hz", "epi_year", "epi_week"))
    write.csv(out$balance_stats, 
      paste(dir_output, "/out_ps_bal_inc_p_success.csv", sep = ""),
      row.names = FALSE)
    plot <- out$balance_plot
    ggsave(paste(dir_output, "/out_ps_bal_inc_p_success.png", sep = ""),
      dpi = "print", width = 15, height = 10, units = "cm")
              
    # Reproduction number
    out <- f_gps(tsw_rn, exposures[4], confounders_rn, palette_cb[4])
    x <- out$dataset    
    colnames(x) <- gsub("gps", "gps_p_success", colnames(x))
    tsw_rn <- merge(tsw_rn, x[, c("hz", "epi_year", "epi_week", 
      grep("gps_p_success", colnames(x), value = TRUE))],
      by = c("hz", "epi_year", "epi_week"))
    write.csv(out$balance_stats, 
      paste(dir_output, "/out_ps_bal_rn_p_success.csv", sep = ""),
      row.names = FALSE)
    plot <- out$balance_plot
    ggsave(paste(dir_output, "/out_ps_bal_rn_p_success.png", sep = ""),
      dpi = "print", width = 15, height = 10, units = "cm")

    
  # #...................................
  # ## Proportion of successful AND timely SDBs (continuous)
  #     
  #   # Incidence change
  #   out <- f_gps(tsw_inc, exposures[3], confounders_inc, palette_cb[6])
  #   x <- out$dataset
  #   colnames(x) <- gsub("gps", "gps_p_both", colnames(x))
  #   tsw_inc <- merge(tsw_inc, x[, c("hz", "epi_year", "epi_week", 
  #     grep("gps_p_both", colnames(x), value = TRUE))],
  #     by = c("hz", "epi_year", "epi_week"))
  #   write.csv(out$balance_stats, 
  #     paste(dir_output, "/out_ps_bal_inc_p_both.csv", sep = ""),
  #     row.names = FALSE)
  #   plot <- out$balance_plot
  #   ggsave(paste(dir_output, "/out_ps_bal_inc_p_both.png", sep = ""),
  #     dpi = "print", width = 15, height = 10, units = "cm")
  #             
  #   # Reproduction number
  #   out <- f_gps(tsw_rn, exposures[6], confounders_rn, palette_cb[4])
  #   x <- out$dataset    
  #   colnames(x) <- gsub("gps", "gps_p_both", colnames(x))
  #   tsw_rn <- merge(tsw_rn, x[, c("hz", "epi_year", "epi_week", 
  #     grep("gps_p_both", colnames(x), value = TRUE))],
  #     by = c("hz", "epi_year", "epi_week"))
  #   write.csv(out$balance_stats, 
  #     paste(dir_output, "/out_ps_bal_rn_p_both.csv", sep = ""),
  #     row.names = FALSE)
  #   plot <- out$balance_plot
  #   ggsave(paste(dir_output, "/out_ps_bal_rn_p_both.png", sep = ""),
  #     dpi = "print", width = 15, height = 10, units = "cm")
  # 
    
  #...................................
  ## Proportion of timely SDBs (binary)
    
    # Incidence change
      # define formula for exposure as a function of confounders
      formula_exp <- formula(paste(exposures[2]," ~ ", 
        paste(confounders_inc, collapse = " + ")))
    
      # random forest model of exposure as a function of confounders
      fit_exp <- MatchIt::matchit(formula_exp, data = tsw_inc, 
        method = "nearest", ratio = 2, caliper = 0.2, distance = "randomforest")
      
      # extract PS values
      tsw_inc[, "ps_p_timely"] <- fit_exp$distance
        # distribution of PS
        f_hist("ps_p_timely", tsw_inc, c(NA, NA))

      # inspect covariate balance
      out <- cobalt::bal.tab(fit_exp, m.threshold = 0.10, un = TRUE)
      out
      write.csv(out$Balance, paste(dir_output, "/out_ps_bal_inc_p_timely.csv", 
        sep = ""), row.names = TRUE)
      cobalt::love.plot(fit_exp, binary = "std", continuous = "std", 
        m.threshold = 0.10, colors = palette_cb[c(7, 6)])
      ggsave(paste(dir_output, "/out_ps_bal_inc_p_timely.png", sep = ""), 
        dpi = "print", width = 20, height = 10, units = "cm")

      # extract matched dataset
      tsw_inc_p_timely <- MatchIt::match.data(fit_exp)     
          
    # Reproduction number
      # define formula for exposure as a function of confounders
      formula_exp <- formula(paste(exposures[5]," ~ ", 
        paste(confounders_rn, collapse = " + ")))
    
      # random forest model of exposure as a function of confounders
      fit_exp <- MatchIt::matchit(formula_exp, data = tsw_rn, 
        method = "nearest", ratio = 2, caliper = 0.2, distance = "randomforest")
      
      # extract PS values
      tsw_rn[, "ps_p_timely"] <- fit_exp$distance
        # distribution of PS
        f_hist("ps_p_timely", tsw_rn, c(NA, NA))

      # inspect covariate balance
      out <- cobalt::bal.tab(fit_exp, m.threshold = 0.10, un = TRUE)
      out
      write.csv(out$Balance, paste(dir_output, "/out_ps_bal_rn_p_timely.csv", 
        sep = ""), row.names = TRUE)
      cobalt::love.plot(fit_exp, binary = "std", continuous = "std", 
        m.threshold = 0.10, colors = palette_cb[c(7, 4)])
      ggsave(paste(dir_output, "/out_ps_bal_rn_p_timely.png", sep = ""), 
        dpi = "print", width = 20, height = 10, units = "cm")
    
      # extract matched dataset
      tsw_rn_p_timely <- MatchIt::match.data(fit_exp)     



#.........................................................................................
### ENDS
#.........................................................................................
  
