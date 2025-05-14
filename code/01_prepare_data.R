#...............................................................................
### +++ EFFECT OF SAFE AND DIGNIFIED BURIALS ON EBOV TRANSMISSION IN DRC +++ ###
#...............................................................................

#...............................................................................
## --------- R SCRIPT TO PREPARE DATA AND GRAPH / TABULATE PATTERNS --------- ##
#................................................................................

                              # Written by Francesco Checchi, LSHTM (July 2020)
                              # francesco.checchi@lshtm.ac.uk 


#...............................................................................                           
### Generating time series of health zones and time
#...............................................................................
    
  #...................................       
  ## Determine start and end of analysis period (based on extent of SDB dataset)
    # start
    year_start <- min(sdb_dataset$epi_year, na.rm = T)
    x <- sdb_dataset[which(sdb_dataset$epi_year == year_start), "epi_week"]
    week_start <- min(x, na.rm = T)
    date_start <- min(sdb_dataset$date, na.rm = T)
  
    # end
    year_end <- max(sdb_dataset$epi_year, na.rm = T)
    x <- sdb_dataset[which(sdb_dataset$epi_year == year_end), "epi_week"]
    week_end <- max(x, na.rm = T)    
    date_end <- max(sdb_dataset$date, na.rm = T)

  #...................................    
  ## Create data frame of geographic administrative units
  admin_units <- unique(admin_units[, c("province", "territory", "hz")])
  
  #...................................    
  ## Time series of health zones and days
    # (as date and sequential day numbers; needed for R_n estimation)
  hz <- unique(admin_units[, "hz"])
  dates <- seq(from = as.Date(date_start), to = as.Date(date_end) , by = 1)
  days <- data.frame(dates, c(1:length(dates)))
  colnames(days) <- c("date", "day")
  tsd <- expand.grid(hz, dates)
  colnames(tsd) <- c("hz", "date")
  tsd <- merge(tsd, days, by="date")
  tsd <- tsd[order(tsd[, "hz"], tsd[, "date"]), ]
  tsd[, "hz"] <- as.character(tsd[, "hz"])

  #...................................    
  ## Time series of health zones and epidemiological 1 week- or 2-week increment
  hz <- unique(admin_units[, "hz"])
  weeks <- c(c((year_start * 100 + week_start) : (year_start * 100 + 52) ), 
    c((year_end*100 + 1) : (year_end*100 + week_end)) )
  x <- weeks[seq(2,length(weeks), 2)]
  intervals <- data.frame(weeks1 = weeks, weeks2 = sort(rep(x, 2)) )
  if (t_overall == 1) {tsw <- expand.grid(hz, intervals$weeks1)}
  if (t_overall == 2) {tsw <- expand.grid(hz, x)}
  colnames(tsw) <- c("hz", "weeks")
  tsw$epi_year <- floor(tsw$weeks / 100)
  tsw$epi_week <- tsw$weeks %% 100
  tsw <- tsw[order(tsw$hz, tsw$weeks), ]
  tsw$hz <- as.character(tsw$hz)
  tsw$date <- as.Date(paste(tsw$epi_year, tsw$epi_week, 1), 
    format = "%Y %U %u")
  tsw <- tsw[, colnames(tsw) != "weeks"]
    

#...............................................................................     
### Describing trends in EVD cases
#...............................................................................
  
  #...................................   
  ## Prepare EVD case dataset
    
    # Prepare daily data
    cases_d <- merge(tsd, evd_cases_day, by = c("hz", "date"), all.x = T)

    # Aggregate by analysis interval
    cases_w <- aggregate(cases_d[ , c("confirmed", "probable", "cases", "dead", 
      "nosocomial", "contact")], 
      by = cases_d[, c("hz", "epi_year", "epi_week")], FUN = sum )
      
    # Make sure all HZs and analysis intervals are featured
    cases_w <- merge(tsw, cases_w, by = c("hz", "epi_year", "epi_week"), 
        all.x = T)
      
    # Set cases to zero if value is NA
    cases_w[is.na(cases_w)] <- 0
      
  #...................................   
  ## Plot trends in EVD cases over time, by case type, interval and health zone
      # version 1
    
    # Preparatory steps
      # use weeks data
      df <- cases_w[, c("date", "hz", "confirmed", "probable")]
      
      # identify province
      df <- merge(df, admin_units, by = "hz", all.x = T)
      df[, "hz"] <- paste(df$hz, " (", df$province, ")", sep = "")
      
      # eliminate HZs with no cases
      x <- aggregate(df[, c("confirmed", "probable")], by = list(hz = df$hz), 
        FUN = sum)
      df <- subset(
        df, ! hz %in% x[which(x$confirmed == 0 & x$probable == 0), "hz"] )
      
      # reshape data long
      df <- reshape2::melt(df, id.vars=c("hz", "date", "province", "territory"),
        variable.name = "case_type", value.name = "n_cases")
      df[which(df$n_cases == 0), "n_cases"] <- NA
      
    # Draw plot
    plot <- ggplot(data = df, aes(colour = case_type, fill = case_type, 
      x = date, y = n_cases)) +
      geom_bar(alpha = 0.4, position = "stack", stat = "identity") +
      scale_colour_manual(values = palette_cb[c(8,14)]) +
      scale_fill_manual(values = palette_cb[c(8,14)]) +
      scale_y_continuous("incident reported EVD cases per week") +
      theme_bw() +
      facet_wrap(~hz, ncol = 5) +
      guides(colour = "none") +
      theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1) ) +
      scale_x_date("", expand = c(0, 0) , minor_breaks = NULL, 
        date_breaks = "2 months", date_labels = "%b-%Y" ) +
      labs(fill = "case type:  ")

    # Save plot
    ggsave(paste0(dir_output, "/out_hz_epicurves1.pdf"), height = 28, 
      width = 22, units = "cm", dpi = "print")

  # #...................................   
  # ## Plot trends in EVD cases over time, by case type, interval and health zone
  #     # version 2
  #   
  #   # Preparatory steps
  #     # use weeks data
  #     df <- cases_w[, c("date", "hz", "cases", "probable")]
  #     
  #     # identify province
  #     df <- merge(df, admin_units, by = "hz", all.x = T)
  #     df$hz <- paste0(df$hz, " (", df$province, ")")
  #     
  #     # eliminate HZs with no cases
  #     x <- aggregate(df[, c("cases", "probable")], by = list(hz = df$hz), 
  #       FUN = sum)
  #     df <- subset(
  #       df, ! hz %in% x[which(x$cases == 0), "hz"] )
  #     
  #     # reshape data long
  #     df <- reshape2::melt(df, id.vars=c("hz", "date", "province", "territory"),
  #       variable.name = "case_type", value.name = "n_cases")
  #     df[which(df$n_cases == 0), "n_cases"] <- NA
  #     df$case_type <- as.character(df$case_type)
  #     df[which(df$case_type == "cases"), "case_type"] <-"confirmed and probable"
  #     df[which(is.na(df$n_cases)), "n_cases"] <- 0
  #     
  #   # Draw plot
  #   plot <- ggplot(data = df, aes(colour = case_type, 
  #     x = date, y = n_cases)) +
  #     geom_step(linewidth = 0.75) +
  #     scale_colour_manual(values = palette_cb[c(8, 14)]) +
  #     scale_y_continuous("incident reported EVD cases per week", 
  #       expand = c(0, 0)) +
  #     theme_bw() +
  #     facet_wrap(~hz, ncol = 5) +
  #     theme(legend.position = "bottom",
  #       axis.text.x = element_text(angle = 45, hjust = 1) ) +
  #     scale_x_date("", expand = c(0, 0) , minor_breaks = NULL, 
  #       date_breaks = "2 months", date_labels = "%b-%Y" ) +
  #     labs(colour = "")
  # 
  #   # Save plot
  #   ggsave(paste0(dir_output, "/out_hz_epicurves2.pdf"), height = 28, 
  #     width = 22, units = "cm", dpi = "print")
    
    
  #...................................   
  ## Tabulate trends in EVD cases over time, by health zone, and screen out 
    # health zones with overly sparse data
    
    # Create table
      # aggregate
      out <- aggregate(cases_w[, c("confirmed", "probable", "cases")], 
        by = list(hz = cases_w$hz), FUN = sum)
      out <- merge(out, unique(admin_units[, c("province", "hz")]), by = "hz",
        all.x = T)      
      out <- out[order(out[, "province"], out[, "hz"]), ]
      out <- out[, c("province", "hz", "confirmed", "probable", "cases")]
      
      # apply decision rule: at least n cumulative cases
      out[, "eligible"] <- ifelse(out$cases >= cases_threshold, "yes", "no")

      # compute column totals and percentages
      out <- rbind(out, 
        c("Totals", colSums(out[, c("confirmed", "probable", "cases")]),"", ""))
      
      # eliminate health zones with no cases
      out <- subset(out, cases > 0)
      out[, "percent"] <- round(100 * as.numeric(out$cases) / 
          sum(as.numeric(out$cases[1:(nrow(out)-1)])), digits = 1)

    # Identify health zones that are eligible
    hz_ok <- out[which(out$eligible == "yes"), "hz"]
          
    # Save as csv    
    write.csv(out, paste0(dir_output, "/out_table_cases_by_hz.csv"), 
      row.names = F)
        
    # Save as Word table
    colnames(out) <- c("Province", "Health zone", "Confirmed cases", 
      "Probable cases", "Total cases", "Eligible", "Percent of all cases")
    out <- flextable(out)      
    out <- autofit(out)
    flextable::save_as_rtf(out, 
      path = paste(dir_output, "/out_table_cases_by_hz.rtf", sep =""))

#...............................................................................     
### Computing transmission outcome 1: Effective reproduction number
#...............................................................................
          
  #...................................
  ## Estimate effective reproduction number R_n for each HZ, by interval
    
    # Prepare daily data
    cases_d <- merge(tsd, evd_cases_day, by = c("hz", "date"), all.x = T)
      
      # select variables of interest
      cases_d <- cases_d[, c("hz", "date", "cases")]
      cases_d <- cases_d[order(cases_d$hz, cases_d$date), ]
    
      # set cases to zero if value is NA
      cases_d[which(is.na(cases_d$cases)), "cases"] <- 0

    # Initialise output  
    r_est_d <- data.frame()

    # For each health zone eligible for analysis...
    for (i in hz_ok) {

      # progress statement
      print(paste("now estimating R_n for ", i, " health zone", sep = "") )
            
      # select hz time series
      cases_d_i <- subset(cases_d, hz == i)$cases
      
      # identify position of first non-zero value
      x <- min(which(cases_d_i != 0))
      
      # set time windows
      t_start <- seq(x + 1, length(cases_d_i) - t_overall * 7)
      t_end <- t_start + t_overall * 7
      
      # estimate R_n (Wallinga-Teunis method)
       r_est_i <- wallinga_teunis(
        incid = cases_d_i,
        method = "parametric_si", 
        config = list(
          mean_si = si_mean,
          std_si = si_sd,
          method = "parametric_si",
          n_sim = boot_rt,
          t_start = t_start,
          t_end = t_end
        )
      )
      
      # extract useful parameters
      r_est_i <- as.data.frame(r_est_i[["R"]])[, c("t_start","Mean(R)","Std(R)", 
        "Quantile.0.025(R)", "Quantile.0.975(R)")]
      
      # add to results
      r_est_i[, "hz"] <- i
      r_est_d <- rbind(r_est_d, r_est_i)
    }

    # Rename columns
    colnames(r_est_d) <- c("day", "r_mean", "r_se", "r_lci", "r_uci", "hz")
    
    # Map consecutive day numbers to dates
    r_est_d <- merge(r_est_d, subset(tsd, hz %in% hz_ok), by=c("hz", "day"), 
      all = T)
  
  #...................................
  ## Eliminate estimates at the start of each local outbreak, if desired
  if (omit_initial) {  
    # Sort R_n estimates
    r_est_d <- r_est_d[order(r_est_d$hz, r_est_d$date), ]
    
    # Set R estimates to null for 2 weeks (~ 1 serial interval) from start
      # of each non-missing segment of time series
    for (i in (nrow(r_est_d)-1):2) {
      if ((! is.na(r_est_d[i, "r_mean"])) &
        (is.na(r_est_d[i-1, "r_mean"])) & 
        (r_est_d[i-1, "hz"] == r_est_d[i, "hz"])
      ) {
        r_est_d[i:(i+round(si_mean, 0)), 
          c("r_mean", "r_se", "r_lci", "r_uci")] <- NA
      }       
    }
  }

  #...................................
  ## Convert to interval timescale and merge with cases

    # Convert to values over interval
    r_est_w <- merge(subset(tsw, hz %in% hz_ok), r_est_d, by = c("hz", "date"), 
      all.x = T )

    # Merge R_n estimates with cases
    r_est_w <- merge(r_est_w, 
      subset(cases_w, hz %in% hz_ok)[, c("hz", "date", "cases")], 
      by = c("hz", "date"), all.x = T)  

    # Save estimates
    write.csv(r_est_w, paste0(dir_output, "/out_rn_w.csv"), row.names = F)

  #...................................
  ## Plot median and range of R_n point estimates for each health zone, 
    # along with number of days estimates are available for
    # Prepare data  
    df <- r_est_d
    df <- subset(df, ! is.na(r_mean))
    
    # Compute number of days with estimate per HZ
    df[, "n_obs"] <- 1
    x <- aggregate(df$n_obs, by = list(df$hz), FUN = sum)
    colnames(x) <- c("hz", "n_obs")
    df <- df[, colnames(df) != "n_obs"]
    df <- merge(df, x, by = "hz", all.x = T)
          
    # Plot
    plot <- ggplot(data = df, aes(hz, r_mean)) +
      geom_boxplot(colour = palette_cb[11], fill = palette_cb[11], alpha =0.5) +
      geom_hline(aes(yintercept = 1), colour=palette_cb[4], linetype="dashed") +
      scale_y_continuous("estimated daily effective reproduction number", 
        limits = c(-0.2, 4)) +
      scale_x_discrete("health zone") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank()) +
      geom_text(data = x, aes(x = hz, y = rep(0, nrow(x)), 
        label = as.character(n_obs) ), nudge_y = -0.2, 
        colour = palette_cb[5] )

    # Save plot
    ggsave(paste(dir_output, "/out_rn_dist_by_hz.pdf", sep = ""), 
      height = 15, width = 23, units = "cm", dpi = "print")
    
  #...................................
  ## Plot point estimate and CI of R_n point estimates for each HZ, by day 
    # Prepare data  
    df <- r_est_d
    df <- merge(df, admin_units, by = "hz", all.x = T)
    df$hz <- paste(df$hz, " (", df$province, ")", sep = "")
    
    # Draw plot
    plot <- ggplot(data = df, aes(x = date, y = r_mean)) +
      geom_step(linewidth = 0.75, colour = palette_cb[11]) +
      geom_stepribbon(aes(ymin  = r_lci, ymax = r_uci), fill = palette_cb[11],
        alpha = 0.3) +
      geom_hline(aes(yintercept = 1), colour=palette_cb[4], linetype = "21") +
      scale_y_continuous("estimated daily effective reproduction number", 
        expand = c(0, 0), trans = scales::pseudo_log_trans(), 
        breaks = c(0, 1, 2, 5, 10)) +
      theme_bw() +
      facet_wrap(~hz, ncol = 5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1) ) +
      scale_x_date("", expand = c(0, 0) , minor_breaks = NULL, 
        date_breaks = "2 months", date_labels = "%b-%Y" )

    # Save plot
    ggsave(paste(dir_output, "/out_rn_hz_evolution_day.pdf", sep = ""), 
      height = 22, width = 22, units = "cm", dpi = "print")

  #...................................
  ## Plot point estimate and CI of R_n point estimates for each HZ, by interval 
    # Prepare data  
    df <- r_est_w
    df <- merge(df, admin_units, by = "hz", all.x = T)
    df$hz <- paste(df$hz, " (", df$province, ")", sep = "")
    
    # Draw plot
    plot <- ggplot(data = df, aes(x = date, y = r_mean)) +
      geom_step(linewidth = 0.75, colour = palette_cb[11]) +
      geom_stepribbon(aes(ymin  = r_lci, ymax = r_uci), fill = palette_cb[11],
        alpha = 0.3) +
      geom_hline(aes(yintercept = 1), colour=palette_cb[4], linetype = "21") +
      scale_y_continuous(
        paste0("estimated effective reproduction number (", 
          t_overall, "-week intervals)"), 
        expand = c(0, 0), trans = scales::pseudo_log_trans(), 
        breaks = c(0, 1, 2, 5, 10)) +
      theme_bw() +
      facet_wrap(~hz, ncol = 5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1) ) +
      scale_x_date("", expand = c(0, 0) , minor_breaks = NULL, 
        date_breaks = "2 months", date_labels = "%b-%Y" )

    # Save plot
    ggsave(paste(dir_output, "/out_rn_hz_evolution_interval.pdf", sep = ""), 
      height = 22, width = 22, units = "cm", dpi = "print")
        
#...............................................................................     
### Computing transmission outcome 2: Incidence change between time windows
    # (relative period incidence), as a different transmission outcome 
    # by time window - centered two weeks apart or more, to avoid overlap
    # two weeks distance between windows because that is approximately 
    # the length of a serial interval
    # low number of cases in current window compared to previous window 
    # indicates transmission is slowing, and vice versa
#...............................................................................
        
  #...................................
  ## Compute numbers of cases in current and previous time window 
    
    # Compute current and previous window incidence
    cases_w[, c("cases_prev", "cases_curr")] <- 
      f_window_transm(window_size_f = window_transmission)
    
    # Also create weights (sqrt of sum of the 2 time windows), to try later
    cases_w$cases_wt <- (cases_w$cases_prev + cases_w$cases_curr)^0.5


  #...................................
  ## Plot evolution of incidence change for each HZ, by time intervals 
    # Prepare data  
    df <- cases_w
    df <- merge(df, admin_units, by = "hz", all.x = T)
    df <- subset(df, df$hz %in% hz_ok)    
    df$hz <- paste(df$hz, " (", df$province, ")", sep = "")
    df$pi_it <- NA
    df$cases_2w <- df$cases_prev + df$cases_curr
    x <- which(df$cases_2w > 0 & ! is.na(df$cases_2w) )
    df[x, "pi_it"] <- df[x, "cases_curr"] / df[x, "cases_2w"]
    
    # Draw plot
    plot <- ggplot(data = df, aes(x = date, y = pi_it)) +
      geom_step(linewidth = 0.75, colour = palette_cb[5]) +
      geom_hline(aes(yintercept = 0.5), colour=palette_cb[14], linetype = "21") +
      scale_y_continuous(name = expression(pi[i][t]), 
        expand = c(0.02, 0.02), trans = scales::pseudo_log_trans(), 
        breaks = c(0, 0.5, 1)) +
      theme_bw() +
      facet_wrap(~hz, ncol = 5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1) ) +
      scale_x_date("", expand = c(0, 0) , minor_breaks = NULL, 
        date_breaks = "2 months", date_labels = "%b-%Y" )

    # Save plot
    ggsave(paste(dir_output, "/out_inc_hz_evolution.pdf", sep = ""), height =25, 
      width = 22, units = "cm", dpi = "print")

    
    
#...............................................................................                         
### Preparing any other datasets that need further data management at this stage
#............................................................................... 
    
  #...................................   
  ## SDB performance (exposure)
    # Exclude instances with EBOV status = negative (definite non-case)
    sdb_dataset <- subset(sdb_dataset, evd_status != "negative")
    
    # Successful SDBs
      # identify number of SDBs (any done)
      sdb_dataset$n_sdb_all <- 1
      
      # identify number of SDBs with outcome of interest
        # (success or failure excluding SDBs not needed)
      sdb_dataset$n_outcome <- ifelse(sdb_dataset$outcome_lshtm %in% 
          c("failure", "success"), 1, 0)

      # identify successful SDBs
      sdb_dataset$n_success <- NA
      sdb_dataset[which(sdb_dataset$n_outcome == 1 & 
        sdb_dataset$outcome_lshtm == "success"), "n_success"] <- 1 
      sdb_dataset[which(sdb_dataset$n_outcome == 1 & 
        sdb_dataset$outcome_lshtm == "failure"), "n_success"] <- 0 

      # identify components of success
      sdb_dataset$n_secured <- NA
      sdb_dataset[which(sdb_dataset$n_outcome == 1 & 
        sdb_dataset$secured_body == "yes"), "n_secured"] <- 1 
      sdb_dataset[which(sdb_dataset$n_outcome == 1 & 
        sdb_dataset$secured_body != "yes"), "n_secured"] <- 0 

      sdb_dataset$n_disinfected <- NA
      sdb_dataset[which(sdb_dataset$n_outcome == 1 & 
        sdb_dataset$disinfected_site == "yes"), "n_disinfected"] <- 1 
      sdb_dataset[which(sdb_dataset$n_outcome == 1 & 
        sdb_dataset$disinfected_site != "yes"), "n_disinfected"] <- 0 

      sdb_dataset$n_specimen <- NA
      sdb_dataset[which(sdb_dataset$n_outcome == 1 & 
        sdb_dataset$specimen_collected == "yes"), "n_specimen"] <- 1 
      sdb_dataset[which(sdb_dataset$n_outcome == 1 & 
        sdb_dataset$specimen_collected != "yes"), "n_specimen"] <- 0 
                  
    # Timely SDBs (delay < 24h)
      # identify number of SDBs with information on delay
      sdb_dataset$n_delay <- ifelse(sdb_dataset$delay_sdb %in% 
          c("<= 24h", "> 24h"), 1, 0)

      # identify timely SDBs
      sdb_dataset$n_timely <- NA
      sdb_dataset[which(sdb_dataset$n_delay == 1 & 
        sdb_dataset$delay_sdb == "> 24h"), "n_timely"] <- 0    
      sdb_dataset[which(sdb_dataset$n_delay == 1 & 
        sdb_dataset$delay_sdb == "<= 24h"), "n_timely"] <- 1    

    # # Successful AND timely SDBs
    #   # identify number of SDBs with information on outcome and delay
    #   sdb_dataset$n_both <- 0
    #   sdb_dataset[which(sdb_dataset$n_outcome == 1 & sdb_dataset$n_delay == 1), 
    #     "n_both"] <- 1

      # # identify SDBs that were successful and timely
      # sdb_dataset$n_both_ok <- NA
      # sdb_dataset[which(sdb_dataset$n_outcome == 1 & sdb_dataset$n_delay == 1), 
      #   "n_both_ok"] <- 0    
      # sdb_dataset[which(sdb_dataset$n_outcome == 1 & sdb_dataset$n_delay == 1 
      #   & sdb_dataset$n_success == 1 & sdb_dataset$n_timely == 1), 
      #   "n_both_ok"] <- 1
    
    # Describe overlap of key SDB success components
    df <- subset(sdb_dataset, outcome_lshtm %in% c("success", "failure"))
    df <- df[, c("secured_body", "disinfected_site", "buried_sdb")]
    df[df[] == "yes"] <- T
    df[df[] == "no"] <- F
    df[df[] == "unclear"] <- NA
    df <- na.omit(df)
    for (i in 1:ncol(df)) { df[, i] <- as.logical(df[, i])}
    colnames(df) <- c("body secured", "site disinfected", 
      "decedent buried")
    plot <- ggvenn(df, fill_color = palette_cb[c(2,8,14)], fill_alpha = 0.3)
    ggsave(paste(dir_output, "/out_venn_sdb_comps.pdf", sep = ""), height = 20, 
      width = 25, units = "cm", dpi = "print")
    
    # Aggregate by time interval
    sdb_w <- aggregate(sdb_dataset[ , c("n_sdb_all", "n_outcome", "n_success", 
      "n_delay", "n_timely", "n_secured", "n_disinfected", "n_specimen")], 
      by = sdb_dataset[, c("hz", "epi_year", "epi_week")], FUN = sum, na.rm = T)
    if (t_overall == 2) {
      sdb_w$weeks1 <- sdb_w$epi_year * 100 + sdb_w$epi_week
      sdb_w <- merge(sdb_w, intervals, by = "weeks1", all.x = T)
      sdb_w <- aggregate(sdb_w[ , c("n_sdb_all", "n_outcome", "n_success", 
        "n_delay", "n_timely", "n_secured", "n_disinfected", "n_specimen")], 
        by = sdb_w[, c("hz", "weeks2")], FUN = sum, na.rm = T)
      colnames(sdb_w)[colnames(sdb_w) == "weeks2"] <- "weeks"
      sdb_w$epi_year <- floor(sdb_w$weeks / 100)
      sdb_w$epi_week <- sdb_w$weeks %% 100
    }
    sdb_w <- sdb_w[order(sdb_w$hz, sdb_w$epi_year, sdb_w$epi_week), ]
    sdb_w <- sdb_w[, colnames(sdb_w) != "weeks"]
    
  #...................................   
  ## Vaccination status of cases
    
    # Aggregate to time interval
    vacc_w <- vaccinated_week
    if (t_overall == 2) {
      vacc_w$weeks1 <- vacc_w$epi_year * 100 + vacc_w$epi_week
      vacc_w <- merge(vacc_w, intervals, by = "weeks1", all.x = T)
      vacc_w <- aggregate(vacc_w[ , c("n_vaccinated", "n_vacc_status")], 
        by = vacc_w[, c("hz", "weeks2")], FUN = sum, na.rm = T)
      colnames(vacc_w)[colnames(vacc_w) == "weeks2"] <- "weeks"      
      vacc_w$epi_year <- floor(vacc_w$weeks / 100)
      vacc_w$epi_week <- vacc_w$weeks %% 100
    }
      
    # Make sure all HZs and dates are featured
    vacc_w <- vacc_w[, colnames(vacc_w) != "weeks"]
    vacc_w <- merge(tsw, vacc_w, by = c("hz", "epi_year", "epi_week"), all.x= T)

    # Set cases to zero if value is NA
    vacc_w[is.na(vacc_w)] <- 0


  #...................................   
  ## Insecurity datasets: change time interval if needed
  if (t_overall == 2) {
    
    # Attacks against EVD response
    attacks_evd_weeks$weeks1 <- attacks_evd_weeks$epi_year * 100 + 
      attacks_evd_weeks$epi_week
    attacks_evd_weeks <- merge(attacks_evd_weeks, intervals, 
      by = "weeks1", all.x = T)
    attacks_evd_weeks <- aggregate(attacks_evd_weeks[ , 
      c("against_evd", "suspend_activities")], 
      by = attacks_evd_weeks[, c("hz", "weeks2")], 
      FUN = sum, na.rm = T)
    colnames(attacks_evd_weeks)[colnames(attacks_evd_weeks) == 
        "weeks2"] <-"weeks"      
    attacks_evd_weeks$epi_year <- floor(attacks_evd_weeks$weeks / 100)
    attacks_evd_weeks$epi_week <- attacks_evd_weeks$weeks %% 100
    attacks_evd_weeks <- 
      attacks_evd_weeks[, colnames(attacks_evd_weeks) != "weeks"]
    
    # ACLED insecurity
    insecurity_acled_weeks$weeks1 <- insecurity_acled_weeks$epi_year * 100 + 
      insecurity_acled_weeks$epi_week
    insecurity_acled_weeks <- merge(insecurity_acled_weeks, intervals, 
      by = "weeks1", all.x = T)
    insecurity_acled_weeks <- aggregate(insecurity_acled_weeks[ , 
      c("n_events", "fatalities")], 
      by = insecurity_acled_weeks[, c("territory", "weeks2")], 
      FUN = sum, na.rm = T)
    colnames(insecurity_acled_weeks)[colnames(insecurity_acled_weeks) == 
        "weeks2"] <-"weeks"      
    insecurity_acled_weeks$epi_year <- floor(insecurity_acled_weeks$weeks / 100)
    insecurity_acled_weeks$epi_week <- insecurity_acled_weeks$weeks %% 100
    insecurity_acled_weeks <- 
      insecurity_acled_weeks[, colnames(insecurity_acled_weeks) != "weeks"]
    
    # ETC presence
    etc_presence$weeks1 <- etc_presence$epi_year * 100 + 
      etc_presence$epi_week
    etc_presence <- merge(etc_presence, intervals, 
      by = "weeks1", all.x = T)
    colnames(etc_presence)[colnames(etc_presence) == 
        "weeks2"] <-"weeks"      
    etc_presence$epi_year <- floor(etc_presence$weeks / 100)
    etc_presence$epi_week <- etc_presence$weeks %% 100
    etc_presence <- etc_presence[, ! colnames(etc_presence) 
      %in% c("weeks", "weeks1")]
  }
    
                 
#...............................................................................                            
### Merging transmission outcomes, exposures, possible confounders and 
  # population denominators onto time series
#............................................................................... 
    
  #...................................   
  ## Merge admin units and population
  tsw <- merge(tsw, admin_units, by = "hz", all.x = T)
  tsw <- merge(tsw, population, by = "hz", all.x = T)

  #...................................   
  ## Merge exposure dataset
  tsw <- merge(tsw, sdb_w, by = c("hz", "epi_year", "epi_week"), all.x = T)
    
  #...................................   
  ## Merge confounder datasets
  tsw <- merge(tsw, vacc_w[, colnames(vacc_w) != "date"], 
    by = c("hz", "epi_year", "epi_week"), all.x = T)
  tsw <- merge(tsw, attacks_evd_weeks, by = c("hz", "epi_year", "epi_week"), 
    all.x = T)
  tsw <- merge(tsw, insecurity_acled_weeks, 
    by = c("territory", "epi_year", "epi_week"), all.x = T)
  tsw <- merge(tsw, 
    etc_presence[, ! colnames(etc_presence) %in% c("month", "year")], 
    by = c("hz", "epi_year", "epi_week"), all.x = T)
  tsw <- merge(tsw, health_facilities, by = "hz", all.x = T)
  tsw <- merge(tsw, cell_coverage, by = "hz", all.x = T)
  tsw <- merge(tsw, radio_coverage, by = "hz", all.x = T)
  tsw <- merge(tsw, road_coverage, by = "hz", all.x = T)
  tsw <- merge(tsw, mines, by = "hz", all.x = T)
    
  #...................................   
  ## Merge transmission outcome data
    # Net reproduction number
    tsw <- merge(tsw, r_est_w[, ! colnames(r_est_w) %in% c("cases", "date")], 
      by = c("hz", "epi_year", "epi_week"), all.x = T)
    
    # EVD cases including relative period incidence 
      # (also contains some explanatory variables - see below)
    tsw <- merge(tsw, cases_w[, ! colnames(cases_w) %in% c("date")], 
      by = c("hz", "epi_year", "epi_week"), all.x = T)

  #...................................   
  ## Restrict observations to eligible health zones
  tsw <- subset(tsw, hz %in% hz_ok)

       
#...............................................................................                            
### Preparing SDB-related explanatory variables
#...............................................................................
    
  #...................................   
  ## Compute SDB exposures
    
    # Prepare data
    tsw <- tsw[order(tsw$hz, tsw$epi_year, tsw$epi_week), ]
    hz_ok <- sort(hz_ok)    
    tsw$n_sdb_all <- ifelse(is.na(tsw$n_sdb_all), 0, tsw$n_sdb_all)
    # tsw$n_sdb_ok_rate <- tsw$n_both_ok * 100000 / tsw$pop
  
    # Proportion of successful SDBs
    tsw[, c("p_success_prev", "p_success_curr")] <-
      f_window_transm(data_f = tsw, hz_f = hz_ok, vars_f = 
        c("n_success", "n_outcome"), window_size_f = window_transmission)
    
    # proportion of timely SDBs  
    tsw[, c("p_timely_prev", "p_timely_curr")] <-
      f_window_transm(data_f = tsw, hz_f = hz_ok, vars_f = 
        c("n_timely", "n_delay"), window_size_f = window_transmission)

    # weekly ratio of any SDBs done to population per 100,000 person-weeks,
      # by time window; here use all SDBs (coverage proxy)
    tsw[, c("r_sdb_to_pop_prev", "r_sdb_to_pop_curr")] <-
      f_window_transm(data_f = tsw, hz_f = hz_ok, vars_f = 
        c("n_sdb_all", "pop"), window_size_f = window_transmission) * 100000

  #...................................   
  ## Visualise SDB coverage by comparing SDBs for confirmed/unknown EBOV deaths
    # to reported EBOV deaths
    
    # Prepare dataset of SDBs on confirmed cases by time interval
    sdb_conf <- subset(sdb_dataset, evd_status == "positive")
    sdb_w_conf <- aggregate(sdb_conf$n_sdb_all, 
      by = sdb_conf[, c("hz", "epi_year", "epi_week")], FUN = sum, na.rm = T)
    colnames(sdb_w_conf) <- c("hz", "epi_year", "epi_week", "n_sdb_conf")
    if (t_overall == 2) {
      sdb_w_conf$weeks1 <- sdb_w_conf$epi_year * 100 + 
        sdb_w_conf$epi_week
      sdb_w_conf <- merge(sdb_w_conf, intervals, 
        by = "weeks1", all.x = T)
      sdb_w_conf <- aggregate(list(n_sdb_conf = sdb_w_conf$n_sdb_conf), 
        by = sdb_w_conf[, c("hz", "weeks2")], 
        FUN = sum, na.rm = T)
      colnames(sdb_w_conf)[colnames(sdb_w_conf) == 
          "weeks2"] <-"weeks"      
      sdb_w_conf$epi_year <- floor(sdb_w_conf$weeks / 100)
      sdb_w_conf$epi_week <- sdb_w_conf$weeks %% 100
      sdb_w_conf <- 
        sdb_w_conf[, colnames(sdb_w_conf) != "weeks"]
    }

    # Prepare dataset of SDBs on unknown-status cases by interval
    sdb_unk <- subset(sdb_dataset, evd_status == "unclear")
    sdb_w_unk <- aggregate(sdb_unk$n_sdb_all, 
      by = sdb_unk[, c("hz", "epi_year", "epi_week")], FUN = sum, na.rm = T)
    colnames(sdb_w_unk) <- c("hz", "epi_year", "epi_week", "n_sdb_unk")
    if (t_overall == 2) {
      sdb_w_unk$weeks1 <- sdb_w_unk$epi_year * 100 + 
        sdb_w_unk$epi_week
      sdb_w_unk <- merge(sdb_w_unk, intervals, 
        by = "weeks1", all.x = T)
      sdb_w_unk <- aggregate(list(n_sdb_unk = sdb_w_unk$n_sdb_unk), 
        by = sdb_w_unk[, c("hz", "weeks2")], 
        FUN = sum, na.rm = T)
      colnames(sdb_w_unk)[colnames(sdb_w_unk) == 
          "weeks2"] <-"weeks"      
      sdb_w_unk$epi_year <- floor(sdb_w_unk$weeks / 100)
      sdb_w_unk$epi_week <- sdb_w_unk$weeks %% 100
      sdb_w_unk <- 
        sdb_w_unk[, colnames(sdb_w_unk) != "weeks"]
    }
        
    # Prepare dataset of confirmed EBOV deaths
    cases_w_dead <- cases_w[, c("hz", "epi_year", "epi_week", "dead")]
    cases_w_dead <- subset(cases_w_dead, dead > 0)
    
    # Merge all into one
    conf <- merge(sdb_w_conf, sdb_w_unk, 
      by = c("hz", "epi_year", "epi_week"), all = T)
    conf <- merge(conf, cases_w_dead, 
      by = c("hz", "epi_year", "epi_week"), all = T)
    # conf <- subset(conf, hz %in% hz_ok)
    
    # Add dates for plotting and prepare a bit more
    conf <- merge(tsw[, c("province", "hz", "epi_year", "epi_week", "date")],
      conf, by = c("hz", "epi_year", "epi_week"), all.x = T) 
    conf$date <- as.Date(conf$date)
    conf[which(is.na(conf$n_sdb_conf)), "n_sdb_conf"] <- 0
    conf[which(is.na(conf$n_sdb_unk)), "n_sdb_unk"] <- 0
    conf[which(is.na(conf$dead)), "dead"] <- 0
    conf$hz <- paste(conf$hz, " (", conf$province, ")", sep = "")

    # Make annotations for each facet (n dead, n SDBs)
      # reported deaths
      txt_dead <- aggregate(conf$dead, by = list(conf$hz), FUN = sum)
      colnames(txt_dead) <- c("hz", "dead")
      txt_dead$dead <- paste(txt_dead$dead, "deaths in line list")
      
      # SDBs on confirmed cases
      txt_n_sdb_conf <- aggregate(conf$n_sdb_conf, by = list(conf$hz),FUN = sum)
      colnames(txt_n_sdb_conf) <- c("hz", "n_sdb_conf")
      txt_n_sdb_conf$n_sdb_conf <- paste(txt_n_sdb_conf$n_sdb_conf, 
        "SDBs - confirmed EBOV+")

      # SDBs on unknown-status cases
      txt_n_sdb_unk <- aggregate(conf$n_sdb_unk, by = list(conf$hz), FUN = sum)
      colnames(txt_n_sdb_unk) <- c("hz", "n_sdb_unk")
      txt_n_sdb_unk$n_sdb_unk <- paste(txt_n_sdb_unk$n_sdb_unk, 
        "SDBs - unconfirmed")
  
    # Plot
    plot <- ggplot(data = conf, aes(x = date)) +
      geom_bar(aes(y = dead), stat = "identity", colour = palette_cb[6], 
        fill = palette_cb[1], alpha = 0.2) +
      geom_step(aes(y = n_sdb_conf), colour = palette_cb[2]) +
      geom_step(aes(y = n_sdb_unk), colour = palette_cb[12]) +
      theme_bw() +
      facet_wrap(~hz, ncol = 5) +
      theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1) ) +
      scale_x_date("", expand = c(0, 0) , minor_breaks = NULL, 
        date_breaks = "2 months", date_labels = "%b-%Y" ) +
      scale_y_continuous(expand = c(0, 0), name =
        "confirmed EVD deaths / SDBs with confirmed or unconfirmed EBOV status") +
      geom_text(x = as.Date("2018-09-01"), y = 70, aes(label = dead), 
        data = txt_dead, colour = palette_cb[6], hjust = 0, size = 2.5) +
      geom_text(x = as.Date("2018-09-01"), y = 60, aes(label = n_sdb_conf), 
        data = txt_n_sdb_conf, colour = palette_cb[2], hjust = 0, size = 2.5) +
      geom_text(x = as.Date("2018-09-01"), y = 50, aes(label = n_sdb_unk), 
        data = txt_n_sdb_unk, colour = palette_cb[12], hjust = 0, size = 2.5)
    
    # Save plot
    ggsave(paste(dir_output, "/out_sdb_coverage.pdf", sep = ""), height = 23, 
      width = 22, units = "cm", dpi = "print")


  #...................................   
  ## Plot correlation between weekly cases and weekly SDBs done
    
    # Create plot
    plot <- ggplot(data = tsw, aes(x = cases_curr, y = n_sdb_all,
      colour = province, fill = province)) + 
      geom_jitter(alpha = 0.5, width = 0.1, height = 0.1) +
      geom_smooth(method = "gam") +             
      scale_color_manual(values = palette_cb[c(3,12)]) +
      scale_fill_manual(values = palette_cb[c(3,12)]) +
      scale_y_continuous("weekly number of SDBs performed",
        trans = pseudo_log_trans(base = 2), breaks = seq(0, 100, by = 10) ) +
      scale_x_continuous("weekly number of EVD cases",
        trans = pseudo_log_trans(base = 2), breaks = seq(0, 60, by = 10) ) +
      expand_limits(x = 0, y = 0) +
      theme_bw() +
      theme(legend.position = "none") +
      facet_wrap(.~province, ncol = 2)
    
    # Save plot
    ggsave(paste(dir_output, "/out_sdb_cases_corr.pdf", sep = ""), height = 15, 
      width = 25, units = "cm", dpi = "print")
      

#...............................................................................                            
### Preparing explanatory variables based on EVD line list
#...............................................................................
    
  #...................................   
  ## Percentage of cases with nosocomial transmission, by time window
    
    if (! smooth_p_nosocomial) {
      # Compute variable
      tsw[, c("p_nosocomial_prev", "p_nosocomial_curr")] <-
        f_window_transm(data_f = tsw, window_size_f = window_transmission,
          vars_f = c("nosocomial", "cases"), hz_f = hz_ok)
    }
    
    if (smooth_p_nosocomial) {
      # Compute and smooth variable
      tsw <- tsw[order(tsw$hz, tsw$epi_year, tsw$epi_week), ]
      tsw$p_nosocomial <- NA
      for (i in hz_ok) {
  
        # select data needed
        x <- tsw[which(tsw$hz == i), c("nosocomial", "cases")]
        x$week <- 1:nrow(x)
  
        # smooth variable by fitting a binomial additive model
        x$other <- x$cases - x$nosocomial
        x <- x[, c("nosocomial", "other", "cases", "week")]
        form <- as.formula(paste0("cbind(", "nosocomial", ",",
          "other",")", "~", "pb(week)"))
        fit <- gamlss(form, family = BI, data = na.omit(x))
        tsw[which(tsw$hz == i), "p_nosocomial"] <- inv.logit(predict(fit))
      }
  
      # Visualise
      ggplot(data = tsw, aes(x = date)) +
        geom_point(aes(y = nosocomial/cases, size = cases), alpha = 0.25,
          colour = palette_cb[6]) +
        geom_line(aes(y = p_nosocomial), linewidth = 1, colour = palette_cb[6]) +
        theme_bw() +
        facet_wrap(~hz, ncol = 5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
        scale_x_date("", minor_breaks = NULL,
          date_breaks = "2 months", date_labels = "%b-%Y" ) +
        scale_y_continuous("proportion of nosocomial cases",
          labels = percent) +
        scale_size_continuous("number of cases")
  
      # Compute variable by time window
      tsw[, c("p_nosocomial_prev", "p_nosocomial_curr")] <-
        f_window_transm(data_f = tsw, window_size_f = window_transmission,
          vars_f = "p_nosocomial", hz_f = hz_ok, op_f = "mean")
    }
   
  #...................................   
  ## Percentage of cases with any known link to another known case, by window
    
    if (! smooth_p_epi_link) {
      # Compute variable
      tsw[, c("p_epi_link_prev", "p_epi_link_curr")] <-
        f_window_transm(data_f = tsw, window_size_f = window_transmission,
          vars_f = c("contact", "cases"), hz_f = hz_ok)
    }
    
    if (smooth_p_epi_link) {
      # Compute and smooth variable
      tsw <- tsw[order(tsw$hz, tsw$epi_year, tsw$epi_week), ]
      tsw$p_epi_link <- NA
      for (i in hz_ok) {
  
        # select data needed
        x <- tsw[which(tsw$hz == i), c("contact", "cases")]
        x$week <- 1:nrow(x)
  
        # smooth variable by fitting a binomial additive model
        x$other <- x$cases - x$contact
        x <- x[, c("contact", "other", "cases", "week")]
        form <- as.formula(paste0("cbind(", "contact", ",",
          "other",")", "~", "pb(week)"))
        fit <- gamlss(form, family = BI, data = na.omit(x))
        tsw[which(tsw$hz == i), "p_epi_link"] <- inv.logit(predict(fit))
      }
  
      # Visualise
      ggplot(data = tsw, aes(x = date)) +
        geom_point(aes(y = contact/cases, size = cases), alpha = 0.25,
          colour = palette_cb[12]) +
        geom_line(aes(y = p_epi_link), linewidth = 1, colour = palette_cb[12]) +
        theme_bw() +
        facet_wrap(~hz, ncol = 5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
        scale_x_date("", minor_breaks = NULL,
          date_breaks = "2 months", date_labels = "%b-%Y" ) +
        scale_y_continuous("proportion of cases with an epidemiological contact",
          labels = percent) +
        scale_size_continuous("number of cases")
  
      # Compute variable by time window
      tsw[, c("p_epi_link_prev", "p_epi_link_curr")] <-
        f_window_transm(data_f = tsw, window_size_f = window_transmission,
          vars_f = "p_epi_link", hz_f = hz_ok, op_f = "mean")
    }

  #...................................   
  ## Vaccination coverage based on case-coverage formula and 
    # assumed vaccine effectiveness
    
    # # Setting all NA values to zero (makes no difference)
    # tsw[which(is.na(tsw$n_vaccinated)), "n_vaccinated"] <- 0
    # tsw[which(is.na(tsw$n_vacc_status)), "n_vacc_status"] <- 0

    if (!smooth_vcp) {
      # Vaccination coverage, by time window and lag of interest
      for (i in windows_other) {
  
        for (j in lags) {
  
          # first calculate vaccination coverage among cases
          out_ij <- f_window_other(var_f = "n_vaccinated", window_size = i,
            lag_f = j) / f_window_other(var_f = "n_vacc_status",
              window_size = i, lag_f = j)
  
          # then estimate vaccination coverage in population
          x <- out_ij / (1 - (1 - out_ij) * ve)
          x[which(x %in% NaN)] <- NA
          tsw[, paste("vcp_", i, "w_lag", j, sep="")] <- x
        }
      }
    }
    
    if (smooth_vcp) {
      # Initialise new variables
      tsw[, c("vcc", "vcc_sm", "vcp")] <- NA
      tsw <- tsw[order(tsw[, "hz"], tsw[, "epi_year"], tsw[, "epi_week"]), ]
  
      # Smooth vaccination coverage among cases (VCC) by Health Zone
      for (i in hz_ok) {
  
        # select data needed
        x <- tsw[which(tsw$hz == i), c("n_vaccinated", "n_vacc_status")]
        x$week <- 1:nrow(x)
  
        # add fake data to week 1 (VCC must be zero at the start!)
        x[1, "n_vaccinated"] <- 0
        x[1, "n_vacc_status"] <- 100
  
        # calculate vaccination coverage among cases
        x$vcc <- x$n_vaccinated / x$n_vacc_status
        x[which(is.na(x$vcc)), "vcc"] <- NA
        tsw[which(tsw$hz == i), "vcc"] <- x$vcc
  
        # smooth VCC by fitting a binomial additive model
        x$n_unvaccinated <- x$n_vacc_status - x$n_vaccinated
        x <- x[, c("n_vaccinated", "n_unvaccinated", "n_vacc_status", "week")]
        form <- as.formula(paste0("cbind(", "n_vaccinated", ",",
          "n_unvaccinated",")", "~", "pbm(week)"))
        fit <- gamlss(form, family = BI, data = na.omit(x),
          weights = n_vacc_status)
        tsw[which(tsw$hz == i), "vcc_sm"] <- inv.logit(predict(fit))
      }
  
      # Estimate vaccination coverage in population based on smoothed VCC
      tsw$vcp <- tsw$vcc_sm / (1 - (1 - tsw$vcc_sm) * ve)
  
      # Visualise weekly vaccination coverage over time, by HZ
      plot <- ggplot(data = tsw, aes(x = date)) +
        geom_point(aes(y = vcc, size = n_vacc_status), alpha = 0.25,
          colour = palette_cb[6]) +
        geom_line(aes(y = vcc_sm), linewidth = 1, colour = palette_cb[6]) +
        geom_line(aes(y = vcp), linetype = "11", linewidth = 1,
          colour = palette_cb[12]) +
        theme_bw() +
        facet_wrap(~hz, ncol = 5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
        scale_x_date("", minor_breaks = NULL,
          date_breaks = "2 months", date_labels = "%b-%Y" ) +
        scale_y_continuous(labels = percent) +
        scale_size_continuous("number of cases with known vaccination status")
      ggsave(paste(dir_output, "/out_vaxcov_by_hz.pdf", sep = ""),
        height = 23, width = 22, units = "cm", dpi = "print")
  
      # Compute VCP by time window and lag of interest
      for (i in windows_other) {
        for (j in lags) {
          tsw[, paste0("vcp_", i, "w_lag", j)] <-
            f_window_other(var_f = "vcp", window_size = i, lag_f = j,
              op_f = "mean")
        }
      }
    }
    
  #...................................   
  ## Visualise weekly vaccination coverage over time, by province
    
    # Prepare data
    df <- aggregate(tsw[, c("cases", "n_vaccinated", "n_vacc_status")], 
      by = tsw[, c("province", "date")], FUN = sum)
    df <- df[order(df$province, df$date), ]
    
    # Compute VCC
    df$vcc <- df$n_vaccinated / df$n_vacc_status
    
    # Smooth VCC by fitting a binomial growth model, by province
    df$vcc_sm <- NA
    df$vcc_sm_lci <- NA
    df$vcc_sm_uci <- NA
    
    for (i in unique(df$province)) {
      x <- subset(df, province == i)
      x$week <- 1:nrow(x)
      form <- as.formula(paste0("cbind(", "n_vaccinated", ",", 
        "n_vacc_status",")", "~", "s(week, bs = 'cr', pc = 0)"))
      
      # smooth VCC by fitting a binomial additive model
      x$n_unvaccinated <- x$n_vacc_status - x$n_vaccinated
      x <- x[, c("n_vaccinated", "n_unvaccinated", "n_vacc_status", "week")]
      form <- as.formula(paste0("cbind(", "n_vaccinated", ",", 
        "n_unvaccinated",")", "~", "pbm(week)"))
      x$wt <- c(10, rep(1, nrow(x) - 1))
      fit <- gamlss(form, family = BI, data = na.omit(x), 
        weights = wt) 
      df[which(df$province == i), "vcc_sm"] <- inv.logit(predict(fit))
      x <- predict(fit, type = "link", se.fit = T)
      df[which(df$province == i), "vcc_sm_lci"] <- 
        inv.logit(x$fit - 1.96 * x$se.fit)
      df[which(df$province == i), "vcc_sm_uci"] <- 
        inv.logit(x$fit + 1.96 * x$se.fit)
    }

    # Estimate vaccination coverage in population based on smoothed VCC
    df$vcp <- df$vcc_sm / (1 - (1 - df$vcc_sm) * ve)
    df$vcp_lci <- df$vcc_sm_lci / (1 - (1 - df$vcc_sm_lci) * ve)
    df$vcp_uci <- df$vcc_sm_uci / (1 - (1 - df$vcc_sm_uci) * ve)

    # Prepare for plotting
    scale_factor <- max(df$cases, na.rm = T) / max(df$vcp, na.rm = T)

    # Plot
    plot <- ggplot(data = df, aes(x = date)) + 
      geom_bar(aes(y = cases), stat = "identity", fill = palette_cb[2], 
        colour = palette_cb[2], alpha = 0.3) +
      geom_line(aes(y = vcc_sm * scale_factor), colour = palette_cb[6],
        linetype = "solid") +
      geom_ribbon(aes(ymin = vcc_sm_lci*scale_factor, 
        ymax = vcc_sm_uci*scale_factor),
        alpha = 0.2, fill = palette_cb[6]) +
      geom_line(aes(y = vcp * scale_factor), colour = palette_cb[12],
        linetype = "22") +
      geom_ribbon(aes(ymin = vcp_lci*scale_factor, ymax = vcp_uci*scale_factor),
        alpha = 0.2, fill = palette_cb[12]) +
      scale_y_continuous("new EVD cases per week", limits = c(0, NA),
        sec.axis = sec_axis(~./scale_factor, labels = label_percent(),
          name = "estimated EVD vaccination coverage (%)")) +
      scale_x_date("", expand=c(0,0) , minor_breaks=NULL, date_breaks="1 month",
        date_labels = "%b-%Y") +
      theme_bw() +
      facet_wrap(~province, nrow = 2) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Save plot
    ggsave(paste(dir_output, "/out_vaxcov_by_province.pdf", sep = ""), 
      height = 15, width = 20, units = "cm", dpi = "print")

    
#...............................................................................                            
### Preparing explanatory variables based on other datasets
#...............................................................................
    
  #...................................   
  ## Attacks against EVD response
    # Set all NA case totals to zero
    tsw[which(is.na(tsw$against_evd)), "against_evd"] <- 0
    tsw[which(is.na(tsw$suspend_activities)), "suspend_activities"] <- 0
    
    # Number of attacks and suspensions, by time window and lag of interest
    for (i in windows_other) {
      
      for (j in lags) {
        
        tsw[, paste("n_against_evd_", i, "w_lag", j, sep = "")] <- 
          f_window_other(var_f = "against_evd", window_size = i, lag_f = j)
        
        tsw[, paste("n_suspensions_", i, "w_lag", j, sep = "")] <- 
          f_window_other(var_f = "suspend_activities", window_size = i,lag_f= j)
      }
    }

  #...................................   
  ## Insecurity events
    # Set all NA case totals to zero
    tsw[which(is.na(tsw$n_events)), "n_events"] <- 0
    tsw[which(is.na(tsw$fatalities)), "fatalities"] <- 0

    # Number of events and fatalities (as a rate per 100,000 population), 
      # by time window and lag of interest
    for (i in windows_other) {
      
      for (j in lags) {
        tsw[, paste("rate_events_", i, "w_lag", j, sep = "")] <- 
          f_window_other(var_f = "n_events", window_size = i, lag_f = j) * 
          100000 / tsw[, "pop"]
        
        tsw[, paste("rate_fatalities_", i, "w_lag", j, sep = "")] <- 
          f_window_other(var_f = "fatalities", window_size = i, lag_f = j) * 
          100000 / tsw[, "pop"]
      }
    }
    

  #...................................   
  ## Presence of an ETC/transit centre
    # Whether an ETC or TC was open within the health zone at the given time
    tsw <- tsw[order(tsw[, "hz"], tsw[, "epi_year"], tsw[, "epi_week"]), ]
    out <- c()
    for (i in hz_ok) {
      
      # select ETC open column for health zone time series
      tsw_i <- subset(tsw, hz == i)[, "open"]
      
      # base is that an ETC/TC was not open
      out_i <- rep(0, length(tsw_i))
      
      # if an ETC/TC was ever open, set value from that time point to the end of
      # the time series as 'yes'
      if ("yes" %in% tsw_i) {out_i[match("yes", tsw_i):length(out_i)] <- 1}
      
      # add results
      out <- c(out, out_i)
    }
    # add results
    tsw[, "etc_open"] <- out # where 0 = no and 1 = yes
      
    # Whether an ETC or TC was open within the health zone during most of the 
    # time window, by time window and lag of interest
    for (i in windows_other) {
      for (j in lags) {
        tsw[, paste("etc_open_", i, "w_lag", j, sep="")] <- 
          f_window_other(var_f = "etc_open", window_size = i, lag_f = j, 
            op_f = "mean")
        
        tsw[, paste("etc_open_", i, "w_lag", j, "_cat", sep="")] <- 
          ifelse (tsw[, paste("etc_open_", i, "w_lag", j, sep="")] >= 0.50, 
            "yes", "no")
      }
    }      
      
        
  #...................................   
  ## Convert other explanatory variables to population rates per 100,000  
    # Number of health facilities
    tsw$rate_hf <- tsw$n_hf * 100000 / tsw$pop  

    # Road coverage
    tsw$rate_road_length <- tsw$road_length * 100000 / tsw$pop  

    # Number of mines
    tsw$rate_mines <- tsw$n_mines * 100000 / tsw$pop  

        
#.........................................................................................                            
### Exploring variable distributions and creating categorical variables
#.........................................................................................                            

  #...................................   
  ## Explore SDB-related variable distributions and create categorical variables 

    # Proportion of successful SDBs
    f_hist("p_success_prev", tsw, c(NA, NA))
    for (i in grep("p_success", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.24999, 0.49999, 0.74999, 1),
        labels = c("< 0.25", "0.25 to 0.49", "0.50 to 0.74", ">= 0.75"), 
        include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }
    
    # Proportion of timely SDBs
    f_hist("p_timely_prev", tsw, c(NA, NA))
    for (i in grep("p_timely", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.99999, 1),
        labels = c("< 1.00", "1.00"), include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }

    # Ratio of SDBs (all) to population
    f_hist("r_sdb_to_pop_prev", tsw, c(NA,NA))
    for (i in grep("_to_pop", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.0001, 1.999, 3.999, 5.999, 100),
        labels = c("0", "0.01 to 1.99", "2.00 to 3.99", 
          "4.00 to 5.99", ">= 6.00"), include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }    
    
    
  #...................................   
  ## Explore confounder variable distributions and create categorical variables 
    
    # Percentage of cases with nosocomial transmission
    f_hist("p_nosocomial_prev", tsw, c(NA, NA))
    for (i in grep("p_nosocomial", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.00001, 1), 
        labels = c("0%", "> 0%"), include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }

    # Percentage of cases with a known link to another case
    f_hist("p_epi_link_prev", tsw, c(NA, NA))
    for (i in grep("p_epi_link", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.333, 0.666, 1), 
        labels = c("<33%", "33 to 66%", "> 66%"), include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }

    # Vaccination coverage (population)
    f_hist("vcp_2w_lag2", tsw, c(NA, NA))
    for (i in grep("vcp", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.249, 0.499, 0.749, 1),
        labels = c("< 25%", "25 to 49%", "50 to 74%", ">= 75%"), 
        include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }
    
    # Number of attacks against EVD
    f_hist("n_against_evd_2w_lag2", tsw, c(NA, NA))
    for (i in grep("n_against_evd", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.999, 1.999, 2.999, 100),
        labels = c("0", "1", "2", ">= 3"), include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }

    # Number of EVD response suspensions
    f_hist("n_suspensions_2w_lag2", tsw, c(NA, NA))
    for (i in grep("n_suspensions", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.999, 1.999, 100),
        labels = c("0", "1", ">= 2"), include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }
    
    # Rate of insecurity events
    f_hist("rate_events_2w_lag2", tsw, c(NA, NA))
    for (i in grep("rate_events", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.001, 4.999, 100),
        labels = c("0", "0.1 to 4.9", ">= 5.0"), include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }

    # Rate of insecurity deaths
    f_hist("rate_fatalities_2w_lag2", tsw, c(NA, NA))
    for (i in grep("rate_fatalities", colnames(tsw), value = T)) {
      tsw[, paste0(i, "_cat")] <- cut(tsw[, i], 
        breaks = c(0, 0.001, 9.999, 100),
        labels = c("0", "0.1 to 9.9", ">= 10.0"), include.lowest = T)
      print(paste0("table of ", i, "_cat") )
      print(table(tsw[, paste0(i, "_cat")] ) )
    }    

    # Number of networks
    f_hist("n_networks", tsw, c(NA, NA))
    tsw[, "n_networks_cat"] <- cut(tsw[, "n_networks"], 
      breaks = c(0, 3.999, 100),
      labels = c("< 4 networks", "4 networks"), include.lowest = T)
    table(tsw$n_networks_cat)
    
    # Number of radio frequencies
    f_hist("n_frequencies", tsw, c(NA, NA))
    tsw[, "n_frequencies_cat"] <- cut(tsw[, "n_frequencies"], 
      breaks = c(0, 9.999, 19.999, 100),
      labels = c("< 10 frequencies", "10-19 frequencies", ">= 20 frequencies"), 
      include.lowest = T)
    table(tsw$n_frequencies_cat)
  
    # Number of health facilities per population
    f_hist("rate_hf", tsw, c(NA, NA))
    tsw[, "rate_hf_cat"] <- cut(tsw[, "rate_hf"], 
      breaks = c(0, 24.999, 49.999, 100),
      labels = c("< 25.0", "25.0 to 49.9", ">= 50.0"), include.lowest = T)
    table(tsw$rate_hf_cat)

    # Road length per population
    f_hist("rate_road_length", tsw, c(NA, NA))
    tsw[, "rate_road_length_cat"] <- cut(tsw[, "rate_road_length"], 
      breaks = c(0, 199.99, 399.99, 1000),
      labels = c("< 200 Km", "200-399 Km", ">= 400 Km"), include.lowest = T)
    table(tsw$rate_road_length_cat)
    
    # Number of mines per population
    f_hist("rate_mines", tsw, c(NA, NA))
    tsw[, "rate_mines_cat"] <- cut(tsw[, "rate_mines"], 
      breaks = c(0, 0.001, 1000),
      labels = c("no mining", "some mining"), include.lowest = T)
    table(tsw$rate_mines_cat)
  

  #...................................   
  ## Explore dependent variables' distributions and transform if appropriate
    
    # Explore net reproduction number distribution
      f_hist("r_mean", tsw, c(NA, NA))
      # distribution is skewed, may be better to log it
      for (i in grep("r_mean", colnames(tsw), value = T)) {
        tsw[, paste("ln_", i, sep="")] <- log(tsw[, i])
        }
      f_hist("ln_r_mean", tsw, c(NA, NA))

    # Explore relative period incidence
      # generate proportion of cases in current window, 
      # out of cases in current + previous window (just to visualise)
      tsw$p_curr <- tsw$cases_curr / (tsw$cases_curr + tsw$cases_prev)

      # clean up Na, NaN, Inf
      tsw$p_curr <- ifelse(tsw$p_curr %in% c(NA, NaN, Inf, -Inf), NA,tsw$p_curr) 

      # visualise distributions
      f_hist("p_curr", tsw, c(NA, NA))
    
      
  #...................................   
  ## Explore variable completeness
    
    # List of variables
    vars <- c("cases_curr", "p_curr", "r_mean", 
      "p_success_curr", "p_success_prev", "p_timely_curr", "p_timely_prev", 
      "r_sdb_to_pop_curr", "r_sdb_to_pop_prev", 
      "p_nosocomial_prev","p_nosocomial_curr", 
      "p_epi_link_prev", "p_epi_link_curr",
      "vcp_2w_lag0", "n_against_evd_2w_lag0", "n_suspensions_2w_lag0", 
      "rate_events_2w_lag0", "rate_fatalities_2w_lag0", "etc_open_2w_lag0",
      "rate_hf", "rate_road_length", "rate_mines")
    
    # Prepare completeness-over-time dataset: 0 = missing; 1 = nonmissing
      df <- tsw[, c("date", vars)]
      df[, vars] <- ifelse(is.na(df[, vars]), 0, 1)
      df$n_obs <- 1
    
      # restrict to observations when the outcome is non-missing
      df <- subset(df, p_curr == 1)
      
      # aggregate
      df <- aggregate(df[, c(vars, "n_obs")], by = list("date" = df$date), 
        FUN = sum, na.rm = T)  
      df[, vars] <- df[, vars] / df[, "n_obs"]
      
      # melt
      df <- reshape2::melt(df, id.vars = "date")
      df <- subset(df, variable != "n_obs")  

    # Plot completeness
    plot <- ggplot(data = df, aes(x = date, y = variable, fill = value)) +
      geom_raster()
    plot

  #...................................   
  ## Clean up any NaN and Inf values
    tsw[is.na(tsw)] <- NA  
    tsw[tsw == Inf] <- NA
    tsw[tsw == -Inf] <- NA 
      

  #...................................
  ## Turn categorical variables into factors

    # Data frame to denote which variables are independent variables to include, 
      # and which should be factors
      vars <- data.frame(variable = colnames(tsw), 
        independent = rep("no", length(colnames(tsw) ) ), 
        is_factor = rep("no", length(colnames(tsw) )  ) )
      x <- grep("_cat", colnames(tsw), value = T)
      vars[vars$variable %in% x, c("independent", "is_factor")] <- "yes"
    
    # Transform categorical variables into factors  
    for (i in 1:nrow(vars) ) {
        x <- vars[i, "variable"]
        if (vars[i, "is_factor"] == "yes") {tsw[, x] <- as.factor(tsw[, x]) }
     }
      
    # See breakdown of categories
    for (i in vars[vars$independent == "yes", "variable"]) 
      {print(i); print(table(tsw[, i]))}

    # Specify reference categories for categorical variables
    vars[, "ref"] <- NA
    vars[vars$variable %in% 
        grep("etc_open", vars$variable, value = T), "ref"] <- "no"
    vars[vars$variable %in% 
        grep("p_success", vars$variable, value = T), "ref"] <- "< 0.25"         
    vars[vars$variable %in% 
        grep("p_timely", vars$variable, value = T), "ref"] <- "< 1.00"
    vars[vars$variable %in% 
        grep("r_sdb_to_pop", vars$variable, value = T), "ref"] <- "0"
    vars[vars$variable %in% 
        grep("p_nosocomial", vars$variable, value = T), "ref"] <- "0%"
    vars[vars$variable %in% 
        grep("p_epi_link", vars$variable, value = T), "ref"] <- "> 66%"
    vars[vars$variable %in% 
        grep("vcp", vars$variable, value = T), "ref"] <- "< 25%"
    vars[vars$variable %in% 
        grep("n_against_evd", vars$variable, value = T), "ref"] <- "0"      
    vars[vars$variable %in% 
        grep("n_suspensions", vars$variable, value = T), "ref"] <- "0"
    vars[vars$variable %in% 
        grep("rate_events", vars$variable, value = T), "ref"] <- "0"
    vars[vars$variable %in% 
        grep("rate_fatalities", vars$variable, value = T), "ref"] <- "0"      
    vars[vars$variable %in% 
        grep("n_networks", vars$variable, value = T), "ref"] <- 
      "< 4 networks"
    vars[vars$variable %in% 
        grep("n_frequencies", vars$variable, value = T), "ref"] <- 
      "< 10 frequencies"
    vars[vars$variable %in% 
        grep("rate_hf", vars$variable, value = T), "ref"] <- ">= 50.0"
    vars[vars$variable %in% 
        grep("rate_road_length", vars$variable, value = T), "ref"] <- 
      ">= 400 Km"
    vars[vars$variable %in% 
        grep("rate_mines", vars$variable, value = T), "ref"] <- "no mining"

    # Relevel
    for (i in 1:nrow(vars) ) {
      x <- vars[i, "variable"]
      if (vars[i, "is_factor"] == "yes") {tsw[, x] <- 
        relevel(tsw[, x], vars[i, "ref"]) }
    }


  #...................................
  ## Save dataset
  write.csv(tsw, paste0(dir_output, "/out_tsw.csv"), row.names = F)
             
#.........................................................................................
### ENDS
#.........................................................................................


