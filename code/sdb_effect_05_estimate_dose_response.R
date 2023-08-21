##..............................................................................
### +++ EFFECT OF SAFE AND DIGNIFIED BURIALS ON EBOV TRANSMISSION IN DRC +++ ###
#...............................................................................

#...............................................................................
## -- R SCRIPT TO ESTIMATE DOSE-RESPONSE RELATIONSHIPS B/W SDB AND OUTCOME -- ##
#...............................................................................

                            # Written by Francesco Checchi, LSHTM (March 2023)
                            # francesco.checchi@lshtm.ac.uk 


  # based on Hirano & Imbens (2004), 
    # https://www.math.mcgill.ca/dstephens/PSMMA/Articles/HIrano-Imbens-2004.pdf
  # and Austin (2017), 
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5969262/#sim7615-bib-0003 

#...............................................................................
### Computing dose-response associations between SDB success and incidence
#...............................................................................
    
  #...................................
  ## Hirano & Imbens method
    
    # Specify levels of the exposure that function is evaluated at    
    exp_levels <- seq(0, 1, by = 0.05)    
    
    # Prepare prediction dataset
    df <- tsw_inc

    # Fit GPS exposure model again
    fit_exp_inc <- f_gps_lite(tsw_inc, exposures[1], confounders_inc)
    
    # Prepare output dose-response dataset
    out <- data.frame(exp_levels, out_mean = NA, out_lci = NA, out_uci = NA)
    colnames(out)[1] <- exposures[1]
    
    # For each exposure level...
    for (i in exp_levels) {
      
      # update prediction dataset
      df[, exposures[1]] <- i

      # predict GPS at exposure level
      df$gps_p_success <- dnorm(x = df[, exposures[1]], 
        mean = fitted.values(fit_exp_inc), sd = summary(fit_exp_inc)$sigma)

      # predict outcome at exposure level
      df[, "out_mean"] <- predict(fit_inc_p_success_hi, newdata = df, 
        type = "response")
      x <- as.data.frame(predict(fit_inc_p_success_hi, newdata = df, 
        se.fit = TRUE))
      df[, "out_lci"] <- inv.logit(x[, 1] - 1.96 * x[, 2])
      df[, "out_uci"] <- inv.logit(x[, 1] + 1.96 * x[, 2])
      
      # mean and 95%CI of outcome for this exposure level
      out[out[, exposures[1]] == i, "out_mean"] <- mean(df[, "out_mean"])
      out[out[, exposures[1]] == i, "out_lci"] <- mean(df[, "out_lci"])
      out[out[, exposures[1]] == i, "out_uci"] <- mean(df[, "out_uci"])
      
    }
 
    # Plot dose-response function
    plot_inc_hi <- ggplot(out, aes(x = eval(as.name(exposures[1])))) +
      geom_line(aes(y = out_mean), colour = palette_cb[6]) +
      geom_ribbon(aes(ymin = out_lci, ymax = out_uci), alpha = 0.3, 
        fill = palette_cb[6]) +
      theme_bw() +
      theme(plot.margin = ggplot2::margin(1, 1, 0.25, 0.25, "cm")) +
      scale_y_continuous(name = expression(pi[t]), expand = c(0, 0),
        limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
      scale_x_continuous("percentage of successful SDBs", expand = c(0, 0),
        limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.1) ) +
      geom_hline(yintercept = 0.5, lty = "21", colour = palette_cb[7], 
        linewidth = 1) +
      geom_segment(x = 0.1, xend = 0.1, y = 0, yend = 0.5, 
        arrow = arrow(ends = "both", type = "open", length = unit(0.3, "cm")), 
        colour = palette_cb[7]) +
      annotate("text", x = 0.27, y = 0.25, label = "decreasing transmission", 
        colour = palette_cb[7])
    
    # Save plot and estimates
    ggsave(paste(dir_output, "/out_dose_resp_inc_prop_success_hi.png", sep =""),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out, paste(dir_output, 
      "/out_dose_resp_inc_prop_success_hi.csv", sep = ""), row.names = FALSE)    
    
    
  #...................................
  ## Propensity weights method
    
    # Specify levels of the exposures[1] that function is evaluated at    
    exp_levels <- seq(0, 1, by = 0.05)    
    
    # Prepare prediction dataset
    df <- tsw_inc

    # Fit GPS exposure model again
    fit_exp_inc <- f_gps_lite(tsw_inc, exposures[1], confounders_inc)
    
    # Prepare output dose-response dataset
    out <- data.frame(exp_levels, out_mean = NA, out_lci = NA, out_uci = NA)
    colnames(out)[1] <- exposures[1]
    
    # For each exposure level...
    for (i in exp_levels) {
      
      # update prediction dataset
      df[, exposures[1]] <- i

      # predict outcome at exposure level
      df[, "out_mean"] <- predict(fit_inc_p_success_pw, newdata = df, 
        type = "response")
      x <- as.data.frame(predict(fit_inc_p_success_pw, newdata = df, 
        se.fit = TRUE))
      df[, "out_lci"] <- inv.logit(x[, 1] - 1.96 * x[, 2])
      df[, "out_uci"] <- inv.logit(x[, 1] + 1.96 * x[, 2])
      
      # mean and 95%CI of outcome for this exposure level
      out[out[, exposures[1]] == i, "out_mean"] <- mean(df[, "out_mean"])
      out[out[, exposures[1]] == i, "out_lci"] <- mean(df[, "out_lci"])
      out[out[, exposures[1]] == i, "out_uci"] <- mean(df[, "out_uci"])
      
    }
 
    # Plot dose-response function
    plot_inc_pw <- ggplot(out, aes(x = eval(as.name(exposures[1])))) +
      geom_line(aes(y = out_mean), colour = palette_cb[6]) +
      geom_ribbon(aes(ymin = out_lci, ymax = out_uci), alpha = 0.3, 
        fill = palette_cb[6]) +
      theme_bw() +
      theme(plot.margin = ggplot2::margin(1, 1, 0.25, 0.25, "cm")) +
      scale_y_continuous(name = expression(pi[t]), expand = c(0, 0),
        limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
      scale_x_continuous("percentage of successful SDBs", expand = c(0, 0),
        limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.1) ) +
      geom_hline(yintercept = 0.5, lty = "21", colour = palette_cb[7], 
        linewidth = 1) +
      geom_segment(x = 0.1, xend = 0.1, y = 0, yend = 0.5, 
        arrow = arrow(ends = "both", type = "open", length = unit(0.3, "cm")), 
        colour = palette_cb[7]) +
      annotate("text", x = 0.27, y = 0.25, label = "decreasing transmission", 
        colour = palette_cb[7])
    
    # Save plot and estimates
    ggsave(paste(dir_output, "/out_dose_resp_inc_prop_success_pw.png", sep =""),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out, paste(dir_output, 
      "/out_dose_resp_inc_prop_success_pw.csv", sep = ""), row.names = FALSE)    
    
    
#...............................................................................
### Computing dose-response associations between SDB success and reproduction no
#...............................................................................

  #...................................
  ## Hirano & Imbens method
    
    # Specify levels of the exposure that function is evaluated at    
    exp_levels <- seq(0, 1, by = 0.05)    
    
    # Prepare prediction dataset
    df <- tsw_rn

    # Fit GPS exposure model again
    fit_exp_rn <- f_gps_lite(tsw_rn, exposures[4], confounders_rn)
    
    # Prepare output dose-response dataset
    out <- data.frame(exp_levels, out_mean = NA, out_lci = NA, out_uci = NA)
    colnames(out)[1] <- exposures[4]
    
    # For each exposure level...
    for (i in exp_levels) {
      
      # update prediction dataset
      df[, exposures[4]] <- i

      # predict GPS at exposure level
      df$gps_p_success <- dnorm(x = df[, exposures[4]], 
        mean = fitted.values(fit_exp_rn), sd = summary(fit_exp_rn)$sigma)

      # predict outcome at exposure level
      df[, "out_mean"] <- predict(fit_rn_p_success_hi, newdata = df, 
        type = "response")
      x <- as.data.frame(predict(fit_rn_p_success_hi, newdata = df, 
        se.fit = TRUE))
      df[, "out_lci"] <- x[, 1] - 1.96 * x[, 2]
      df[, "out_uci"] <- x[, 1] + 1.96 * x[, 2]
      
      # mean and 95%CI of outcome for this exposure level
      out[out[, exposures[4]] == i, "out_mean"] <- mean(df[, "out_mean"])
      out[out[, exposures[4]] == i, "out_lci"] <- mean(df[, "out_lci"])
      out[out[, exposures[4]] == i, "out_uci"] <- mean(df[, "out_uci"])
      
    }
 
    # Plot dose-response function
    plot_rn_hi <- ggplot(out, aes(x = eval(as.name(exposures[4])))) +
      geom_line(aes(y = out_mean), colour = palette_cb[4]) +
      geom_ribbon(aes(ymin = out_lci, ymax = out_uci), alpha = 0.3, 
        fill = palette_cb[4]) +
      theme_bw() +
      theme(plot.margin = ggplot2::margin(1, 1, 0.25, 0.25, "cm")) +
      scale_y_continuous("effective reproduction number", expand = c(0, 0),
        limits = c(0, 2.1), breaks = seq(0, 2, by = 0.2) ) +
      scale_x_continuous("percentage of successful SDBs", expand = c(0, 0),
        limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.1) ) +
      geom_hline(yintercept = 1, lty = "21", colour = palette_cb[7], 
        linewidth = 1) +
      geom_segment(x = 0.1, xend = 0.1, y = 0, yend = 1, 
        arrow = arrow(ends = "both", type = "open", length = unit(0.3, "cm")), 
        colour = palette_cb[7]) +
      annotate("text", x = 0.27, y = 0.5, label = "decreasing transmission", 
        colour = palette_cb[7])
    
    # Save plot and estimates
    ggsave(paste(dir_output, "/out_dose_resp_rn_prop_success_hi.png", sep =""),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out, paste(dir_output, 
      "/out_dose_resp_rn_prop_success_hi.csv", sep = ""), row.names = FALSE)    
    
    
  #...................................
  ## Propensity weights method
    
    # Specify levels of the exposure that function is evaluated at    
    exp_levels <- seq(0, 1, by = 0.05)    
    
    # Prepare prediction dataset
    df <- tsw_rn

    # Fit GPS exposure model again
    fit_exp_rn <- f_gps_lite(tsw_rn, exposures[4], confounders_rn)
    
    # Prepare output dose-response dataset
    out <- data.frame(exp_levels, out_mean = NA, out_lci = NA, out_uci = NA)
    colnames(out)[1] <- exposures[4]
    
    # For each exposure level...
    for (i in exp_levels) {
      
      # update prediction dataset
      df[, exposures[4]] <- i

      # predict outcome at exposure level
      df[, "out_mean"] <- predict(fit_rn_p_success_pw, newdata = df, 
        type = "response")
      x <- as.data.frame(predict(fit_rn_p_success_pw, newdata = df, 
        se.fit = TRUE))
      df[, "out_lci"] <- x[, 1] - 1.96 * x[, 2]
      df[, "out_uci"] <- x[, 1] + 1.96 * x[, 2]
      
      # mean and 95%CI of outcome for this exposure level
      out[out[, exposures[4]] == i, "out_mean"] <- mean(df[, "out_mean"])
      out[out[, exposures[4]] == i, "out_lci"] <- mean(df[, "out_lci"])
      out[out[, exposures[4]] == i, "out_uci"] <- mean(df[, "out_uci"])
      
    }
 
    # Plot dose-response function
    plot_rn_pw <- ggplot(out, aes(x = eval(as.name(exposures[4])))) +
      geom_line(aes(y = out_mean), colour = palette_cb[4]) +
      geom_ribbon(aes(ymin = out_lci, ymax = out_uci), alpha = 0.3, 
        fill = palette_cb[4]) +
      theme_bw() +
      theme(plot.margin = ggplot2::margin(1, 1, 0.25, 0.25, "cm")) +
      scale_y_continuous("effective reproduction number", expand = c(0, 0),
        limits = c(0, 2.1), breaks = seq(0, 2, by = 0.2) ) +
      scale_x_continuous("percentage of successful SDBs", expand = c(0, 0),
        limits = c(0, 1), labels = scales::percent, breaks = seq(0, 1, 0.1) ) +
      geom_hline(yintercept = 1, lty = "21", colour = palette_cb[7], 
        linewidth = 1) +
      geom_segment(x = 0.1, xend = 0.1, y = 0, yend = 1, 
        arrow = arrow(ends = "both", type = "open", length = unit(0.3, "cm")), 
        colour = palette_cb[7]) +
      annotate("text", x = 0.27, y = 0.5, label = "decreasing transmission", 
        colour = palette_cb[7])
    
    # Save plot and estimates
    ggsave(paste(dir_output, "/out_dose_resp_rn_prop_success_pw.png", sep =""),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out, paste(dir_output, 
      "/out_dose_resp_rn_prop_success_pw.csv", sep = ""), row.names = FALSE)    
    
 
#...............................................................................
### Preparing a combination graph
#...............................................................................
    
    # Prepare graph
    plot <- ggarrange(plot_inc_hi, plot_rn_hi, plot_inc_pw, plot_rn_pw,
      ncol = 2, nrow = 2, labels = c(
        "Effect on next-window incidence (Hirano-Imbens method)",
        "Effect on reproduction number (Hirano-Imbens method)",
        "Effect on next-window incidence (propensity weights method)",
        "Effect on reproduction number (propensity weights method)"),
      label.x = 0, hjust = 0, font.label = list(size = 11, colour = "grey20"),
      align = "hv")
    
    # Save plot
    ggsave(paste(dir_output, "/out_dose_resp_combi.png", sep =""),
      dpi = "print", width = 30, height = 30, units = "cm")
    

             
#...............................................................................
### ENDS
#...............................................................................

