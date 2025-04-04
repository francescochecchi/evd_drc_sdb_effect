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
    
    # Compute dose-response function
    out <- f_dose(fit_f = fit_inc_p_success_hi, data_f = tsw_inc_p_success,
      exposure_f = exposures[1], confounders_f = confounders_inc_p_success,
      outcome_f = "inc", method_adj_f = "hi", method_gps_f = "lm")
    
    # Save plot and estimates
    plot_inc_hi <- out$plot
    ggsave(paste0(dir_output, "/out_dose_resp_inc_p_success_hi.png"),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out$out, paste(dir_output, 
      "/out_dose_resp_inc_p_success_hi.csv", sep = ""), row.names = F)    
    
 
  #...................................
  ## Hirano & Imbens method - with interaction
    
    # Compute dose-response function
    out <- f_dose(fit_f = fit_inc_p_success_hi_int, data_f = tsw_inc_p_success,
      exposure_f = exposures[1], confounders_f = confounders_inc_p_success,
      outcome_f = "inc", method_adj_f = "hi")
    
    # Save plot and estimates
    plot_inc_hi_int <- out$plot
    ggsave(paste0(dir_output, "/out_dose_resp_inc_p_success_hi_int.png"),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out$out, paste(dir_output, 
      "/out_dose_resp_inc_p_success_hi_int.csv", sep = ""), row.names = F)    
    
       
  #...................................
  ## Propensity weights method
    
    # Compute dose-response function
    out <- f_dose(fit_f = fit_inc_p_success_pw, data_f = tsw_inc_p_success,
      exposure_f = exposures[1], confounders_f = confounders_inc_p_success,
      outcome_f = "inc", method_adj_f = "pw")
    
    # Save plot and estimates
    plot_inc_pw <- out$plot
    ggsave(paste0(dir_output, "/out_dose_resp_inc_p_success_pw.png"),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out$out, paste(dir_output, 
      "/out_dose_resp_inc_p_success_pw.csv", sep = ""), row.names = F)    
     
  #...................................
  ## Propensity weights method - with interaction
    
    # Compute dose-response function
    out <- f_dose(fit_f = fit_inc_p_success_pw_int, data_f = tsw_inc_p_success,
      exposure_f = exposures[1], confounders_f = confounders_inc_p_success,
      outcome_f = "inc", method_adj_f = "pw")
    
    # Save plot and estimates
    plot_inc_pw_int <- out$plot
    ggsave(paste0(dir_output, "/out_dose_resp_inc_p_success_pw_int.png"),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out$out, paste(dir_output, 
      "/out_dose_resp_inc_p_success_pw_int.csv", sep = ""), row.names = F)    
   
      
#...............................................................................
### Computing dose-response associations between SDB success and reproduction no
#...............................................................................

  #...................................
  ## Hirano & Imbens method
    
    # Compute dose-response function
    out <- f_dose(fit_f = fit_rn_p_success_hi, data_f = tsw_rn_p_success,
      exposure_f = exposures[4], confounders_f = confounders_rn_p_success,
      outcome_f = "rn", method_adj_f = "hi")
    
    # Save plot and estimates
    plot_rn_hi <- out$plot
    ggsave(paste0(dir_output, "/out_dose_resp_rn_p_success_hi.png"),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out$out, paste(dir_output, 
      "/out_dose_resp_rn_p_success_hi.csv", sep = ""), row.names = F)    
   
  #...................................
  ## Hirano & Imbens method - with interaction
    
    # Compute dose-response function
    out <- f_dose(fit_f = fit_rn_p_success_hi_int, data_f = tsw_rn_p_success,
      exposure_f = exposures[4], confounders_f = confounders_rn_p_success,
      outcome_f = "rn", method_adj_f = "hi")
    
    # Save plot and estimates
    plot_rn_hi_int <- out$plot
    ggsave(paste0(dir_output, "/out_dose_resp_rn_p_success_hi_int.png"),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out$out, paste(dir_output, 
      "/out_dose_resp_rn_p_success_hi_int.csv", sep = ""), row.names = F)    
    

  #...................................
  ## Propensity weights method
    
    # Compute dose-response function
    out <- f_dose(fit_f = fit_rn_p_success_pw, data_f = tsw_rn_p_success,
      exposure_f = exposures[4], confounders_f = confounders_rn_p_success,
      outcome_f = "rn", method_adj_f = "pw")
    
    # Save plot and estimates
    plot_rn_pw <- out$plot
    ggsave(paste0(dir_output, "/out_dose_resp_rn_p_success_pw.png"),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out$out, paste(dir_output, 
      "/out_dose_resp_rn_p_success_pw.csv", sep = ""), row.names = F)    
 
  #...................................
  ## Propensity weights method - with interaction
    
    # Compute dose-response function
    out <- f_dose(fit_f = fit_rn_p_success_pw_int, data_f = tsw_rn_p_success,
      exposure_f = exposures[4], confounders_f = confounders_rn_p_success,
      outcome_f = "rn", method_adj_f = "pw")
    
    # Save plot and estimates
    plot_rn_pw_int <- out$plot
    ggsave(paste0(dir_output, "/out_dose_resp_rn_p_success_pw_int.png"),
      dpi = "print", width = 15, height = 10, units = "cm")
    
    write.csv(out$out, paste(dir_output, 
      "/out_dose_resp_rn_p_success_pw_int.csv", sep = ""), row.names = F)    
 
   
     
#...............................................................................
### Preparing combination graphs
#...............................................................................

  #...................................
  ## Without interaction
        
    # Prepare graph
    plot <- ggarrange(plot_inc_hi, plot_rn_hi, plot_inc_pw, plot_rn_pw,
      ncol = 2, nrow = 2, labels = c(
        "Effect on next-window incidence (Hirano-Imbens method)",
        "Effect on net reproduction number (Hirano-Imbens method)",
        "Effect on next-window incidence (propensity weights method)",
        "Effect on net reproduction number (propensity weights method)"),
      label.x = 0, hjust = 0, font.label = list(size = 11, colour = "grey20"),
      align = "hv")
    
    # Save plot
    ggsave(paste0(dir_output, "/out_dose_resp_combi.png"),
      dpi = "print", width = 30, height = 30, units = "cm")
    
  #...................................
  ## With interaction
        
    # Prepare graph
    plot <- ggarrange(plot_inc_hi_int, plot_rn_hi_int, plot_inc_pw_int, 
      plot_rn_pw_int, ncol = 2, nrow = 2, labels = c(
        "Effect on next-window incidence (Hirano-Imbens method)",
        "Effect on net reproduction number (Hirano-Imbens method)",
        "Effect on next-window incidence (propensity weights method)",
        "Effect on net reproduction number (propensity weights method)"),
      label.x = 0, hjust = 0, font.label = list(size = 11, colour = "grey20"),
      align = "hv")
    
    # Save plot
    ggsave(paste0(dir_output, "/out_dose_resp_combi_int.png"),
      dpi = "print", width = 30, height = 30, units = "cm")
    

             
#...............................................................................
### ENDS
#...............................................................................

