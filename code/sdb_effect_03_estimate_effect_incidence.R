#...............................................................................
### +++ EFFECT OF SAFE AND DIGNIFIED BURIALS ON EBOV TRANSMISSION IN DRC +++ ###
#...............................................................................

#...............................................................................
## -------- R SCRIPT TO ESTIMATE EFFECTS OF SDB ON INCIDENCE CHANGE --------- ##
#...............................................................................

                            # Written by Francesco Checchi, LSHTM (March 2023)
                            # francesco.checchi@lshtm.ac.uk 


#...............................................................................
### Estimating effect of SDB success on incidence change
#...............................................................................

  #...................................
  ## Hirano-Imbens method
      
    # Define formula for outcome as a function of exposure and GPS
    form_curr <- as.formula(paste("cbind(", outcomes[1], ",", outcomes[2],")", 
      "~", 
      paste(c(exposures[1], "gps_p_success", "(1|hz)"), collapse = "+") ) )

    # Observe correlation matrix of confounders
    corrplot(cor(tsw_inc[, c(confounders_inc, exposures[2], "gps_p_success")]), 
      type = "upper")
      
    # Fit GLMM for outcome model
    fit_curr <- glmmTMB(form_curr, data = tsw_inc, family = binomial() )
    model_parameters(fit_curr, exponentiate = TRUE)
    AIC(fit_curr)

    # Adjust for confounders (add one at a time, retain if it improves AIC / 
      # modifies effect size by > 10%)
      # proportion of timely SDBs
      form_new <- update(form_curr, ~ . + p_timely_prev_3w_cat )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new
      
      # proportion of nosocomial cases
      form_new <- update(form_curr, ~ . + p_nosocomial_prev_3w )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

      # proportion of cases with a known epidemiological link
      form_new <- update(form_curr, ~ . + p_epi_link_prev_3w )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

      # Ebola vaccination coverage
      form_new <- update(form_curr, ~ . + vc_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

      # presence of an open Ebola treatment centre
      form_new <- update(form_curr, ~ . + etc_open_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new
      
      # attacks against EVD response
      form_new <- update(form_curr, ~ . + n_against_evd_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # rate of insecurity events
      form_new <- update(form_curr, ~ . + rate_events_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # number of phone networks
      form_new <- update(form_curr, ~ . + n_networks )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
            
      # number of radio frequencies
      form_new <- update(form_curr, ~ . + n_frequencies )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # presence of health facilities
      form_new <- update(form_curr, ~ . + rate_hf )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # road network
      form_new <- update(form_curr, ~ . + rate_road_length )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude      
  
      # ratio of SDBs to population
      form_new <- update(form_curr, ~ . + r_sdb_to_pop_prev_3w )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

    # Explore plausible interactions 
      # exposure with proportion timely
      form_new <- update(form_curr, ~ . + 
          p_success_prev_3w : p_timely_prev_3w_cat)
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude

      # exposure with ratio of SDB burials per population
      form_new <- update(form_curr, ~ . + 
          p_success_prev_3w : r_sdb_to_pop_prev_3w)
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new
     
    # Save model
    fit_inc_p_success_hi <- fit_curr
    x <- model_parameters(fit_curr, exponentiate = TRUE) 
    out <- flextable(x[, c("Parameter", "Coefficient", "CI_low", "CI_high", "p",
      "Effects")])
    out <- set_formatter(out, 
      "Coefficient" = function (x) {sprintf("%.02f", x)},
      "CI_low" = function (x) {sprintf("%.02f", x)},
      "CI_high" = function (x) {sprintf("%.02f", x)},
      "p" = function (x) {sprintf("%.03f", x)}
    )  
    out <- autofit(out)
    flextable::save_as_docx(out, 
      path = paste(dir_output, "/out_glmm_inc_p_success_hi.docx", sep =""))

    
  #...................................
  ## Propensity weighting method
      
    # Define formula for outcome as a function of exposure and GPS
    form_curr <- as.formula(paste("cbind(", outcomes[1], ",", outcomes[2],")", 
      "~", 
      paste(c(exposures[1], "(1|hz)"), collapse = "+") ) )

    # Fit GLMM for outcome model
    fit_curr <- glmmTMB(form_curr, data = tsw_inc, family = binomial(),
      weights = gps_p_success_wt )
    model_parameters(fit_curr, exponentiate = TRUE)
    AIC(fit_curr)

    # Adjust for confounders (add one at a time, retain if it improves AIC / 
      # modifies effect size by > 10%)
      # proportion of timely SDBs
      form_new <- update(form_curr, ~ . + p_timely_prev_3w_cat )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new
      
      # proportion of nosocomial cases
      form_new <- update(form_curr, ~ . + p_nosocomial_prev_3w )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

      # proportion of cases with a known epidemiological link
      form_new <- update(form_curr, ~ . + p_epi_link_prev_3w )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

      # Ebola vaccination coverage
      form_new <- update(form_curr, ~ . + vc_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

      # presence of an open Ebola treatment centre
      form_new <- update(form_curr, ~ . + etc_open_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new
      
      # attacks against EVD response
      form_new <- update(form_curr, ~ . + n_against_evd_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # rate of insecurity events
      form_new <- update(form_curr, ~ . + rate_events_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # number of phone networks
      form_new <- update(form_curr, ~ . + n_networks )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new
            
      # number of radio frequencies
      form_new <- update(form_curr, ~ . + n_frequencies )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # presence of health facilities
      form_new <- update(form_curr, ~ . + rate_hf )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # road network
      form_new <- update(form_curr, ~ . + rate_road_length )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude      
  
      # ratio of SDBs to population
      form_new <- update(form_curr, ~ . + r_sdb_to_pop_prev_3w )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude

    # Explore plausible interactions 
      # exposure with proportion timely
      form_new <- update(form_curr, ~ . + 
          p_success_prev_3w : p_timely_prev_3w_cat)
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude, difficult to interpret coefficients

      # exposure with ratio of SDB burials per population
      form_new <- update(form_curr, ~ . + 
          p_success_prev_3w : r_sdb_to_pop_prev_3w)
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude

    # Save model
    fit_inc_p_success_pw <- fit_curr
    x <- model_parameters(fit_curr, exponentiate = TRUE) 
    out <- flextable(x[, c("Parameter", "Coefficient", "CI_low", "CI_high", "p",
      "Effects")])
    out <- set_formatter(out, 
      "Coefficient" = function (x) {sprintf("%.02f", x)},
      "CI_low" = function (x) {sprintf("%.02f", x)},
      "CI_high" = function (x) {sprintf("%.02f", x)},
      "p" = function (x) {sprintf("%.03f", x)}
    )  
    out <- autofit(out)
    flextable::save_as_docx(out, 
      path = paste(dir_output, "/out_glmm_inc_p_success_pw.docx", sep =""))

    
    
#...............................................................................
### Estimating effect of SDB timeliness on incidence change
#...............................................................................

  #...................................
  ## Weighting method (only option)

    # Define formula for outcome as a function of exposure
    form_curr <- as.formula(paste("cbind(", outcomes[1], ",", outcomes[2],")", 
      "~", 
      paste(c(exposures[2], "(1|hz)"), collapse = "+") ) )

    # Fit GLMM for outcome model
    fit_curr <- glmmTMB(form_curr, data = tsw_inc_p_timely, family = binomial(),
      weights = weights)
    model_parameters(fit_curr, exponentiate = TRUE)
    AIC(fit_curr)

    # Adjust for confounders (add one at a time, retain if it improves AIC / 
      # modifies effect size by > 10%)
      # proportion of successful SDBs
      form_new <- update(form_curr, ~ . + p_success_prev_3w)
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new
      
      # proportion of nosocomial cases
      form_new <- update(form_curr, ~ . + p_nosocomial_prev_3w )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

      # proportion of cases with a known epidemiological link
      form_new <- update(form_curr, ~ . + p_epi_link_prev_3w )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

      # Ebola vaccination coverage
      form_new <- update(form_curr, ~ . + vc_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

      # presence of an open Ebola treatment centre
      form_new <- update(form_curr, ~ . + etc_open_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude

      # attacks against EVD response
      form_new <- update(form_curr, ~ . + n_against_evd_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # rate of insecurity events
      form_new <- update(form_curr, ~ . + rate_events_3w_lag2 )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # number of phone networks
      form_new <- update(form_curr, ~ . + n_networks )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
            
      # number of radio frequencies
      form_new <- update(form_curr, ~ . + n_frequencies )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # presence of health facilities
      form_new <- update(form_curr, ~ . + rate_hf )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude
      
      # road network
      form_new <- update(form_curr, ~ . + rate_road_length )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude      
  
      # ratio of SDBs to population
      form_new <- update(form_curr, ~ . + r_sdb_to_pop_prev_3w )
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude

    # Explore plausible interactions 
      # exposure with proportion success
      form_new <- update(form_curr, ~ . + 
          p_success_prev_3w : p_timely_prev_3w_cat)
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # exclude: changes point estimate and lowers AIC, but coefficients
          # hard to interpret

      # exposure with ratio of SDB burials per population
      form_new <- update(form_curr, ~ . + 
          p_timely_prev_3w_cat : r_sdb_to_pop_prev_3w)
      fit_new <- update(fit_curr, formula = form_new)
      model_parameters(fit_new, exponentiate = TRUE)
      AIC(fit_new)
        # keep in model
        fit_curr <- fit_new
        form_curr <- form_new

    # Save model
    fit_inc_p_timely_pw <- fit_curr
    x <- model_parameters(fit_curr, exponentiate = TRUE) 
    out <- flextable(x[, c("Parameter", "Coefficient", "CI_low", "CI_high", "p",
      "Effects")])
    out <- set_formatter(out, 
      "Coefficient" = function (x) {sprintf("%.02f", x)},
      "CI_low" = function (x) {sprintf("%.02f", x)},
      "CI_high" = function (x) {sprintf("%.02f", x)},
      "p" = function (x) {sprintf("%.03f", x)}
    )  
    out <- autofit(out)
    flextable::save_as_docx(out, 
      path = paste(dir_output, "/out_glmm_inc_p_timely_pw.docx", sep =""))


#...............................................................................
### ENDS
#...............................................................................
  
