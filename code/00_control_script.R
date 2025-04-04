#...............................................................................
### +++ EFFECT OF SAFE AND DIGNIFIED BURIALS ON EBOV TRANSMISSION IN DRC +++ ###
#...............................................................................

#...............................................................................
## ------- 'CONTROL' R SCRIPT TO SOURCE INPUTS AND CALL OTHER SCRIPTS ------- ##
#...............................................................................

                              # Written by Francesco Checchi, LSHTM (July 2020)
                              # francesco.checchi@lshtm.ac.uk 


#...............................................................................
### Preparatory steps
#...............................................................................

  #...................................      
  ## Install or load required R packages

    # pacman and EpiEstim
    if (!"pacman" %in% rownames(installed.packages()))
      {install.packages("pacman")}
    if (!"EpiEstim" %in% rownames(installed.packages()))
      {install.packages('EpiEstim', 
        repos =c('https://mrc-ide.r-universe.dev', 
          'https://cloud.r-project.org'))
    }
    
    # Other packages
    pacman::p_load(boot, broom.mixed, cobalt, corrplot, EpiEstim, 
      flextable, gamlss, ggvenn, glmmTMB, ggplot2, ggpubr, influence.ME, 
      lattice, lme4, lubridate, MatchIt, pammtools, parameters, 
      readxl, scales, tidyverse, viridis, weights, zoo)
    
  #...................................      
  ## Starting setup

    # Clean up from previous code / runs
    rm(list=ls(all=T) )
  
    # Set font
    windowsFonts(Arial=windowsFont("Arial"))

    # Set working directories
    dir_code <- dirname(rstudioapi::getActiveDocumentContext()$path )
    setwd(dir_code)
    print( getwd() )
    dir_input <- gsub("/code", "/input", dir_code)
    dir_output <- gsub("/code", "/output", dir_code)
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
    palette_cb <- viridis(16)
    show_col(palette_cb)


  #...................................      
  ## Specify analysis parameters

    # Minimum number of confirmed + probable cases over the time series that a 
        # health zone must have to be eligible for analysis
    cases_threshold <- 5
  
    # Time interval (in weeks) for overall analysis
    t_overall <- 1
    
    # Number of time intervals over which to quantify SDB interventions, as well
      # as other epidemiological parameters
    window_transmission <- 2

    # Number of time intervals over which to evaluate effects of other
      # potentially confounding time-varying factors
    windows_other <- c(1, 2, 4)
    if (t_overall == 2) {windows_other <- 
      windows_other[! windows_other %in% seq(1, 9, 2)] / 2}

    # Lags (in time intervals) to explore for time-varying confounders
    lags <- c(0, 2)
    
    # Serial interval parameters (Jombart et al. Eurosurveillance)
    si_mean <- 15.3 # days
    si_sd <- 7.3 # days
    
    # Number of boostrap runs over which to estimate the net reproduction number 
      # using the Wallinga-Teunis method
    boot_rt <- 500
    
    # Whether to omit an initial period of 15 days (~ 1 serial interval) 
      # within each continuous stretch of reproduction number time series 
      # to avoid bias due to over-estimation of R at start of local outbreaks
    omit_initial <- T # default analysis used F
    
    # Assumed Ebola vaccine effectiveness (to estimate vaccination coverage 
      # based on case-coverage formula)
    ve <- 0.85

    # Whether the following variables should be smoothed across the time series:
    smooth_vcp <- F
    smooth_p_nosocomial <- F
    smooth_p_epi_link <- F
    
    
#...............................................................................
### Sourcing bespoke functions
#...............................................................................
    
source(paste0(dir_code, "/00_user_functions.R"))    
    
        
#...............................................................................
### Reading in required files
#...............................................................................

  #...................................      
  ## Variable dictionary
  filename <- paste0(dir_input, "/evd_drc_sdb_effect_datasets_pub.xlsx")
  dictionary <- as.data.frame(read_excel(filename, sheet= "dictionary"))
    
    # Dataset names
    dnames <- unique(dictionary[, "tab"])
      
  #...................................      
  ## Read in all the datasets
    # For each dataset...
    for (i in dnames) {
      # progress
      print(paste0("now reading in dataset of ", i))
      
      # read in
      x <- as.data.frame(read_excel(filename, sheet = i))
      assign(i,  x)
      
      # only keep needed columns
      cols_needed <- subset(dictionary, tab == i)[, "use"]
      cols_needed <- which(! cols_needed %in% c("no") )
      x <- get(i)[, cols_needed]
      assign(i, x)
    }
    

#...............................................................................                            
### Preparing data for analysis
#...............................................................................
    
source(paste0(dir_code, "/01_prepare_data.R"))    

    
#...............................................................................                            
### Exploring associations and computing propensity scores
#...............................................................................
    
source(paste0(dir_code, "/02_compute_propensity_scores.R"))    

    
#...............................................................................                            
### Estimating effect of SDB on change in incidence
#...............................................................................

source(paste0(dir_code, "/03_estimate_effect_incidence.R"))  

        
#...............................................................................                            
### Estimating effect of SDB on the reproduction number
#...............................................................................
    
source(paste0(dir_code, "/04_estimate_effect_reff.R"))         

    
#...............................................................................                            
### Estimating the dose-response effect of SDB success
#...............................................................................
    
source(paste0(dir_code, "/05_estimate_dose_response.R"))     

     
#...............................................................................
### ENDS
#...............................................................................


