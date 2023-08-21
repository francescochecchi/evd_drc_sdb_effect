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

    # List of required packages
    pacs_needed <- c("boot", "broom.mixed", "cobalt", "corrplot", "EpiEstim", 
      "flextable", "gamlss", "glmmTMB", "ggplot2", "ggpubr", "influence.ME", 
      "lattice", "lme4", "lubridate", "parameters", "RColorBrewer", 
      "readxl", "scales", "tidyverse", "weights", "zoo")
    
    # Install any packages not yet installed
    already <- pacs_needed %in% row.names(installed.packages())
    if (any(already == FALSE)) { install.packages(pacs_needed[! already]) }

    # Load all packages    
    lapply(pacs_needed, library, character.only = TRUE)
    
  #...................................      
  ## Starting setup

    # Clean up from previous code / runs
    rm(list=ls(all=TRUE) )
  
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
    palette_cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
      "#0072B2", "#D55E00", "#CC79A7")
    show_col(palette_cb)


  #...................................      
  ## Specify analysis parameters

    # Minimum number of confirmed + probable cases over the time series that a 
        # health zone must have to be eligible for analysis
    cases_threshold <- 5
  
    # Time windows (in weeks) over which to quantify SDB interventions, as well
      # as other epidemiological parameters 
    windows_transmission <- c(1, 3)

    # Time windows (in weeks) over which to evaluate effects of other 
      # potentially confounding time-varying factors
    windows_other <- c(3, 5)
    
    # Lags (weeks) to explore for potentially confounding time-varying factors
    lags <- c(0, 2)
    
    # Serial interval parameters (Jombart et al. Eurosurveillance)
    si_mean <- 15.3 # days
    si_sd <- 7.3 # days
    
    # Number of boostrap runs over which to estimate the net reproduction number 
      # using the Wallinga-Teunis method
    boot_rt <- 500

    # Assumed Ebola vaccine effectiveness (to estimate vaccination coverage 
      # based on case-coverage formula)
    ve <- 0.85


#...............................................................................
### Sourcing bespoke functions
#...............................................................................
    
source(paste(dir_code, "/sdb_effect_00_user_functions.R", sep = ""), echo =TRUE)    
    
        
#...............................................................................
### Reading in required files
#...............................................................................

  #...................................      
  ## Variable dictionary
  filename <- paste(dir_input, "/evd_drc_sdb_effect_datasets_pub.xlsx", sep ="")
  dictionary <- as.data.frame(read_excel(filename, sheet= "dictionary"))
    
    # Dataset names
    dnames <- unique(dictionary[, "tab"])
      
  #...................................      
  ## Read in all the datasets
    # For each dataset...
    for (i in dnames) {
      # read in
      assign(i, as.data.frame(read_excel(filename, sheet = i)) )
      
      # only keep needed columns
      cols_needed <- subset(dictionary, tab == i)[, "use"]
      cols_needed <- which(! cols_needed %in% c("no") )
      x <- get(i)[, cols_needed]
      assign(i, x)
    }
    

#...............................................................................                            
### Preparing data for analysis
#...............................................................................
    
source(paste(dir_code, "/sdb_effect_01_prepare_data.R", sep = ""), echo = TRUE)    

    
#...............................................................................                            
### Exploring associations and computing propensity scores
#...............................................................................
    
source(paste(dir_code, "/sdb_effect_02_compute_propensity_scores.R", sep = ""), 
  echo = TRUE)    

    
#...............................................................................                            
### Estimating effect of SDB on change in incidence
#...............................................................................

source(paste(dir_code, "/sdb_effect_03_estimate_effect_incidence.R", sep = ""), 
  echo = TRUE)  

        
#...............................................................................                            
### Estimating effect of SDB on the reproduction number
#...............................................................................
    
source(paste(dir_code, "/sdb_effect_04_estimate_effect_reff.R", sep = ""), 
  echo = TRUE)         

    
#...............................................................................                            
### Estimating the dose-response effect of SDB success
#...............................................................................
    
source(paste(dir_code, "/sdb_effect_05_estimate_dose_response.R", sep = ""), 
  echo = TRUE)     

     
#...............................................................................
### ENDS
#...............................................................................


