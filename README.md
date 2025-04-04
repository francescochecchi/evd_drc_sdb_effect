## Effect of a safe and dignified burial intervention on Ebola virus transmission in the eastern Democratic Republic of the Congo, 2018-2019
## Repository of datasets and R scripts
March 2025

Main author: Francesco Checchi (London School of Hygiene and Tropical Medicine), francesco.checchi@lshtm.ac.uk

Partners: International Federation of Red Cross Red Crescent Societies (IFRC), Ministry of Health, Congolese Red Cross

Donor: Research for Health in Humanitarian Crises (R2HC) Programme (ref. 43670) of the Enhancing Learning and Research for Humanitarian Assistance (ELRHA)

### Project description
In this study, we analysed a variety of previously collected datasets to estimate the effect of the Safe and Dignified Burial (SDB) service implemented by the Congolese Red Cross, Civil protection units and other actors, and supported by the IFRC, during the Ebola virus epidemic in the eastern Democratic Republic of Congo (2018-2020; our analysis covers > 90% of all cases but ends in September 2019). The study outcomes were two measures of Ebola transmission, the change in incidence over successive time windows and the effective reproduction number. We used propensity score analysis to approach quasi-experimental conditions of allocation of the SDB exposure, thereby reducing confounding bias.

### Repository structure
The repository consists of three folders:
* folder `/input` contains a single Excel file that contains all data needed for the analysis. Each worksheet (tab) in this file contains a separate dataset. There is also a worksheet with a description of all the datasets, and a variable dictionary worksheet. The R script will automatically read each of the datasets.
* folder `/code` contains the R scripts needed to run the analysis (see below).
* folder `/output` is where the scripts will save output datasets, documents and figures. The folder has been pre-loaded with the outputs, but when the analysis is re-run these will automatically be replaced with the new outputs.

### Description of the R scripts
* `00_control_script.R` is the only script required to run the analysis without changes. This script installs and loads required R packages, sets needed parameters, sources functions, reads datasets and then sources all the other scripts in order. Altogether, the entire analysis should complete within 5-10 minutes on a standard laptop. You should interact with the other R scripts only if they need to be changed, or if a specific section of the analysis needs to be done.
* `00_user_functions.R` contains all of the bespoke functions written for this analysis, which will be then called by the other R scripts.
* `01_prepare_data.R` processes and merges together all the datasets, creating exposures (SDB timeliness and successfulness), transmission outcomes (change in incidence over adjacent time windows and effective reproduction number, the latter estimated within this script) and confounder variables. It also generates a few descriptive tables and graphs.
* `02_compute_propensity_scores.R` explores different lags, time windows and categorical versus continuous forms of key variables. For each transmission outcome and exposure, it then computes propensity scores, analyses and graphs resulting confounder balance, and readies dataset for PS-adjusted analysis.
* `03_estimate_effect_incidence.R` fits binomial and linear models to estimate the effect of the two exposures on the change in incidence outcome, using two methods of propensity score adjustment and also adding other confounders to the models. The final model coefficients are the main output.
* `04_estimate_effect_reff.R` does the same as the above code, but for the other transmission outcome, namely the effective reproduction number.
* `05_estimate_dose_response.R` computes and graphs the dose-response associations between SDB successfulness and each of the transmission outcomes.


### How to replicate the analysis
Analysts should make sure all the three folders are saved to the same directory, not altering the folder structure. The control script `sdb_effect_00_control.script.R`, when run, will automatically detect this drectory. An updated version of R software should be installed (https://www.r-project.org/). It is also recommended to run the code from the open-source RStudio interface (https://www.rstudio.com/products/rstudio/download/). Both R and RStudio are free and open-source.
