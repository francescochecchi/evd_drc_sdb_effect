#...............................................................................
### +++ EFFECT OF SAFE AND DIGNIFIED BURIALS ON EBOV TRANSMISSION IN DRC +++ ###
#...............................................................................

#...............................................................................
## ------------------ UNUSED FUNCTIONS - NOT FULLY TESTED ------------------- ##
#...............................................................................

                              # Written by Francesco Checchi, LSHTM (July 2020)
                              # francesco.checchi@lshtm.ac.uk


#...............................................................................
### Function to examine model diagnostics of OLS fit
#...............................................................................

f_diag_lm <- function(fit_f = fit) {

  # Fitted versus residuals (normal distribution of residuals)
  print(hist(residuals(fit_f), xlab="residuals", ylab="percent of total",
    col="grey", cex.main=1, main = paste("distribution of residuals for ",
      "\n", fit_lm$terms[[2]] , sep="") ) )

  # Homoskedasticity
  print (plot(fitted(fit_f), residuals(fit_f), xlab="fitted values",
    ylab="residuals", col="red", cex.main=1,
    main = paste("residuals vs. fitted values for ", "\n",
      fit_lm$terms[[2]] , sep="")) )
  print( abline(a=0, b=0, h=c(-4,4), lty=3, col="grey") )

  # Quantile-quantile plot
  print (qqnorm(residuals(fit_f)/sd(residuals(fit_f)),xlim=c(-4,4),ylim=c(-4,4),
    main = paste("quantile-quantile plot of residuals for ", "\n",
      fit_lm$terms[[2]] , sep=""),
    ylab="standardised residual quantiles", xlab="theoretical quantiles",
    col="red", cex.main=1) )
  print( qqline(residuals(fit_f)/sd(residuals(fit_f)) , lty=3, col="grey" ) )
}



#.........................................................................................
### Function to fit a fixed-effects linear model for longitudinal data and display clean results
#.........................................................................................

f_lm <- function(vars_f, data_f, wt_f, window_transm_f, f_log) {

  # identify dependent variable and weights
  if (f_log == TRUE) {x1 <- c( paste("ln_r_mean_curr_", window_transm_f, "w", sep=""),
                               paste("r_wt_curr_", window_transm_f, "w", sep="") )}
  if (f_log == FALSE) {x1 <- c( paste("r_mean_curr_", window_transm_f, "w", sep=""),
                               paste("r_wt_curr_", window_transm_f, "w", sep="") )}
  # write the model formula
  form <- as.formula( paste(x1[1], "~", paste(vars_f, collapse= " + "), sep="")  )

  # fit linear model, with or without weights
  if (wt_f == TRUE) { fit <- lm(form, data = data_f, weights = data_f[, x1[2]]) }
  if (wt_f == FALSE) { fit <- lm(form, data = data_f) }
  return(fit)
}



#.........................................................................................
### Function to fit a mixed linear model for longitudinal data and display clean results
#.........................................................................................

f_lmm <- function(vars_f, data_f, wt_f, window_transm_f, f_log) {
  # identify dependent variable and weights
  if (f_log == TRUE) {x1 <- c( paste("ln_r_mean_curr_", window_transm_f, "w", sep=""),
                               paste("r_wt_curr_", window_transm_f, "w", sep="") )}
  if (f_log == FALSE) {x1 <- c( paste("r_mean_curr_", window_transm_f, "w", sep=""),
                               paste("r_wt_curr_", window_transm_f, "w", sep="") )}
  # write the model formula
  form <- as.formula( paste(x1[1], "~", paste(vars_f, collapse= " + "), "+ (1|hz)", sep="")  )

  # fit linear model with random effect, with or without weights
  if (wt_f == TRUE) { fit <- lmer(form, data = data_f, REML = FALSE, weights = data_f[, x1[2]]) }
  if (wt_f == FALSE) { fit <- lmer(form, data = data_f, REML = FALSE) }
  return(fit)
}









...............................................................................
### Functions to compute multivariate generalised propensity scores (GPS)
#.............................................................................
  # Credit: Justin Williams (see https://arxiv.org/pdf/2008.13767.pdf)
  # https://github.com/williazo/mvGPS/blob/master/R/mvGPS.R
  # dependent package WeightIt did not load successfully so adapted functions from Github repo


  #...................................
  ## Main function

mvGPS <- function(D, C, common=FALSE, trim_w=FALSE, trim_quantile=0.99){
    check_result <- D_C_check(D, C, common)
    assign("D", check_result$D)
    assign("C", check_result$C)

    m <- ncol(D)

    for(i in seq_len(m)){
        if(i==1){
            #marginal densities factorized
            d_1 <- lm(D[, i] ~ 1)
            d_1_mu <- coef(d_1)
            d_1_sigma <- summary(d_1)$sigma
            f_d_1 <- dnorm(D[, i], mean=d_1_mu, sd=d_1_sigma)

            #generalized propensity score
            gps_d_1 <- lm(D[, i] ~ C[[i]] - 1)
            gps_1_beta <- coef(gps_d_1)
            gps_1_Xb <- model.matrix(gps_d_1) %*% gps_1_beta
            gps_1_sigma <- summary(gps_d_1)$sigma
            f_gps_1 <- dnorm(D[, i], mean=gps_1_Xb, sd=gps_1_sigma)
        } else {
            cond_dens <- lapply(seq_len(m - 1) + 1, function(x){
                #full conditional marginal densities
                d_x <- lm(D[, x] ~ D[, seq_len(x-1)])
                d_x_beta <- coef(d_x)
                d_x_Xb <- model.matrix(d_x) %*% d_x_beta
                d_x_sigma <- summary(d_x)$sigma
                f_d_x <- dnorm(D[, x], mean=d_x_Xb, sd=d_x_sigma)

                #full conditional generalized propensity scores
                gps_x <- lm(D[, x] ~ D[, seq_len(x-1)] + C[[x]] - 1)
                gps_x_beta <- coef(gps_x)
                gps_x_Xb <- model.matrix(gps_x) %*% gps_x_beta
                gps_x_sigma <- summary(gps_x)$sigma
                f_gps_x <- dnorm(D[, x], gps_x_Xb, gps_x_sigma)

                return(list(marg=f_d_x, gps=f_gps_x))
            })
        }
    }
    cond_results <- unlist(cond_dens, recursive=FALSE)
    num_args <- cond_results[which(names(cond_results)=="marg")]
    num_args[["marg_1"]] <- f_d_1
    denom_args <- cond_results[which(names(cond_results)=="gps")]
    denom_args[["gps_1"]] <- f_gps_1

    score <- Reduce("*", denom_args)
    w <- Reduce("*", num_args)/score
    if(trim_w==TRUE){
        #trimming the large weights
        w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
    }
    return(list(score=score, wts=w))
}


  #...................................
  ## Internal function to check that inputs to mvGPS function are suitable

D_C_check <- function(D, C, common){
    D <- as.matrix(D)
    m <- ncol(D)
    if(is.null(colnames(D))) colnames(D) <- paste0("D", seq_len(m))
    n <- nrow(D)
    if(common){
        C <- as.matrix(C)
        if(is.null(colnames(C))) colnames(C) <- paste0("C", seq_len(ncol(C)))
        if(is.list(C)) stop("common=TRUE, expecting C to be single matrix of common confounders", call.=FALSE)
        C <- rep(list(C), m)
    } else {
        if(!is.list(C)) stop("common=FALSE, C must be list of length m", call.=FALSE)
        C <- lapply(seq_len(length(C)), function(x){
            X <- as.matrix(C[[x]])

            if(is.null(colnames(X))){
                p_j <- ncol(X)
                colnames(X) <- paste0("C_",x, "_", seq_len(p_j))
            }
            return(X)
            })

    }
    C_k <- unlist(lapply(C, ncol))
    C_n <- unlist(lapply(C, nrow))
    if(m < 2) stop("Exposure must be multivariate. See details to ensure formula is properly specified", call.=FALSE)
    if(!all(C_n==n)) stop("Each matrix in C must have same number of observations, n, as D", call.=FALSE)
    if(length(C)!=m) stop("Set of confounders not equal to number of exposures, m.")
    return(list(D=D, C=C))
}


................................................................................
### Functions to output components to estimate bivariate mvGPS at any new exposure value, given data
#.............................................................................
  # Adapted from Justin Williams (see https://arxiv.org/pdf/2008.13767.pdf)
  # https://github.com/williazo/mvGPS/blob/master/R/mvGPS.R


  #...................................
  ## Main function

f_mvGPS_fits <- function(D, C, common=FALSE, trim_w=FALSE, trim_quantile=0.99){
  check_result <- f_mod_mvGPS_chk(D, C, common)
  assign("D", check_result$D)
  assign("C", check_result$C)

  # Exposure 1
    #marginal densities factorized
    d_1 <- lm(D[, 1] ~ 1)
    d_1_mu <- coef(d_1)
    d_1_sigma <- summary(d_1)$sigma

    #generalized propensity score
    gps_d_1 <- lm(D[, 1] ~ C[[1]] - 1)
    gps_1_beta <- coef(gps_d_1)
    gps_1_Xb <- model.matrix(gps_d_1) %*% gps_1_beta
    gps_1_sigma <- summary(gps_d_1)$sigma

  # Exposure 2
    #full conditional marginal densities
    d_2 <- lm(D[, 2] ~ D[, 1])
    d_2_beta <- coef(d_2)
    d_2_Xb <- model.matrix(d_2) %*% d_2_beta
    d_2_sigma <- summary(d_2)$sigma

    #full conditional generalized propensity scores
    gps_d_2 <- lm(D[, 2] ~ D[, 1] + C[[2]] - 1)
    gps_2_beta <- coef(gps_d_2)
    gps_2_Xb <- model.matrix(gps_d_2) %*% gps_2_beta
    gps_2_sigma <- summary(gps_d_2)$sigma

return(list(
  "d_1_mu" = d_1_mu, "d_1_sigma" = d_1_sigma, "gps_1_Xb" = gps_1_Xb, "gps_1_sigma" = gps_1_sigma,
  "d_2_Xb" = d_2_Xb, "d_2_sigma" = d_2_sigma, "gps_2_Xb" = gps_2_Xb, "gps_2_sigma" = gps_2_sigma
  ))

}

  #...................................
  ## Internal function to check that inputs to mvGPS function are suitable

f_mvGPS_chk <- function(D, C, common){
    D <- as.matrix(D)
    m <- ncol(D)
    if(is.null(colnames(D))) colnames(D) <- paste0("D", seq_len(m))
    n <- nrow(D)
    if(common){
        C <- as.matrix(C)
        if(is.null(colnames(C))) colnames(C) <- paste0("C", seq_len(ncol(C)))
        if(is.list(C)) stop("common=TRUE, expecting C to be single matrix of common confounders", call.=FALSE)
        C <- rep(list(C), m)
    } else {
        if(!is.list(C)) stop("common=FALSE, C must be list of length m", call.=FALSE)
        C <- lapply(seq_len(length(C)), function(x){
            X <- as.matrix(C[[x]])

            if(is.null(colnames(X))){
                p_j <- ncol(X)
                colnames(X) <- paste0("C_",x, "_", seq_len(p_j))
            }
            return(X)
            })

    }
    C_k <- unlist(lapply(C, ncol))
    C_n <- unlist(lapply(C, nrow))
    if(m < 2) stop("Exposure must be multivariate. See details to ensure formula is properly specified", call.=FALSE)
    if(!all(C_n==n)) stop("Each matrix in C must have same number of observations, n, as D", call.=FALSE)
    if(length(C)!=m) stop("Set of confounders not equal to number of exposures, m.")
    return(list(D=D, C=C))
}



#.............................................................................
### Function to compute multivariate generalised propensity scores (GPS) for any exposure combination
#.............................................................................
  # Adapted from Justin Williams (see https://arxiv.org/pdf/2008.13767.pdf)
  # https://github.com/williazo/mvGPS/blob/master/R/mvGPS.R


f_mvGPS_comp <- function(exp_value1, exp_value2, fits_mvgps){
  # exp_value_1  and _2 are the values of the two exposures
  # fits_mvgps is the output of function f_mvGPS_fits

  # Exposure 1
    # marginal densities factorized
    f_d_1 <- dnorm(exp_value_1, mean = fits_mvgps[["d_1_mu"]], sd = fits_mvgps[["d_1_sigma"]])

    # generalized propensity score
    f_gps_1 <- dnorm(exp_value_1, mean = fits_mvgps[["gps_1_Xb"]], sd = fits_mvgps[["gps_1_sigma"]])

  # Exposure 2
    # full conditional marginal densities
    f_d_2 <- dnorm(exp_value_2, mean = fits_mvgps[["d_2_Xb"]], sd = fits_mvgps[["d_2_sigma"]])

    # full conditional generalized propensity scores
    f_gps_2 <- dnorm(exp_value_2, fits_mvgps[["gps_2_Xb"]], fits_mvgps[["gps_2_sigma"]])

    # GPS
    gps <- Reduce('*', list(f_gps_1, f_gps_2) )

    # GPS weight
    gps_wt <- Reduce('*', list(f_d_1, f_d_2) ) / gps

  return(data.frame(gps, gps_wt))
}



#...............................................................................
### ENDS
#...............................................................................


