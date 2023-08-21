#...............................................................................
### +++ EFFECT OF SAFE AND DIGNIFIED BURIALS ON EBOV TRANSMISSION IN DRC +++ ###
#...............................................................................

#...............................................................................
## -------------------------- FUNCTIONS FOR ANALYSIS ------------------------ ##
#...............................................................................

                              # Written by Francesco Checchi, LSHTM (July 2020)
                              # francesco.checchi@lshtm.ac.uk



#...............................................................................
### Function to examine model diagnostics of LMM fit
#...............................................................................

f_diag_lmm <- function(fit_f = fit) {
  
  # Fitted versus residuals (normal distribution of residuals)
  print(hist(resid(fit_f), xlab="residuals", ylab="percent of total",
    col=palette_cb[1], cex.main=1,
    main = paste("distribution of residuals for ", "\n", formula(fit_f)[2] ,
      sep="") ) )

  # Homoskedasticity
  print (plot(fitted(fit_f), resid(fit_f),xlab="fitted values",ylab="residuals",
    col=palette_cb[7], cex.main=1,
    main = paste("residuals vs. fitted values for ", "\n", formula(fit_f)[2] ,
    sep="")) )
  print( abline(a=0, b=0, h=c(-4,4), lty=3, col=palette_cb[1]) )

  # Quantile-quantile plot
  print (qqnorm(resid(fit_f)/sd(resid(fit_f)), xlim=c(-4,4), ylim=c(-4,4) ,
    main = paste("quantile-quantile plot of residuals for ", "\n",
    formula(fit_f)[2] , sep=""),
    ylab="standardised residual quantiles", xlab="theoretical quantiles",
    col=palette_cb[7], cex.main=1) )
  print( qqline(resid(fit_f)/sd(resid(fit_f)) , lty=3, col=palette_cb[1] ) )
}


#...............................................................................
### Function to fit a fixed effects general linear model for group logit data 
  # and display clean results
#...............................................................................

f_glm <- function(vars_f, data_f, wt_f, window_transm_f) {
  
  # Identify dependent variables and weights
  x <- c(paste("cases_curr_", window_transm_f, "w", sep=""), 
    paste("cases_prev_", window_transm_f, "w", sep=""),
    paste("cases_wt_", window_transm_f, "w", sep=""))

  # Write the model formula
  form <- as.formula( paste("cbind(", x[1], ", ", x[2], ")", "~", 
    paste(vars_f, collapse= " + "), sep="")  )

  # Fit model with fixed effects only, with or without weights
  if (wt_f == TRUE) { fit <- glm(form, data = data_f, family = binomial(), 
    weights =  data_f[, x[3]]) }
  if (wt_f == FALSE) { fit <- glm(form, data = data_f, family = binomial() ) }
  return(fit)
}



#...............................................................................
### Function to fit a mixed general linear model for longitudinal group logit 
  # data and display clean results
#...............................................................................

f_glmm <- function(vars_f, data_f, wt_f, window_transm_f) {
  
  # Identify dependent variables and weights
  x <- c(paste("cases_curr_", window_transm_f, "w", sep=""), 
    paste("cases_prev_", window_transm_f, "w", sep=""),
    paste("cases_wt_", window_transm_f, "w", sep=""))
  
  # Write the model formula
  form <- as.formula( paste("cbind(", x[1], ", ", x[2], ")", "~", 
    paste(vars_f, collapse= " + "), "+ (1|hz)", sep="")  )
  
  # Fit model with random effect, with or without weights
  if (wt_f == TRUE) { fit <- glmer(form, data = data_f, family = binomial(), 
    weights = data_f[, x[3]]) }
  if (wt_f == FALSE) { fit <- glmer(form, data = data_f, family = binomial()) }
  return(fit)
}



#...............................................................................
### Function to compute Generalised Propensity Scores and confounder balance
  # based on Hirano & Imbens (2004), 
    # https://www.math.mcgill.ca/dstephens/PSMMA/Articles/HIrano-Imbens-2004.pdf
  # and Austin (2017), 
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5969262/#sim7615-bib-0003 
#...............................................................................

f_gps <- function(data_f, exposure_f, confounders_f, colour_f) {

  #...................................
  ## Compute Generalised Propensity Scores
    # Define formula for exposure as a function of confounders
    formula_exp <- formula(paste(exposure_f," ~ ", 
      paste(confounders_f, collapse = " + "), "+ (1 | hz)"))

    # Linear mixed model of exposure as a function of confounders
      # (HZ = random effect)
    fit_exp <- lmer(formula = formula_exp, data = data_f)

    # Compute GPS values
    data_f$gps <- dnorm(x = data_f[, exposure_f], 
      mean = fitted.values(fit_exp), sd = summary(fit_exp)$sigma)
      
    # Distribution of GPS
    f_hist("gps", data_f, c(NA,NA))

  #...................................
  ## Compute GPS stabilised weights

    # Numerators
    data_f$wt_num <- dnorm(x = (data_f[, exposure_f] - 
      mean(data_f[, exposure_f]) ) / sd(data_f[, exposure_f]), mean = 0, sd = 1)

    # Weights
    data_f$gps_wt <- data_f$wt_num / data_f$gps
      
    # Truncate to avoid extreme weights
    data_f$gps_wt <- ifelse(data_f$gps_wt > quantile(data_f$gps_wt, probs = 0.99, 
      na.rm = TRUE),
      quantile(data_f$gps_wt, probs = 0.99, na.rm = TRUE), data_f$gps_wt )
    data_f$gps_wt <- ifelse(data_f$gps_wt < quantile(data_f$gps_wt, probs = 0.01, 
      na.rm = TRUE),
      quantile(data_f$gps_wt, probs = 0.01, na.rm = TRUE), data_f$gps_wt )
    

  #...................................
  ## Check balance of confounders with adjustment by GPS
    # with adjustment by GPS, correlations of confounders with outcome 
    # should be lower (ideally all < 0.10)

    # Divide data into quantiles by exposure, compute medians of each quantile
    exp_q <- quantile(data_f[, exposure_f], probs = c(0.33333, 0.66667, 1))
      # 3 groups arbitrarily, to avoid data sparsity
    data_f$exp_q <- cut(data_f[, exposure_f], breaks = c(0, exp_q), 
      labels = 1:(length(exp_q)),
      ordered_result = TRUE, include.lowest = TRUE)
    exp_mids <- aggregate(data_f[, exposure_f], by = list(exp_q = data_f$exp_q), 
      FUN = median)
    colnames(exp_mids) <- c("exp_q", "median")

    # Prepare output for balance-with-adjustment statistics
    out <- expand.grid(as.character(levels(data_f$exp_q)), as.character(1:4), 
      confounders_f)
    colnames(out) <- c("exp_q", "gps_q", "confounder")
    out[, c("mean_diff", "se_diff", "n_obs", "n_obs_k", "n_obs_nonk")] <- NA
    out$exp_q <- as.ordered(out$exp_q)
    out$gps_q <- as.ordered(out$gps_q)

    # For each exposure_f quantile...
    for (k in levels(data_f$exp_q)) {
      # compute GPS using exposure model by fixing exposure to quantile median
      data_f$gps_mid <- dnorm(x = rep(exp_mids[exp_mids$exp_q == k, "median"], 
        nrow(data_f)),mean = fitted.values(fit_exp),sd = summary(fit_exp)$sigma)

      # divide data into quantiles by GPS as estimated above
      gps_q <- quantile(data_f$gps_mid, probs = seq(0, 1, by = 0.25))
      data_f$gps_q <- cut(data_f$gps_mid, breaks = gps_q, 
        labels = 1:(length(gps_q) - 1), ordered_result = TRUE, 
        include.lowest = TRUE)

      # now, for each GPS quantile...
      for (l in levels(data_f$gps_q)) {
        # select all data in this level
        x1 <- subset(data_f, gps_q == l)

        # split data by whether the exposure quantile == k or not
        x2 <- subset(x1, exp_q == k)
        x3 <- subset(x1, exp_q != k)

        # for each confounder...
        for (j in confounders_f) {
          # find position in output vector
          x4 <- which(out$exp_q == k & out$gps_q == l & out$confounder == j)

          # compute difference of means and its standard error
          out[x4, "mean_diff"] <- mean(x2[, j]) - mean(x3[, j])
          out[x4, "se_diff"] <- sqrt( (sd(x2[, j])^2 / nrow(x2)) +  
            (sd(x3[, j])^2 / nrow(x3))  )
        }

        # capture number of observations within l level
        out[which(out$exp_q == k & out$gps_q == l), "n_obs"] <- nrow(x1)
        out[which(out$exp_q == k & out$gps_q == l), "n_obs_k"] <- nrow(x2)
        out[which(out$exp_q == k & out$gps_q == l), "n_obs_nonk"] <- nrow(x3)
      }
    }

    # Compute t-statistics for association of confounder with exposure, 
      # after GPS adjustment

      # calculate weighted mean difference and SE of difference for each 
        # exposure quantile and confounder
      out$mean_diff_wts <- out$mean_diff * out$n_obs
      out$se_diff_wts <- out$se_diff * out$n_obs
      out <- aggregate(out[, c("n_obs", "n_obs_k", "n_obs_nonk", 
        "mean_diff_wts", "se_diff_wts")],
        by = out[, c("exp_q", "confounder")], FUN = sum)
      out$wt_mean_diff <- out$mean_diff_wts / out$n_obs
      out$wt_se_diff <- out$se_diff_wts / out$n_obs

      # compute t-statistic for each exposure quantile compared to other 
        # quantiles, and for each confounder
      out$t_value <- out$wt_mean_diff / out$wt_se_diff
      out$t_value <- round(out$t_value, 2)

      # reshape wide and make other changes
      out_adj <- reshape2::dcast(out[, c("exp_q", "confounder", "t_value")], 
        formula = confounder ~ exp_q)
      colnames(out_adj)[colnames(out_adj) != "confounder"] <- 
        paste("adj", levels(data_f$exp_q), sep = "_")
      out_adj$confounder <- as.character(out_adj$confounder)
      out_adj <- out_adj[order(out_adj$confounder), ]

  #...................................
  ## Check balance of confounders without adjustment by GPS

    # Prepare output for balance-without-adjustment statistics
    out <- expand.grid(as.character(levels(data_f$exp_q)), confounders_f)
    colnames(out) <- c("exp_q", "confounder")
    out[, c("mean_diff", "se_diff", "n_obs", "n_obs_k", "n_obs_nonk")] <- NA
    out$exp_q <- as.ordered(out$exp_q)

    # For each exposure quantile...
    for (k in levels(data_f$exp_q)) {

      # split data by whether the exposure quantile == k or not
      x2 <- subset(data_f, exp_q == k)
      x3 <- subset(data_f, exp_q != k)

      # for each confounder...
      for (j in confounders_f) {
        # find position in output vector
        x4 <- which(out$exp_q == k & out$confounder == j)

        # compute difference of means and its standard error
        out[x4, "mean_diff"] <- mean(x2[, j]) - mean(x3[, j])
        out[x4, "se_diff"] <- sqrt( (sd(x2[, j])^2 / nrow(x2)) +  
            (sd(x3[, j])^2 / nrow(x3))  )
      }

      # capture number of observations within k quantile
      out[which(out$exp_q == k), "n_obs"] <- nrow(data_f)
        out[which(out$exp_q == k ), "n_obs_k"] <- nrow(x2)
        out[which(out$exp_q == k ), "n_obs_nonk"] <- nrow(x3)
    }

    # Compute t-statistics for association of confounder with exposure, 
      # without GPS adjustment

      # compute t-statistic for each exposure quantile compared to other 
        # quantiles, and for each confounder
      out$t_value <- out$mean_diff / out$se_diff
      out$t_value <- round(out$t_value, 2)

      # reshape wide and make other changes
      out_unadj <- reshape2::dcast(out[, c("exp_q", "confounder", "t_value")], 
        formula = confounder ~ exp_q)
      colnames(out_unadj)[colnames(out_unadj) != "confounder"] <- 
        paste("unadj", levels(data_f$exp_q), sep = "_")
      out_unadj$confounder <- as.character(out_unadj$confounder)
      out_unadj <- out_unadj[order(out_unadj$confounder), ]


  #...................................
  ## Compare balance without and with adjustment by GPS

    # Merge outputs
    out <- merge(out_unadj, out_adj, by = "confounder")

    # Plot change in t-values for each confounder, before and after adjustment
      # prepare long version of output
      out_unadj$method <- "unadjusted"
      colnames(out_unadj) <- gsub("unadj_", "", colnames(out_unadj))
      out_adj$method <- "GPS-adjusted"
      colnames(out_adj) <- gsub("adj_", "", colnames(out_adj))
      df <- rbind(out_unadj, out_adj)
      df <- reshape2::melt(data = df, id.vars = c("confounder", "method"), 
        variable.name = "quantile", value.name = "t_statistic")
      df$method <- factor(df$method, levels = c("unadjusted", "GPS-adjusted"))
      df$quantile <- as.numeric(df$quantile)
      df[which(df$quantile == 1), "quantile"] <- "1st quantile"
      df[which(df$quantile == 2), "quantile"] <- "2nd quantile"
      df[which(df$quantile == 3), "quantile"] <- "3rd quantile"

      # plot
      plot <- ggplot(df, aes(y = t_statistic, x = method, group = confounder)) +
        scale_y_continuous(name = "t-statistic", 
          breaks = c(-5, -1.96, 0, 1.96, 5)) +
        scale_x_discrete("", expand = c(0.1, 0.1)) +
        theme_bw() +
        geom_line(colour = colour_f, alpha = 0.7) +
        annotate(geom = "rect", fill = colour_f, alpha = 0.3,
          xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96) +
        facet_wrap(~quantile, ncol = 2, nrow = 2)

  #...................................
  ## Return outputs
  return(list(dataset = data_f, balance_stats = out, balance_plot = plot) )

}


#...............................................................................
### Function to fit and output Generalised Propensity Scores model only
#...............................................................................

f_gps_lite <- function(data_f, exposure_f, confounders_f) {

  # Define formula for exposure as a function of confounders
  formula_exp <- formula(paste(exposure_f," ~ ", 
    paste(confounders_f, collapse = " + "), "+ (1 | hz)"))

  # Linear mixed model of exposure_f as a function of confounders_f 
    # (HZ = random effect)
  fit_exp <- lmer(formula = formula_exp, data = data_f)

  # Output model fit
  return(fit_exp)
}


#...............................................................................
### Function to plot histograms of any variable
#...............................................................................

f_hist <- function(var_f, data_f, lims_f) {

  # start plot
  plot <- ggplot(data_f)

    # if the variable has >= 20 unique values...
    if (length(unique(na.omit(data_f[, var_f]))) >= 20) {
      plot <- plot + 
        geom_histogram(aes(x = as.numeric(data_f[, var_f]) ), 
          color = palette_cb[4], fill = palette_cb[4],
          alpha = 0.5) + 
        theme_bw() + 
        xlab(var_f) + 
        scale_x_continuous(expand = c(0, 0), limits = lims_f )
    }

    # otherwise...
    if (length(unique(na.omit(data_f[, var_f]))) < 20) {
      plot <- plot + 
        geom_histogram(aes(x = as.numeric(data_f[, var_f]) ), stat = "count", 
          color = palette_cb[4], fill = palette_cb[4], alpha = 0.5) + 
        theme_bw() + 
        xlab(var_f) + 
        scale_x_continuous(expand = c(0, 0), limits = lims_f )
    }

  print(plot)
}


#...............................................................................
### Function to compute prediction confidence intervals or prediction profiles 
  # for bootstrapping, from a GAMLSS fit, by posterior simulation: see
  # https://r.789695.n4.nabble.com/Prediction-interval-with-GAM-td3460175.html 
#...............................................................................

f_interval <- function(fit_f = fit, n_boot_f = 1000, data_f, profile_f = FALSE){
  
  # Extract coefficient estimates and variance-covariance matrix from fit
  beta <- na.omit( coef(fit_f) ) # only mu coefficients, i.e. 
    # for estimating linear predictor; omit NAs as random effects come out as NA
    # convert to inverse (for later operation)
    beta_inv <- ginv(matrix(beta))
  vcov_matrix <- vcov(fit_f)[1:length(beta), 1:length(beta)] # only mu coeff's
    # Cholesky decomposition
    cholesky_dec <- chol(vcov_matrix)
  
  # Simulate a large number of random beta coefficient values 
    # drawn from an assumed normal error distribution
    # around their posterior estimates
  beta_sim <- t(cholesky_dec) %*% matrix(rnorm(n_boot_f * length(beta)), 
    length(beta), n_boot_f) + as.vector(beta) 
  
  # Predict both mu and sigma on the new data
  mu_pred <- predict(fit_f, newdata = data_f[, all.vars(formula(fit_f))])
  sigma_pred <- predict(fit_f, newdata = data_f[, all.vars(formula(fit_f))], 
    what = "sigma") # for later
  
  # Produce a 'linear prediction matrix'
    # since lp_matrix %*% coef(fit_f) = mu_pred , the code below uses the 
    # inverse (matrix 'division') to get lp_matrix, which gamlss doesn't output
  lp_matrix <- mu_pred %*% beta_inv 
  
  # Compute point estimate linear predictions and back-transform them 
    # if appropriate, for all the random beta values
    # only log and logit back-transform implemented here for now
  pred_sim <- lp_matrix %*% beta_sim
  if ( fit_f$mu.link == "logit" ) {pred_sim <- inv.logit(pred_sim)}
  if ( fit_f$mu.link == "log" ) {pred_sim <- exp(pred_sim)}
  
  # Lastly, generate random predictions by combining the point estimates 
    # with the other estimated distributional parameters
    # only a few distributions implemented here so far
  if ("PO" %in% family(fit_f)) 
    { rand_sim <- matrix( rPO(n = prod(dim(pred_sim)), mu = pred_sim), 
      nrow(pred_sim), ncol(pred_sim) ) }
  if ("NBI" %in% family(fit_f)) 
    { rand_sim <- matrix( rNBI(n = prod(dim(pred_sim)), mu = pred_sim, 
      sigma = rep(exp(sigma_pred), n_boot_f) ), 
      nrow(pred_sim), ncol(pred_sim) ) }
  if ("NBII" %in% family(fit_f)) 
    { rand_sim <- matrix( rNBII(n = prod(dim(pred_sim)), mu = pred_sim, 
      sigma = rep(exp(sigma_pred), n_boot_f) ), 
      nrow(pred_sim), ncol(pred_sim) ) }
  if ("BE" %in% family(fit_f)) 
    { rand_sim <- matrix( rBE(n = prod(dim(pred_sim)), mu = pred_sim, 
      sigma = rep(inv.logit(sigma_pred), n_boot_f) ), 
      nrow(pred_sim), ncol(pred_sim) ) }
  
  # If desire 80% and 95% confidence interval...
  if (profile_f == FALSE) {
    point_est <- mu_pred
    if ( fit_f$mu.link == "log") { point_est <- exp(mu_pred) }
    if ( fit_f$mu.link == "logit") { point_est <- inv.logit(mu_pred) }
    
    out <- cbind(data_f,
      point_est,
      t(apply(rand_sim, 1, quantile, prob = c(0.025,0.975), na.rm = TRUE)),
      t(apply(rand_sim, 1, quantile, prob = c(0.10,0.90), na.rm = TRUE))
    )
    out <- as.data.frame(out)
    colnames(out) <- c(colnames(data_f), "pred", "lci025", "uci975", "lci100",
      "uci900")
    return(out)
  }
  
  # ...else output entire prediction profile, sorted ascendingly
  if (profile_f == TRUE) {
    out <- apply(rand_sim, 1, sort)
    return(out)           
  }
  
}



#...............................................................................
### Function to compute indicators over variously lagged time windows of 
    # specified size (weeks), for indicators not directly affecting transmission
#...............................................................................

f_window_other <- function(data_f = tsw, hz_f = hz_ok, var_f, window_size_f, 
  lag_f, op_f = "sum") {

  # Prepare output
  out <- c()

  # For each health zone...
  for (i in hz_f) {
    
    # extract time series of variable of interest, and compute its running sum 
      # or mean over the given window size
    if (op_f == "sum") {out_i <- rollsum(subset(data_f, hz == i)[, var_f], 
      window_size_f, fill = NA, align = "right") }
    if (op_f == "mean") {out_i <- rollmean(subset(data_f, hz == i)[, var_f], 
      window_size_f, fill = NA, align = "right") }
    
    # lag running indicator
    out_i <- c(rep(NA, lag_f), out_i[1:(length(out_i) - lag_f)])
    
    # add output
    out <-c(out, out_i)
  }

  # Clean up NaN and Inf values
  out <- ifelse(out %in% c(NA, NaN, Inf, -Inf), NA, out)

  # Return output
  return(out)
}



#...............................................................................
### Function to compute indicators over previous and current transmission time 
  # window of specified size (in weeks); windows spaced 2 weeks apart or 
  # as many weeks as are needed to ensure windows do not overlap
#...............................................................................

f_window_transm <- function(data_f = cases_w, hz_f = hz, vars_f = "cases", 
  window_size_f = 3, percent_f = FALSE, op_f = "sum") {

  # Figure out spacing between windows
  spacing <- max(2, window_size_f)

  # Prepare output
  out_prev <- c() # previous time window
  out_curr <- c() # current time window

  # For each health zone...
  for (i in hz_f) {
    
    # extract time series of variables of interest into sequentially numbered 
      # vectors, and compute its running sum over the given window size
    for (j in 1:length(vars_f) ) {
      
      if (op_f == "sum") {
      assign(paste("x", j,sep=""), rollsum(subset(data_f, hz == i)[, vars_f[j]], 
        window_size_f, fill = NA, align = "center") )
      }
      
      if (op_f == "mean") {
      assign(paste("x", j,sep=""),rollmean(subset(data_f, hz == i)[, vars_f[j]], 
        window_size_f, fill = NA, align = "center") )
      }
    }
    
    # compute indicator
    if (length(vars_f) == 1 ) { curr_i <- x1 }
    if (length(vars_f) == 2 ) { curr_i <- x1 / x2 ; 
      if (percent_f == TRUE) {curr_i <- curr_i * 100} }

    # indicator in the previous time window
    prev_i <- c(rep(NA, spacing), 
      curr_i[1 : (length(curr_i) - spacing)])
    out_prev <- c(out_prev, prev_i)

    # indicator in the current time window
    out_curr <-c(out_curr, curr_i)
  }

  # Clean up NaN and Inf values
  out_prev <- ifelse(out_prev %in% c(NA, NaN, Inf, -Inf), NA, out_prev)
  out_curr <- ifelse(out_curr %in% c(NA, NaN, Inf, -Inf), NA, out_curr)

  # Return output
  return( cbind(out_prev, out_curr) )
 }



#...............................................................................
### ENDS
#...............................................................................


