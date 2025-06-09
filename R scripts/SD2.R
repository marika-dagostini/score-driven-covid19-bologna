######################################
###   Time Series - Project (SD2)  ###
######################################

# WD ==========================================================================

setwd("C:/Users/computer/Desktop/Corsi/Time Series/Project")

# Libraries ===================================================================
library(readxl)
library(rgdal)
library(spdep)
library(stats4)
library(tidyr)
library(tseries)
library(lubridate)
library(writexl)

# Load Data ===================================================================
rm(list=ls())
set.seed(123)

data_covid = read_excel("dati.xlsx")

y = data_covid$n_cases
dates = data_covid$date

# =============================================================================
# SD2 - PARAMETER OPTIMIZATION
# =============================================================================

# Function Log-likelihood SD2 =================================================

neg_log_likelihood_SD2 = function(param, return_ft = F) {
  
  # kappa parameters
  kappa1 = param[1]
  kappa2 = param[2]
  
  # local level component delta_t
  delta = rep(0,length(y))
  delta[1] = param[3]
  
  # trend component beta_t
  beta = rep(0,length(y))
  beta[1] = param[4]
  
  # shape parameter NB distribution
  v = param[5]
  
  # vector of seasonality filter gamma_t
  gamma = matrix(rep(0,7*length(y)), ncol = 7, byrow = TRUE)
  gamma[1,1] = param[6]
  gamma[1,2] = param[7]
  gamma[1,3] = param[8]
  gamma[1,4] = param[9]
  gamma[1,5] = param[10]
  gamma[1,6] = param[11]
  # sum to zero of the seasonality filter
  gamma[1,7] = 
    -(gamma[1,1]+gamma[1,2]+gamma[1,3]+gamma[1,4]+gamma[1,5]+gamma[1,6])
  
  # kappa_j = kappa for each j=1, ..., 7
  kappa = param[12]
  
  # D matrix
  first_day_week = lubridate::wday(dates[1], week_start = 1)
  D = rep(c(1,0,0,0,0,0,0,0),7)
  D = c(rep(0,first_day_week-1), D)
  D = D[1:(length(D)-first_day_week+1)]
  D = matrix(D, ncol = 7, byrow = T)
  # Get the index of are all zeros row and remove it
  zero_row_index = which(apply(D, 1, function(row) all(row == 0)))
  # obtain D for one week
  D = D[-zero_row_index,]
  # Determine the number of times to replicate the matrix
  num_replicates = round(length(y)/7,0)
  D = rep(c(t(D)), num_replicates)
  D = matrix(D, ncol = 7, byrow = TRUE)
  
  # seasonality component s_t
  s = rep(0,length(y))
  s[1] = D[1,] %*% gamma[1,]
  
  # score driven parameter ft: conditional mean of new COVID-19 cases
  ft = rep(0,length(y))
  ft[1] = exp(delta[1] + s[1])
  
  # scaled score function updating term ut
  u = rep(0,length(y))
  u[1] = y[1]/ft[1] - 1
  
  # log of the NB DGP
  logp = rep(0,length(y))
  logp[1] = dnbinom(x = y[1], mu = ft[1] , size = v, log = TRUE)
  
  # Loop to estimate the SD1 parameters
  for (t in 2:length(y)) {
    delta[t] = delta[t-1] + beta[t-1] + kappa1*u[t-1]
    beta[t]= beta[t-1] + kappa2*u[t-1]
    
    # find j for which Dj,t = 1
    j = which(D[t-1, ] == 1)
    
    # Define kappa_t
    kappa_t = rep(-kappa/(7-1), 7)
    kappa_t[j] = kappa
    
    # update gamma
    gamma[t,] = gamma[t-1,] + kappa_t * u[t-1]
    
    # update s
    s[t] = D[t,] %*% gamma[t,]
    
    # update ft
    ft[t] = exp(delta[t] + s[t])
    
    # update ut
    u[t] = y[t]/ft[t] - 1
    
    # compute the NB density
    logp[t] = dnbinom( x = y[t], mu = ft[t] , size = v, log = TRUE)
  }
  
  if(return_ft == F){
    return(-sum(logp))} # negative since we use nlminb(): to find the maximum, 
  else if (return_ft == T){ # we must minimize the negative objective function
    
    plot.ts(y, main = "Filtered estimates of new cases \nof infection with COVID-19 (SD2)")
    lines(ft, col = "red", lwd = 1)
    legend("topright", legend = c("Observed", "Filtered estimates"),
           col=c("black", "red"), lty=1)
    
    return(ft)
  }
}

# Initial Values ==============================================================

iv_kappa1 = 0.01
bound_kappa1 = c(0,10)

iv_kappa2 = 0.01
bound_kappa2 = c(0,10)

iv_delta1 = log(mean(y))
bound_delta1 = c(-1e20, 1e20)

iv_beta1 = 0 
bound_beta1 = c(-1e20, 1e20) 

iv_v = 5 #Dispersion
bound_v = c(1e-20, 1e20)

iv_gamma_1_1 = 0
bound_gamma_1_1 = c(-1e20, 1e20)

iv_gamma_1_2 = 0
bound_gamma_1_2 = c(-1e20, 1e20)

iv_gamma_1_3 = 0
bound_gamma_1_3 = c(-1e20, 1e20)

iv_gamma_1_4 = 0
bound_gamma_1_4 = c(-1e20, 1e20)

iv_gamma_1_5 = 0
bound_gamma_1_5 = c(-1e20, 1e20)

iv_gamma_1_6 = 0
bound_gamma_1_6 = c(-1e20, 1e20)

iv_kappa = 0.01
bound_kappa = c(0,10)

par_names = c("kappa1", "kappa2", "delta", "beta", "v", 
              "gamma1", "gamma2", "gamma3", 
              "gamma4", "gamma5", "gamma6", "kappa")

initial_values = c(iv_kappa1, iv_kappa2,
                   iv_delta1, iv_beta1, iv_v,
                   iv_gamma_1_1, iv_gamma_1_2, iv_gamma_1_3, 
                   iv_gamma_1_4, iv_gamma_1_5, iv_gamma_1_6,
                   iv_kappa)

lower_bounds = c(bound_kappa1[1], bound_kappa2[1],
                 bound_delta1[1] ,bound_beta1[1], bound_v[1],
                 bound_gamma_1_1[1], bound_gamma_1_2[1], bound_gamma_1_3[1],
                 bound_gamma_1_4[1], bound_gamma_1_5[1], bound_gamma_1_6[1],
                 bound_kappa[1])

upper_bounds = c(bound_kappa1[2], bound_kappa2[2],
                 bound_delta1[2], bound_beta1[2], bound_v[2], 
                 bound_gamma_1_1[2], bound_gamma_1_2[2], bound_gamma_1_3[2],
                 bound_gamma_1_4[2], bound_gamma_1_5[2], bound_gamma_1_6[2],
                 bound_kappa[2])

# Optimization ================================================================
opt_initial = nlminb(objective = neg_log_likelihood_SD2, 
                     start = initial_values,
                     lower = lower_bounds,
                     upper = upper_bounds)

opt_initial$par

# Optimize again with optimal parameters as initial values ====================

opt = nlminb(objective = neg_log_likelihood_SD2, 
             start = opt_initial$par,
             lower = lower_bounds,
             upper = upper_bounds)
opt$par # Basically no changes

# Significance of estimated optimal parameters ================================

# calculate a numerical approximation to the Hessian matrix 
hessian_matrix = hessian(neg_log_likelihood_SD2, opt$par)

# Invert the Hessian matrix to get the covariance matrix:
#   - If the Hessian matrix is not singular, the inverse matrix is computed.
#   - If the Hessian matrix is singular (or computationally singular),
#     the pseudo-inverse matrix is computed instead

cov_matrix = tryCatch(solve(hessian_matrix),
                      error = function(e) return(ginv(hessian_matrix)))

# Get standard errors
std_err = sqrt(diag(cov_matrix))

# Z-scores and p-values for each parameter
z_scores = opt$par / std_err
p_values = 2 * (1 - pnorm(abs(z_scores)))

signif_results = data.frame(Parameter_name = par_names,
                            Parameter = round(opt$par, 4),
                            Std_Error = round(std_err, 4),
                            Z_Score = round(z_scores, 4),
                            P_Value = round(p_values, 4))

# AIC, AICc, and BIC ==========================================================

n_param = length(initial_values)
LL = -opt$objective # maximum value of the likelihood
n = length(y)

AIC = 2*n_param - 2*LL
AICc = AIC + (2*n_param^2 + 2*n_param)/(n-n_param-1)
BIC = n_param*log(n) - 2*LL

metrics = c("LL","AIC","AICc","BIC")
values = c(round(LL,4), round(AIC,4), round(AICc,4), round(BIC,4))
cbind(metrics, values)

# Graph filtered estimates ====================================================
ft_SD2 = neg_log_likelihood_SD2(param = opt$par, return_ft = T)

# save filtered estimates, parameters and metrics
write_xlsx(x = as.data.frame(ft_SD2), 
           path = "SD2_filtered_estimates_full_sample_period.xlsx",
           col_names = TRUE)
write_xlsx(x = as.data.frame(signif_results), 
           path = "SD2_opt_param_full_sample_period.xlsx",
           col_names = TRUE)
write_xlsx(x = as.data.frame(cbind(metrics, values)), 
           path = "SD2_metrics_full_sample_period.xlsx",
           col_names = TRUE)

# =============================================================================
# SD2 - PREDICTIVE CAPACITY
# =============================================================================

# Function SD2 Forecasting ====================================================

prediction_SD2 = function(y, dates, optim_params, window){
  
  # initialize metrics to be monitored
  AIC = numeric()
  AICc = numeric()
  BIC = numeric()
  MSE = numeric()
  MAE = numeric()
  MAPE = numeric()
  
  for (i in round(length(y)/2-window):(length(y)-window)) {
    y_red = y[1:i]
    y_red = append(y_red, rep(NA, window))
    T_ = length(y_red)
    
    # kappa parameters
    kappa1 = optim_params[1]
    kappa2 = optim_params[2]
    
    # local level component delta_t
    delta = rep(0,T_)
    delta[1] = optim_params[3]
    
    # trend component beta_t
    beta = rep(0,T_)
    beta[1] = optim_params[4]
    
    # shape parameter NB distribution
    v = optim_params[5]
    
    # vector of seasonality filter gamma_t
    gamma = matrix(rep(0,7*T_), ncol = 7, byrow = TRUE)
    gamma[1,1] = optim_params[6]
    gamma[1,2] = optim_params[7]
    gamma[1,3] = optim_params[8]
    gamma[1,4] = optim_params[9]
    gamma[1,5] = optim_params[10]
    gamma[1,6] = optim_params[11]
    # sum to zero of the seasonality filter
    gamma[1,7] = 
      -(gamma[1,1]+gamma[1,2]+gamma[1,3]+gamma[1,4]+gamma[1,5]+gamma[1,6])
    
    # vector kappa_t
    kappa = optim_params[12]
    
    # D matrix
    first_day_week = lubridate::wday(dates[1], week_start = 1)
    D = rep(c(1,0,0,0,0,0,0,0),7)
    D = c(rep(0,first_day_week-1), D)
    D = D[1:(length(D)-first_day_week+1)]
    D = matrix(D, ncol = 7, byrow = T)
    # Get the index of are all zeros row and remove it
    zero_row_index = which(apply(D, 1, function(row) all(row == 0)))
    # obtain D for one week
    D = D[-zero_row_index,]
    # Determine the number of times to replicate the matrix
    num_replicates = round(T_/7+1,0)
    D = rep(c(t(D)), num_replicates)
    D = matrix(D, ncol = 7, byrow = TRUE)
    
    # seasonality component s_t
    s = rep(0,T_)
    s[1] = D[1,] %*% gamma[1,]
    
    # score driven parameter ft: conditional mean of new COVID-19 cases
    ft = rep(0,T_)
    ft[1] = exp(delta[1] + s[1])
    
    # scaled score function updating term ut
    u = rep(0,T_)
    u[1] = y[1]/ft[1] - 1
    
    # log of the NB DGP
    logp = rep(0,T_)
    logp[1] = dnbinom(x = y[1], mu = ft[1] , size = v, log = TRUE)
    
    # Forecasting Loop
    for (t in 2:T_) {
      delta[t] = delta[t-1] + beta[t-1] + kappa1*u[t-1]
      beta[t]= beta[t-1] + kappa2*u[t-1]
      
      # find j for which Dj,t = 1
      j = which(D[t-1, ] == 1)
      
      # Define kappa_t
      kappa_t = rep(-kappa/(7-1), 7)
      kappa_t[j] = kappa
      
      # update gamma
      gamma[t,] = gamma[t-1,] + kappa_t * u[t-1]
      
      # update s
      s[t] = D[t,] %*% gamma[t,]
      
      # update ft
      ft[t] = exp(delta[t] + s[t])
      
      # save prediction
      if (is.na(y_red[t])){y_red[t]=ft[t]}
      
      #update log-likelihood
      logp[t] = dnbinom(x = round(y_red[t]), mu = round(ft[t]) , 
                        size = v, log = TRUE)
      
      # update ut
      u[t] = y_red[t]/ft[t] - 1
    }
    
    n_param = length(optim_params)
    LL = sum(logp[1:i]) # maximum value of the log-likelihood 
    # (excluding the prediction window)
    n = length(y[1:i])
    
    #AIC
    AIC_1 = 2*n_param - 2*LL
    AIC = append(AIC, AIC_1)
    
    #AICc
    AICc_1 = AIC_1 + (2*n_param^2 + 2*n_param)/(n-n_param-1)
    AICc = append(AICc, AICc_1)
    
    #BIC
    BIC_1 = n_param*log(n) - 2*LL
    BIC = append(BIC, BIC_1)
    
    y_f = y_red[(length(y_red)-(window-1)):length(y_red)] # forecasts
    y_true = y[(i+1):(i+window)] # true values
    
    #MSE
    MSE_1 = mean((y_f-y_true)^2)
    MSE = append(MSE, MSE_1)
    
    #MAE
    MAE_1 = mean(abs(y_f-y_true))
    MAE = append(MAE, MAE_1)
    
    #MAPE
    MAPE_1 = mean(abs(y_true-y_f)/y_true)
    MAPE = append(MAPE, MAPE_1)
    
  }
  
  metrics = c("Forecasting Window", "AIC", "AICc", "BIC", "MSE", "MAE", "MAPE")
  values = c(window, mean(AIC), mean(AICc), mean(BIC), 
             mean(MSE), mean(MAE), mean(MAPE))
  df_results = as.data.frame(t(values))
  colnames(df_results) = metrics
  return(df_results)
}

# 7-, 14-, and 28-days ahead forecast metrics =================================
SD2_7 = prediction_SD2(y, dates, optim_params=opt$par, window=7)
SD2_14 = prediction_SD2(y, dates, optim_params=opt$par, window=14)
SD2_28 = prediction_SD2(y, dates, optim_params=opt$par, window=28)

# Save AIC, AICc, BIC, and loss functions =====================================
results_SD2 = rbind(SD2_7, SD2_14, SD2_28)

write_xlsx(x = results_SD2, 
           path = "SD2_results_predictive_capacity.xlsx",
           col_names = TRUE)
