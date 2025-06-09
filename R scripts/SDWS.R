########################################
###   Time Series - Project (SDWS)   ###
########################################

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
library(numDeriv)

# Load Data ===================================================================
rm(list=ls())
set.seed(123)

data_covid = read_excel("dati.xlsx")

y = data_covid$n_cases
dates = data_covid$date

# =============================================================================
# SDWS - PARAMETER OPTIMIZATION
# =============================================================================

# Function Log-likelihood SDWS ================================================

neg_log_likelihood_SDWS = function(params, return_ft = F) {
  
  # kappa parameters
  kappa1 = params[1]
  kappa2 = params[2]
  
  # local level component delta_t
  delta = rep(0,length(y))
  delta[1] = params[3]
  
  # trend component beta_t
  beta = rep(0,length(y))
  beta[1] = params[4]
  
  # shape parameter NB distribution
  v = params[5]
  
  # score driven parameter ft: conditional mean of new COVID-19 cases
  ft = rep(0,length(y))
  ft[1] = exp(delta[1])
  
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
    
    # update ft
    ft[t] = exp(delta[t])
    
    # update ut
    u[t] = y[t]/ft[t] - 1
    
    # compute the NB density
    logp[t] = dnbinom( x = y[t], mu = ft[t] , size = v, log = TRUE)
  }
  if(return_ft == F){
    return(-sum(logp))} # negative since we use nlminb(): to find the maximum, 
  else if (return_ft == T){ # we must minimize the negative objective function
    
    plot.ts(y, main = "Filtered estimates of new cases \nof infection with COVID-19 (SDWS)")
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

par_names = c("kappa1", "kappa2", "delta", "beta", "v")

initial_values = c(iv_kappa1, iv_kappa2,
                   iv_delta1, iv_beta1, iv_v)

lower_bounds = c(bound_kappa1[1], bound_kappa2[1],
                 bound_delta1[1] ,bound_beta1[1], bound_v[1])

upper_bounds = c(bound_kappa1[2], bound_kappa2[2],
                 bound_delta1[2], bound_beta1[2], bound_v[2])

# Optimization ================================================================
opt_initial = nlminb(start = initial_values, 
                     objective = neg_log_likelihood_SDWS, 
                     lower = lower_bounds,
                     upper = upper_bounds)

opt_initial$par

# Optimize again with optimal parameters as initial values ====================
opt = nlminb(start = opt_initial$par, 
             objective = neg_log_likelihood_SDWS, 
             lower = lower_bounds,
             upper = upper_bounds)
opt$par # Basically no changes

# Significance of estimated optimal parameters ================================

# calculate a numerical approximation to the Hessian matrix 
hessian_matrix = hessian(neg_log_likelihood_SDWS, opt$par)

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
ft_SDWS = neg_log_likelihood_SDWS(param = opt$par, return_ft = T)

# save filtered estimates, parameters and metrics
write_xlsx(x = as.data.frame(ft_SDWS), 
           path = "SDWS_filtered_estimates_full_sample_period.xlsx",
           col_names = TRUE)
write_xlsx(x = as.data.frame(signif_results), 
           path = "SDWS_opt_param_full_sample_period.xlsx",
           col_names = TRUE)
write_xlsx(x = as.data.frame(cbind(metrics, values)), 
           path = "SDWS_metrics_full_sample_period.xlsx",
           col_names = TRUE)

# =============================================================================
# SDWS - PREDICTIVE CAPACITY
# =============================================================================

# Function SDWS Forecasting ====================================================

prediction_SDWS = function(y, dates, optim_params, window){
  
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
    
    # score driven parameter ft: conditional mean of new COVID-19 cases
    ft = rep(0,T_)
    ft[1] = exp(delta[1])
    
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
      
      # update ft
      ft[t] = exp(delta[t])
      
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
SDWS_7 = prediction_SDWS(y, dates, optim_params = opt$par, window = 7)
SDWS_14 = prediction_SDWS(y, dates, optim_params = opt$par, window = 14)
SDWS_28 = prediction_SDWS(y, dates, optim_params = opt$par, window = 28)

# Save AIC, AICc, BIC, and loss functions =====================================
results_SDWS = rbind(SDWS_7, SDWS_14, SDWS_28)

write_xlsx(x = results_SDWS, 
           path = "SDWS_results_predictive_capacity.xlsx",
           col_names = TRUE)
