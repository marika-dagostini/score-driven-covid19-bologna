#####################################################
###   Time Series - Project (State-Space Model)   ###
#####################################################

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
library(KFAS)
library(numDeriv)
library(SimDesign)
library(Matrix)

# Load Data ===================================================================
rm(list=ls())
set.seed(123)

data_covid = read_excel("dati.xlsx")

y = data_covid$n_cases
dates = data_covid$date

# =============================================================================
# SS - PARAMETER OPTIMIZATION
# =============================================================================

# Function Log-likelihood SS ==================================================

neg_log_likelihood_SS = function(params, return_ft = F) {
  
  # y = (n x 1) matrix containing the observations
  # Z = (p x m x n) system matrix of observation equation
  # T = (m x m x 1) first system matrix of state equation
  # Q = (k x k x 1) covariance matrix of state disturbances eta
  # a1 = (m x 1) matrix containing the expected values of the initial states
  # P1 = (m x m) covariance matrix of the non-diffuse part of 
  #               the initial state vector.
  # u = (n x 1) matrix of an additional parameters in case of non-Gaussian
  #             model ---> dispersion parameter for negative binomial
  
  sigma2_delta = params[1]
  sigma2_beta = params[2]
  sigma2_gamma = params[3]
  v = params[4]
  
  # Design matrix for seasonality component
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
  
  # number of states (9: δt,βt,𝛾_mon,t, 𝛾_tue,t, ..., 𝛾_sun,t)
  m = 9
  n = length(y)
  
  # observation matrix Z
  Zt = array(rep(0, n*1*m), dim=c(1, m, n))
  for(i in 1:n){
    Zt[1, 1, i] = 1
    Zt[1, 3:m, i] = D[i,]
  }
  
  # state transition matrix T
  T.mat = diag(m)
  T.mat[1,2] = 1
  
  # state disturbance matrix R
  R = diag(m)
  
  # covariance matrix of state disturbances Q
  i_7 = matrix(rep(1,7), nrow = 7, ncol = 1)
  Sigma_gamma = sigma2_gamma * (diag(7)-(1/7)*i_7%*%t(i_7))
  
  Q = diag(m)
  Q[1,1] = sigma2_delta
  Q[2,2] = sigma2_beta
  Q[3:m,3:m] = Sigma_gamma
  
  # initial state vector, assuming it starts at 0
  a1 = rep(0, m)
  
  # covariance matrix of the initial state, initialized with some large variance
  P1 = diag(10, m)
  
  model = SSModel(y ~ -1 + SSMcustom(
    Z = Zt, T = T.mat, R = R, Q = Q, a1 = a1, P1 = P1), 
    distribution = "negative binomial", u = v)
  
  fit = fitSSM(model, inits = rep(0, length(params)))


  if(return_ft == F){
    return(-logLik(fit$model))} # negative since we use nlminb(): to find the maximum, 
  else if (return_ft == T){ # we must minimize the negative objective function
    
    model_out = KFS(fit$model, nsim = 1000)
    
    ft = model_out$muhat
    
    plot.ts(y, main = "Filtered estimates of new cases \nof infection with COVID-19 (SS)")
    lines(ft, col = "red", lwd = 1)
    legend("topright", legend = c("Observed", "Filtered estimates"),
           col=c("black", "red"), lty=1)
    return(ft)
  }
}

# Initial Values ==============================================================

iv_sigma2_beta = 0.01
bound_sigma2_beta = c(1e-20, 10)

iv_sigma2_delta = 0.02
bound_sigma2_delta = c(1e-20, 10)

iv_sigma2_gamma = 0.03
bound_sigma2_gamma = c(1e-20, 10)

iv_v = 5 #Dispersion
bound_v = c(1e-20, 1e20)

par_names = c("sigma2_beta", "sigma2_delta", "sigma2_gamma", "v")

initial_values = c(iv_sigma2_beta, iv_sigma2_delta, iv_sigma2_gamma, iv_v)

lower_bounds = c(bound_sigma2_beta[1], bound_sigma2_delta[1],
                 bound_sigma2_gamma[1], bound_v[1])

upper_bounds = c(bound_sigma2_beta[2], bound_sigma2_delta[2],
                 bound_sigma2_gamma[2], bound_v[2])

# Optimization ================================================================
opt_initial = nlminb(start = initial_values, 
                     objective = neg_log_likelihood_SS, 
                     lower = lower_bounds, 
                     upper = upper_bounds)
opt_initial$par

# Optimize again with optimal parameters as initial values ====================

opt2 = nlminb(start = opt_initial$par,
              objective = neg_log_likelihood_SS, 
              lower = lower_bounds, 
              upper = upper_bounds)
opt2$par # Significant changes, perform the optimization another time

opt = nlminb(start = opt2$par, 
             objective = neg_log_likelihood_SS,
             lower = lower_bounds, 
             upper = upper_bounds)
opt$par # Basically no changes

# Significance of estimated optimal parameters ================================

# calculate a numerical approximation to the Hessian matrix 
hessian_matrix = hessian(neg_log_likelihood_SS, opt$par)

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
ft_SS = neg_log_likelihood_SS(param = opt$par, return_ft = T)

# save filtered estimates, parameters and metrics
write_xlsx(x = as.data.frame(ft_SS), 
           path = "SS_filtered_estimates_full_sample_period.xlsx",
           col_names = TRUE)
write_xlsx(x = as.data.frame(signif_results), 
           path = "SS_opt_param_full_sample_period.xlsx",
           col_names = TRUE)
write_xlsx(x = as.data.frame(cbind(metrics, values)), 
           path = "SS_metrics_full_sample_period.xlsx",
           col_names = TRUE)

# =============================================================================
# SS - PREDICTIVE CAPACITY
# =============================================================================

# Function SDWS Forecasting ====================================================

prediction_SS = function(y, dates, optim_params, window){
  
  # initialize metrics to be monitored
  AIC = numeric()
  AICc = numeric()
  BIC = numeric()
  MSE = numeric()
  MAE = numeric()
  MAPE = numeric()
  
  # recall optimal parameters
  sigma2_delta = optim_params[1]
  sigma2_beta = optim_params[2]
  sigma2_gamma = optim_params[3]
  v = optim_params[4]
  
  n = length(y)
  
  # Design matrix for seasonality component
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
  num_replicates = round(n/7,0)
  D = rep(c(t(D)), num_replicates)
  D = matrix(D, ncol = 7, byrow = TRUE)
  
  # number of states (9: δt,βt,𝛾_mon,t, 𝛾_tue,t, ..., 𝛾_sun,t)
  m = 9
  
  # observation matrix Z
  Zt = array(rep(0, n*1*m), dim=c(1, m, n))
  for(i in 1:n){
    Zt[1, 1, i] = 1
    Zt[1, 3:m, i] = D[i,]
  }
  
  # state transition matrix T
  T.mat = diag(m)
  T.mat[1,2] = 1
  
  # state disturbance matrix R
  R = diag(m)
  
  # covariance matrix of state disturbances Q
  i_7 = matrix(rep(1,7), nrow = 7, ncol = 1)
  Sigma_gamma = sigma2_gamma * (diag(7)-(1/7)*i_7%*%t(i_7))
  
  Q = diag(m)
  Q[1,1] = sigma2_delta
  Q[2,2] = sigma2_beta
  Q[3:m,3:m] = Sigma_gamma
  
  # initial state vector, assuming it starts at 0
  a1 = rep(0, m)
  
  # covariance matrix initial states, initialized with some large variance
  P1 = diag(10, m)
  
  for (i in round(length(y)/2-window):(length(y)-window)) {
    y_red = y[1:i]
    
    Zt_red = array(Zt[ , , 1:i], dim = c(1, m, i))
    
    model = SSModel(y_red ~ -1 + SSMcustom(
      Z = Zt_red, T = T.mat, R = R, Q = Q, a1 = a1, P1 = P1), 
      distribution = "negative binomial", u = v)
    
    newdata = SSModel(rep(NA, window) ~ -1 + SSMcustom(
      Z = array(Zt[ , , (i+1):(i+window)], dim = c(1, m, window)), 
      T = T.mat, R = R, Q = Q, a1 = a1, P1 = P1), 
      distribution = "negative binomial", u = v)
    
    pred = predict(model, newdata = newdata,
                   interval = "prediction", level = 0.9, nsim = 100)
    
    n_param = length(optim_params)
    LL = logLik(model) # maximum value of the log-likelihood 
    
    n_red = length(y_red)
    
    #AIC
    AIC_1 = 2*n_param - 2*LL
    AIC = append(AIC, AIC_1)
    
    #AICc
    AICc_1 = AIC_1 + (2*n_param^2 + 2*n_param)/(n_red-n_param-1)
    AICc = append(AICc, AICc_1)
    
    #BIC
    BIC_1 = n_param*log(n_red) - 2*LL
    BIC = append(BIC, BIC_1)
    
    y_f = pred[,1] # forecasts
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
SS_7 = prediction_SS(y, dates, optim_params = opt$par, window = 7)
SS_14 = prediction_SS(y, dates, optim_params = opt$par, window = 14)
SS_28 = prediction_SS(y, dates, optim_params = opt$par, window = 28)

# Save AIC, AICc, BIC, and loss functions =====================================
results_SS = rbind(SS_7, SS_14, SS_28)

write_xlsx(x = results_SS, 
           path = "SS_results_predictive_capacity.xlsx",
           col_names = TRUE)
