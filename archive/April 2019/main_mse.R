# Todos
#   -Add/discuss constraints


# Purpose of tutorial: illustrate DFA (univariate), MDFA (multivariate) MSE (no customization) unconstrained (no regularization)
#   -Customization and regularization will be tackled in separate tutorials
# 1. Provide short overview of main functions in MDFA-package
# 2. Play with (M)DFa
#     -Illustrate explain target, spectrum, amplitude, time-shift, filter coefficients
#     -Illustrate overfitting: univariate and multivariate
# 3. Introduce/discuss potentially useful filter constraints

# Disclaimer/caveat: applications to (currency-)trading are intended for illustrative purposes only 
#   -Filter designs are deliberately 'suboptimal'


rm(list=ls())

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)


# Briev overview of wrappers and main function
head(MDFA_mse)
head(MDFA_mse_constraint)
head(MDFA_cust)
head(MDFA_cust_constraint)
head(MDFA_reg)
head(MDFA_reg_constraint)
# Main estimation function
head(mdfa_analytic)

#-----------------------------------------------------------------------------------------------
# Source common functions

source("Common functions/plot_func.r")
source("Common functions/mdfa_trade_func.r")
source("Common functions/data_load_functions.R")


#-----------------------------------------------------------------------------------------------
# Data: FX 
data_from_IB<-T
hour_of_day<-"16:00"
# Select most liquid pairs: all 6 pairs with USD, EUR, GBP and JPY
i_series_vec<-c(1,2,3,6,7,8)
reload_sp500<-F
path.dat<-paste(getwd(),"/",sep="")

data_load_obj<-data_load_gzd_trading_func(data_from_IB,hour_of_day,reload_sp500,path.dat)

mydata<-data_load_obj$mydata

FX_mat<-data_load_obj$xts_data_mat[,i_series_vec]
log_FX_mat<-log(FX_mat)

colnames(log_FX_mat)

plot_T<-T
anf_plot<-"2000-10-01/"
#---------------------------------------------
# Example 0: Play with DFA
#   MSE
#   Univariate


# Example 0.1 White noise spectrum
# Try
#   -various targets (periodicities); 
#   -various filter-lengths L; 
#   -various Lag (negative,0,positive integers)
# K defines the frequency-grid for estimation: all frequencies omega_j=j*pi/K, j=0,1,2,...,K are considered
#   -In the case of a white noise spectrum we can select arbitrary tightness of frequency grid
#   -If we use the periodogram, then K is specified by the length of the periodogram
#   -Tradeoff: larger K (denser grid) means better MSE (if true process is white noise) but longer computation-time
K<-600
# Spectrum: white noise assumption
#  First column is target, second column is explanatory variable: in a univariate design target and explanatory are the same
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
colnames(weight_func)<-c("target","explanatory")
# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-200

# Estimation based on MDFA-MSE wrapper
mdfa_obj_mse<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 

b<-mdfa_obj_mse$b
plot_estimate_func(mdfa_obj_mse,weight_func,Gamma)

#---------------
# Example 0.2 same as above but with autoregressive spectrum
#   This design replicates classic model-based approaches (if K is 'large')
# Try
#   -various targets (periodicities); 
#   -various filter-lengths L; 
#   -various Lag (negative,0,positive integers)
# K defines the frequency-grid for estimation: all frequencies omega_j=j*pi/K, j=0,1,2,...,K are considered
#   -In the case of a white noise spectrum we can select arbitrary tightness of frequency grid
#   -If we use the periodogram, then K is specified by the length of the periodogram
#   -Tradeoff: larger K (denser grid) means better MSE (if true process is white noise) but longer computation-time
K<-600
# Spectrum: autoregressive spectrum
#  First column is target, second column is explanatory variable: in a univariate design target and explanatory are the same
a1<-0.5
weight_func<-cbind(1/(1-a1*exp(1.i*(0:K)*pi/K)),1/(1-a1*exp(1.i*(0:K)*pi/K)))
colnames(weight_func)<-c("target","explanatory")
# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-200

# Estimation based on MDFA-MSE wrapper
mdfa_obj_mse<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 

b<-mdfa_obj_mse$b
plot_estimate_func(mdfa_obj_mse,weight_func,Gamma)

# Comments
# 1. Fit of target by explanatory improves at/towards the frequencies emphasized/weighted by spectrum
#     Examples: 
#       for a1>0 the fit improves towards frequency zero (compared to previous white noise example)
#       for a1<0 the fit improves towards frequency pi (compared to previous white noise example)
# 2. Arbitrary AR(I)MA- or state-space models could be replicated by 
#   a. Prodining the corresponding target
#   b. Providing the corresponding spectral estimate
#   c. see Wildi/McElroy
# 3. Once replicated, model-based approaches could be enhanced (regularization and/or customization)
# 4. Multivariate models can be replicated by generic MDFA

#---------------------------------
# Example 0.3 same as above but we use autoregressive spectrum based on AR(1)-model fitted to log-returns of EURUSD

K<-600
# Fit AR(1)-model to log-returns of EURUSD
asset<-"EURUSD"
x<-diff(log_FX_mat[,asset])
a1<-arima(x,c(1,0,0))$coef["ar1"]
# Derive autoregressive spectrum
weight_func<-cbind(1/(1-a1*exp(1.i*(0:K)*pi/K)),1/(1-a1*exp(1.i*(0:K)*pi/K)))
colnames(weight_func)<-c("target","explanatory")
# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-200

# Estimation based on MDFA-MSE wrapper
mdfa_obj_mse<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 

b<-mdfa_obj_mse$b
plot_estimate_func(mdfa_obj_mse,weight_func,Gamma)

#-------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

# The following examples apply (M)DFA to currency-trading
#   The function mdfa_mse_reg_trade_func 
#     -computes optimal MSE-filters
#     -computes filter-outputs 
#     -applies filters to trading according to the rule: long/short depending on sign of filter output
#     -computes annualized sharpe ratio

#-------------------------------------------------
#-------------------------------------------------
# Example 1: DFA-lowpass applied to log-returns; white noise spectrum
#   MSE
#   Univariate
# Try various targets (periodicities)
asset<-"EURUSD"
x<-na.exclude(diff(log_FX_mat[,asset]))
# K defines the frequency-grid for estimation: all frequencies omega_j=j*pi/K, j=0,1,2,...,K are considered
# In the case of a white noise spectrum we can select arbitrary tightness of frequency grid
#   If we use the periodogram, then K is specified by the length of the periodogram
# Tradeoff: larger K (denser grid) means better MSE (if true process is white noise) but longer computation-time
K<-600
# Spectrum: white noise assumption
#  First column is target, second column is explanatory variable: in a univariate design target and explanatory are the same
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
colnames(weight_func)<-c("target","explanatory")
# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity<-5
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-200
# Ensure that number of degrees of freedom does not exceed number of available equations (frequencies)
#   Otherwise problem is not feasible (matrix cannot be inverted)
L<-min(2*K,L)
# Select FX-series
asset<-"EURUSD"
# Lag trade execution by at least one day (must be >=1)
lag_fx<-1
# Allow for some plots in the function below
plot_T<-T

# This function is in mdfa_mse_reg_trade_func.r
# MSE estimation and trading: sign rule and signal weighting
mdfa_mse_reg_trade_obj<-mdfa_mse_reg_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func)

#-------------------------------------------------
# Example 2: same as above but we use periodogram (instead of white noise spectrum)
# Design
#   MSE
#   Univariate
#   Periodogram

# Specify in-sample span fpr estimation of spectrum (periodogram)
in_sample_span<-"2017-01-01"
# Select FX-series
asset<-"EURUSD"

x<-na.exclude(diff(log_FX_mat[,asset]))
# Spectrum: 
#   MDFA relies on dft (periodogram is abs(dft)^2)
#   In contrast to periodogram, the dft retains the phase information which is important in a multivariate design
#   First column is dft of target, second column is dft of explanatory: here both series are identical (univariate design)
weight_func<-cbind(per(x[paste("/",in_sample_span,sep="")],T)$DFT,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
colnames(weight_func)<-c("target","explanatory")
# Denseness of frequency-grid is completely specified by length of periodogram (to be constrasted with example 1 where we could specify any K)
K<-nrow(weight_func)-1
# Target
periodicity<-5
# Nowcast
Lag<-0
# Filter length
L<-100
L<-min(2*K,L)
# Delay of trade execution
lag_fx<-1
# Some plots
plot_T<-T

# This function is in mdfa_mse_reg_trade_func.r
# MSE estimation and trading: sign rule and signal weighting
mdfa_mse_reg_trade_obj<-mdfa_mse_reg_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func)


#-----------------------------------------------------------
# Example 3: same as above; play with the filter length L
#   -Since we use the dft large L imply overfitting
#   -Evaluate the extent of overfitting
#   -Select large L (for instance L=900) or 'reasonably large' L (for instance L<-2*periodicity) and observe in-sample vs. out-of-sample performances
# We use all series
in_sample_span<-"2018-01-01"

for (i_series in 1:ncol(log_FX_mat))
{

  
  x<-na.exclude(diff(log_FX_mat[,i_series]))
# Spectrum: periodogram
  weight_func<-cbind(per(x[paste("/",in_sample_span,sep="")],T)$DFT,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
  colnames(weight_func)<-rep(colnames(log_FX_mat)[i_series],2)
  colnames(weight_func)<-c("target","explanatory")
# Length of frequency-grid (always length of spectral estimate)  
  K<-nrow(weight_func)-1
# Target (cutoff, periodicity)
  periodicity<-5
# Nowcast  
  Lag<-0
# Huge L  
  L<-900
# Reasonably large L  
  L<-2*periodicity
# Degrees of freedom must be smaller than number of equations (otherwise problem is singular)  
  L<-min(2*K,L)#L<-900
# Delay trade execution (must be >=1)  
  lag_fx<-1
  plot_T<-T
  
# This function is in mdfa_trade_func.r
# MSE estimation and trading
  
  mdfa_mse_reg_trade_obj<-mdfa_mse_reg_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func)
  
  cum_perf_sign<-mdfa_mse_reg_trade_obj$cum_perf_sign
  
  diff_perf<-diff(cum_perf_sign)
# Save trading performances  
  if (i_series==1)
  {
    diff_perf_mat<-diff_perf  
  } else
  {
    diff_perf_mat<-cbind(diff_perf_mat,diff_perf)
    
  }
  
}

# in sample performances
plot(as.xts(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum)),main="Individual FX-series in-sample")
plot(as.xts(apply(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum),1,mean)),main="Aggregate in-sample")


# out-of-sample performances
plot(as.xts(apply(diff_perf_mat[paste(in_sample_span,"/",sep="")],2,cumsum)),main="Individual FX-series, out-of-sample")
plot(as.xts(apply(apply(diff_perf_mat[paste(in_sample_span,"/",sep="")],2,cumsum),1,mean)),main="Aggregate out-of-sample")



#-----------------------------------------------------------
# Example 4: multivariate (MDFA) with dft (discrete fourier transform)
# We use all series as explanatory variables: 6-dimensional design (this is a bit unwise but we use the example for illustration)
#   Number of degrees of freedom is 6*L   
# Huge overfitting: try L<-200 with periodicity=5 against 'reasonably large' L:
#   1. Coefficients are huge
#   2. Amplitudes are huge and very narrow i.e. far from pi/5
#   3. However: scale of aggregate (output of multivariate filter) is fine and
#     acf of aggregate is cyclical with periodicity 5
# We conclude that wrong scaling cancels accross series and periodicity of 5 is obtained
#   again by cross-sectional cancelling of low-frequency trends
# Large overfitting: don't feel comfortable with such a design (just 'fitting')

# 4.1 Specify in-sample span for estimation
in_sample_span<-"2018-01-01"

# and compute multivariate spectrum
for (i_series in 1:ncol(log_FX_mat))
{
  x<-na.exclude(diff(log_FX_mat[,i_series]))
  # Spectrum: white noise assumption
  
  if (i_series==1)
  {
    weight_func_mat<-per(x[paste("/",in_sample_span,sep="")],T)$DFT  
  } else
  {
    weight_func_mat<-cbind(weight_func_mat,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
  }
}


# 4.2 For each target (FX-series): use all other (FX-) series as explanatory variables: 6-dim designs
#   Loop through each target (FX_series) and compute filtering and trading: collect performances in a matrix

for (i_series in 1:ncol(log_FX_mat))# i_series<-1
{
  
# Data for filtering: 
#   We re-order series such that first column is target variable
#   The columns of the data to be filtered (data_filt) must match the corresponding columns of the spectrum (weight_func) 
#     Otherwise filters and data don't match    
#   Note that in principle we could skip this re-ordering (the only constraint is that x and weight_func must correspond)
#   But we use the re-ordering for the trading-part in mdfa_mse_reg_trade_func: the first column is always the target that must be traded
  series_ordering<-c(i_series,(1:ncol(log_FX_mat))[-i_series])
# Data matrix for filtering: all explanatory data  
  data_filt<-na.exclude(diff(log_FX_mat[,series_ordering]))
# Spectrum: first series is target, all other series are explanatory
#   Ordering of explanatory series must correspond to data_filt (data for filtering)  
  weight_func<-cbind(weight_func_mat[,i_series],weight_func_mat[,series_ordering])
  colnames(weight_func)<-c(colnames(log_FX_mat)[i_series],colnames(log_FX_mat)[series_ordering])
# Length of sepctral grid  
  K<-nrow(weight_func)-1
# Target (cutoff/periodicity)  
  periodicity<-5
# Nowcast  
  Lag<-0
# Large L
  L<-200
# Reasonably large L  
  L<-2*periodicity#
  L<-min(2*K,L)
# Lag of trading execution
  lag_fx<-1
  plot_T<-T
  
  # This function is in mdfa_mse_reg_trade_func.r
  # MSE estimation and trading: sign rule and signal weighting
  
  mdfa_mse_reg_trade_obj<-mdfa_mse_reg_trade_func(K,periodicity,L,Lag,lag_fx,data_filt,plot_T,weight_func)
  
  cum_perf_sign<-mdfa_mse_reg_trade_obj$cum_perf_sign
  yhat<-mdfa_mse_reg_trade_obj$yhat
# Skip first L data points (initialization)  
  yhat<-yhat[L:length(yhat)]
  diff_perf<-diff(cum_perf_sign)
  
  if (i_series==1)
  {
    diff_perf_mat<-diff_perf 
    yhat_mat<-yhat
  } else
  {
    diff_perf_mat<-cbind(diff_perf_mat,diff_perf)
    yhat_mat<-cbind(yhat_mat,yhat)
    
  }
  
}

anf_plot<-paste(in_sample_span,"/",sep="")
anf_plot<-"1918-01-01/"

# Check individual filter outputs as well as aggregate (the latter is the output of MDFA)
plot(yhat_mat[anf_plot],main="Filter outputs: aggregate in bold")
lines(as.xts(apply(yhat_mat[anf_plot],1,mean)),lwd=6)

# Check periodicity of multivariate filter output: should coincide with parameter periodicity set above
acf(as.xts(apply(yhat_mat[anf_plot],1,mean)))





# Plot performances
# in sample
plot(as.xts(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum)),main="In sample")
plot(as.xts(apply(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum),1,mean)),main="Aggregate in-sample")


# out-of-sample
plot(as.xts(apply(diff_perf_mat[paste(in_sample_span,"/",sep="")],2,cumsum)),main="Out-of-sample")
plot(as.xts(apply(apply(diff_perf_mat[paste(in_sample_span,"/",sep="")],2,cumsum),1,mean)),main="Aggregate out-of-sample")



