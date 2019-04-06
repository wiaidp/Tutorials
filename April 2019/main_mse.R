# Todos
#   -Add a model-based example
#   -Amplitude, shift, coeff, spectrum, target
#   -Explain I-MDFA.r: functions available
#     Here we use MSE function 



rm(list=ls())

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)


#-----------------------------------------------------------------------------------------------
# Source common functions

source("Common functions/plot_func.r")
source("Common functions/mdfa_trade_func.r")
source("Common functions/data_load_functions.R")


#-----------------------------------------------------------------------------------------------
# Data: FX and SP500 (the latter has a marked trend)

data_from_IB<-T
hour_of_day<-"16:00"
i_series_vec<-c(1,2,3,6,7,8)
reload_sp500<-F
path.dat<-"C:\\wia_desktop\\2019\\Projekte\\IB\\daily_pick\\Data\\IB\\"

data_load_obj<-data_load_gzd_trading_func(data_from_IB,hour_of_day,reload_sp500,path.dat)

mydata<-data_load_obj$mydata

FX_mat<-data_load_obj$xts_data_mat[,i_series_vec]
log_FX_mat<-log(FX_mat)

colnames(log_FX_mat)

plot_T<-T
anf_plot<-"2000-10-01/"
# Select data: FX or SP500

#-------------------------------------------------
# Example 1: DFA-lowpass applied to log-returns
#   MSE
#   Univariate
#   White noise spectrum
# Try various targets (periodicities)
asset<-"EURUSD"
x<-na.exclude(diff(log_FX_mat[,asset]))
K<-600
# Spectrum: white noise assumption
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
#weight_func[1,]<-1000
periodicity<-5
Lag<--1
L<-200
L<-min(2*K,L)
plot_T<-F
asset<-"EURUSD"
lag_fx<-1
plot_T<-T

# This function is in mdfa_mse_reg_trade_func.r
# MSE estimation and trading: sign rule and signal weighting
mdfa_mse_reg_trade_obj<-mdfa_mse_reg_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func)

#-------------------------------------------------
# Example 2: same as above but we use periodogram (instead of white noise spectrum)
#   MSE
#   Univariate
#   Periodogram

in_sample_span<-"2017-01-01"
i_series<-1

x<-na.exclude(diff(log_FX_mat[,i_series]))
# Spectrum: periodogram
weight_func<-cbind(per(x[paste("/",in_sample_span,sep="")],T)$DFT,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
K<-nrow(weight_func)-1
#weight_func[1,]<-1000
periodicity<-5
Lag<-0
L<-100
L<-min(2*K,L)
plot_T<-F
lag_fx<-1
plot_T<-T

# This function is in mdfa_mse_reg_trade_func.r
# MSE estimation and trading: sign rule and signal weighting
mdfa_mse_reg_trade_obj<-mdfa_mse_reg_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func)


#-----------------------------------------------------------
# Example 3: same as above and playing with the filter length L
#   -Since we use the periodogram large L imply overfitting
#   -Learning to evaluate the extent of overfitting
#   -Select large L (for instance L=900) and 'reasonably large' L (for instance L<-2*periodicity) and observe in-sample vs. out-of-sample performances
# We use all series
in_sample_span<-"2018-01-01"

for (i_series in 1:ncol(log_FX_mat))
{

  
  x<-na.exclude(diff(log_FX_mat[,i_series]))
# Spectrum: periodogram
  weight_func<-cbind(per(x[paste("/",in_sample_span,sep="")],T)$DFT,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
  colnames(weight_func)<-rep(colnames(log_FX_mat)[i_series],2)
# Length of frequency-grid (always length of spectral estimate)  
  K<-nrow(weight_func)-1
# Target (cutoff, periodicity)
  periodicity<-5
# Nowcast  
  Lag<-0
# Reasonably large L  
  L<-2*periodicity
# Huge L  
  L<-900
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
#  L<-2*periodicity#
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
plot(as.xts(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum)))
plot(as.xts(apply(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum),1,mean)))


# out-of-sample
plot(as.xts(apply(diff_perf_mat[paste(in_sample_span,"/",sep="")],2,cumsum)))
plot(as.xts(apply(apply(diff_perf_mat[paste(in_sample_span,"/",sep="")],2,cumsum),1,mean)))



