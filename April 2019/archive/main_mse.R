# Todos
#   -Remove weekday
#   -Remove cross
#   -Add a model-based example

rm(list=ls())

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)




#-----------------------------------------------------------------------------------------------
# Data: FX and SP500 (the latter has a marked trend)
source("data_load_functions.R")

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



#-----------------------------------------------------------------------------------------------
source("plot_func.r")
source("mdfa_mse_trade_func.r")
#----------------------------------------------------------------------------------------------
# DFA-cross: applied to log-FX
# MSE
# Univariate
# White noise spectrum

K<-2400
# Spectrum: white noise assumption
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
#weight_func[1,]<-1000
periodicity_short<-2.5
periodicity_long<-20
Lag<-0
L<-200
L<-min(2*K,L)
plot_T<-F
asset<-"EURUSD"
lag_fx<-1
x<-na.exclude(log_FX_mat[,asset])

plot_T<-T

ma_cross_mse_obj<-ma_cross_mdfa_mse_trade_func(K,periodicity_short,periodicity_long,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)
  
#-------------------------------------------------
# DFA-lowpass: applied to log-returns
# MSE
# Univariate
# White noise spectrum

K<-600
# Spectrum: white noise assumption
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
#weight_func[1,]<-1000
periodicity<-2.5
Lag<--1
L<-20
L<-min(2*K,L)
sign_day_of_week<-rep(1,5)
#sign_day_of_week[5]<--1
plot_T<-F
lag_fx<-1
plot_T<-T

for (i_series in 1:ncol(log_FX_mat))#i_series<-2
{
  x<-na.exclude(diff(log_FX_mat[,i_series]))
  
# This function is in mdfa_mse_trade_func.r
# MSE estimation and trading: sign rule and signal weighting
  mdfa_mse_trade_obj<-mdfa_mse_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)
  
  cum_perf_sign<-mdfa_mse_trade_obj$cum_perf_sign
  
  if (i_series==1)
  {
    diff_perf_mat<-diff(cum_perf_sign)
  } else
  {
    diff_perf_mat<-cbind(diff_perf_mat,diff(cum_perf_sign))
    
  }
}
# in sample
plot(as.xts(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum)))
plot(as.xts(apply(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum),1,mean)))


# out-of-sample
plot(as.xts(apply(na.exclude(diff_perf_mat[paste(in_sample_span,"/",sep="")]),2,cumsum)))
plot(as.xts(apply(apply(na.exclude(diff_perf_mat[paste(in_sample_span,"/",sep="")]),2,cumsum),1,mean)))

if (F)
  acf(na.exclude(diff(as.xts(apply(apply(na.exclude(diff_perf_mat[paste(in_sample_span,"/",sep="")]),2,cumsum),1,mean)))))

#-------------------------------------------------
# DFA-lowpass: applied to log-returns
# MSE
# Univariate
# Periodogram

in_sample_span<-"2017-01-01"

#weight_func[1,]<-1000
periodicity<-2.5
Lag<--1
L<-20
L<-min(2*K,L)
sign_day_of_week<-rep(1,5)
#sign_day_of_week[5]<--1
plot_T<-F
lag_fx<-1
plot_T<-T

for (i_series in 1:ncol(log_FX_mat))
{
  x<-na.exclude(diff(log_FX_mat[,i_series]))
  # Spectrum: periodogram
  weight_func<-cbind(per(x[paste("/",in_sample_span,sep="")],T)$DFT,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
  K<-nrow(weight_func)-1
  
# This function is in mdfa_mse_trade_func.r
# MSE estimation and trading: sign rule and signal weighting
  mdfa_mse_trade_obj<-mdfa_mse_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)

  cum_perf_sign<-mdfa_mse_trade_obj$cum_perf_sign

  if (i_series==1)
  {
    diff_perf_mat<-diff(cum_perf_sign)
  } else
  {
    diff_perf_mat<-cbind(diff_perf_mat,diff(cum_perf_sign))
    
  }
}
# in sample
plot(as.xts(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum)))
plot(as.xts(apply(apply(na.exclude(diff_perf_mat[paste("/",in_sample_span,sep="")]),2,cumsum),1,mean)))


# out-of-sample
plot(as.xts(apply(na.exclude(diff_perf_mat[paste(in_sample_span,"/",sep="")]),2,cumsum)))
plot(as.xts(apply(apply(na.exclude(diff_perf_mat[paste(in_sample_span,"/",sep="")]),2,cumsum),1,mean)))


#-----------------------------------------------------------
# Univariate with periodogram
# Learning to evaluate the extent of overfitting
# Select large L and observe in-sample vs. out-of-sample performances
# All series
in_sample_span<-"2017-01-01"

for (i_series in 1:ncol(log_FX_mat))
{

  
  x<-na.exclude(diff(log_FX_mat[,i_series]))
  # Spectrum: periodogram
  weight_func<-cbind(per(x[paste("/",in_sample_span,sep="")],T)$DFT,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
  colnames(weight_func)<-rep(colnames(log_FX_mat)[i_series],2)
  
  K<-nrow(weight_func)-1
  #weight_func[1,]<-1000
  periodicity<-2.5
  Lag<--1
  L<-2*periodicity
  L<-min(2*K,L)#L<-900
  sign_day_of_week<-rep(1,5)
#  sign_day_of_week[5]<--1
  plot_T<-F
  lag_fx<-1
  plot_T<-T
  
  # This function is in mdfa_mse_trade_func.r
  # MSE estimation and trading: sign rule and signal weighting
  mdfa_mse_trade_obj<-mdfa_mse_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)
  
  cum_perf_sign<-mdfa_mse_trade_obj$cum_perf_sign
  
#  diff_perf<-diff(cum_perf_sign)
  
  if (i_series==1)
  {
    perf_mat<-cum_perf_sign  
  } else
  {
    perf_mat<-cbind(perf_mat,cum_perf_sign)
    
  }
  
}

# in sample
plot(as.xts(apply(perf_mat[paste("/",in_sample_span,sep="")])))
plot(as.xts(apply(perf_mat[paste("/",in_sample_span,sep="")],1,mean)))


# out-of-sample
plot(as.xts(perf_mat[paste(in_sample_span,"/",sep="")]))
plot(as.xts(apply(perf_mat[paste(in_sample_span,"/",sep="")],1,mean)))



#-----------------------------------------------------------
# Multivariate with dft
# All series
# Huge overfitting: try L<-200 with periodicity=5:
#   1. Coefficients are huge
#   2. Amplitudes are huge and very narrow i.e. far from pi/5
#   3. However: scale of aggregate (output of multivariate filter) is fine and
#     acf of aggregate is cyclical with periodicity 5
# We conclude that wrong scaling cancels accross series and periodicity of 5 is obtained
#   again by cross-sectional cancelling of low-frequency trends
in_sample_span<-"2017-01-01"
# Compute multivariate spectrum
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



for (i_series in 1:ncol(log_FX_mat))# i_series<-1
{
  
# Data for filtering: 
#   We re-order series such that first column is target variable
#   The ordering of the data to be filtered (x) and of the xpectrum (weight_func) must correspond
#   Note that in principle we could skip this re-ordering (the only constraint is that x and weight_func must correspond)
#   But we use the re-ordering for the trading-part in mdfa_mse_trade_func: the first column is always the traget that must be traded
  series_ordering<-c(i_series,(1:ncol(log_FX_mat))[-i_series])
  x<-na.exclude(diff(log_FX_mat[,series_ordering]))
# Spectrum: first series is target, all other series are explanatory
# Ordering of series must correspond to x (data for filtering)  
  weight_func<-cbind(weight_func_mat[,i_series],weight_func_mat[,series_ordering])
  colnames(weight_func)<-c(colnames(log_FX_mat)[i_series],colnames(log_FX_mat)[series_ordering])
  K<-nrow(weight_func)-1
  #weight_func[1,]<-1000
  periodicity<-5
  Lag<-0
  L<-2*periodicity#
  L<-min(2*K,L)
  sign_day_of_week<-rep(1,5)
  #  sign_day_of_week[5]<--1
  plot_T<-F
  lag_fx<-1
  plot_T<-T
  
  # This function is in mdfa_mse_trade_func.r
  # MSE estimation and trading: sign rule and signal weighting
  
  mdfa_mse_trade_obj<-mdfa_mse_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)
  
  cum_perf_sign<-mdfa_mse_trade_obj$cum_perf_sign
  yhat<-mdfa_mse_trade_obj$yhat
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


