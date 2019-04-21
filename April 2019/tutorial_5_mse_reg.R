# Todos
#   -link to data
#   -advanced regularization with mdfa_analytic


# Purpose of tutorial: tackle overfitting by imposing shrinkage of the parameter space towards 'universally' meaningful subspace
#   -Customization will be tackled in separate tutorials
# 1. Show that extendend function replicates previous unconstrained MSE design(s)
# 2. Introduce regularization troika: effect on filter coefficients (universally meaningful shrinkage)
#     -Illustrate decay
#     -Illustrate smoothness
#     -Illustrate cross
# 3. Compare unconstrained and constrained designs
# 4. Advanced regularization



rm(list=ls())

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA",force=T)
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)


#-----------------------------------------------------------------------------------------------
# Source common functions

source("Common functions/plot_func.r")
source("Common functions/mdfa_trade_func.r")
source("Common functions/data_load_functions.R")
source("Common functions/arma_spectrum.r")
source("Common functions/ideal_filter.r")
source("Common functions/tic_simulation_experiment.r")


#-----------------------------------------------------------------------------------------------
# Data: FX 

data_from_IB<-T
hour_of_day<-"16:00"
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

#-----------------------------------------------------------------
# Example 0: illustration of regularization features (regularization troika)

# Set all parameters to 'default' values

in_sample_span<-"2017-01-01"
asset<-"EURUSD"

x<-na.exclude(diff(log_FX_mat[,asset]))
# Spectrum: periodogram
weight_func<-cbind(per(x[paste("/",in_sample_span,sep="")],T)$DFT,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
K<-nrow(weight_func)-1
colnames(weight_func)<-c("spectrum target","spectrum explanatory")
# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity<-10
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-50
# MSE
lambda<-eta<-0
# New parameters for regularization
lambda_smooth<-lambda_cross<-0


# Example 0.1 Decay regularization (single most important type of regularization)
# Main idea: data in the remote past should be weighted less heavily or, stated otherwise, filter coefficients should decay towards zero with increasing lag
# Two parameters: we here loke at first one (shape) 

lambda_decay_1<-0.1*1:9
lambda_decay_2<-0.5


b_mat<-matrix(nrow=L,ncol=(ncol(weight_func)-1)*length(lambda_decay_1))
for (i in 1:length(lambda_decay_1))
{
  lambda_decay<-c(lambda_decay_1[i],lambda_decay_2)
  mdfa_obj_decay<-MDFA_reg(L,weight_func,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj

  b_mat[,(i-1)*(ncol(weight_func)-1)+1:(ncol(weight_func)-1)]<-mdfa_obj_decay$b
}

par(mfrow=c(1,1))
mplot <- b_mat[,1+(0:(length(lambda_decay_1)-1)*(ncol(weight_func)-1))]
ax <- Lag + 0 : (L-1)
colo<-rainbow(ncol(mplot))
insamp<-1.e+90
plot_title <- "Series 1"
title_more<-paste("lambda_decay=(",lambda_decay_1,",",lambda_decay_2,")",sep="")
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)





# Example 0.2 Decay regularization (single most important type of regularization)
# Two parameters: we here look at second one (strength)


lambda_decay_1<-0.5
lambda_decay_2hh<-0.1*0:9
lambda_decay_2h<-c(lambda_decay_2hh,lambda_decay_2hh[length(lambda_decay_2hh)]+0.01*0:9,0.995,0.999)
lambda_decay_2<-lambda_decay_2h[c(1+4*c(0:(length(lambda_decay_2h)/4)),length(lambda_decay_2h)-1,
                                  length(lambda_decay_2h))]

b_mat<-matrix(nrow=L,ncol=(ncol(weight_func)-1)*length(lambda_decay_2))
for (i in 1:length(lambda_decay_2))#i<-2
{
  lambda_decay<-c(lambda_decay_1,lambda_decay_2[i])
  
  mdfa_obj_decay<-MDFA_reg(L,weight_func,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj
  
  b_mat[,(i-1)*(ncol(weight_func)-1)+1:(ncol(weight_func)-1)]<-mdfa_obj_decay$b
}

par(mfrow=c(1,1))
mplot <- b_mat[,1+(0:(length(lambda_decay_2)-1)*(ncol(weight_func)-1))]
ax <- Lag + 0 : (L-1)
colo<-rainbow(ncol(mplot))
insamp<-1.e+90
plot_title <- "Filter Coefficients Series 1"
title_more<-paste("lambda_decay=(",lambda_decay_1,",",lambda_decay_2,")",sep="")
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)


# Comments
# The first parameter controls for the shape: 
#   -if lambda_decay[1]=0 all lags are addressed equally
#   -if lambda_decay[1]>0 then lags in the remote past are affected more heavily 
# The second parameter controls for the strength
#   -if lambda_decay[2]>0 then the coefficients are shrunken towards zero
#     -depending on lambda_decay[1] coefficients in the remote past are shrunken more heavily
#   -In the case of extreme regularization (lambda_decay[2]=1) all coefficients vanish
#     -the degrees of freedom shrink from L to 0



# Example 0.3 Smoothness regularization 
# Main idea: filter coefficients should not be too noisy

lambda_decay<-c(0,0)
lambda_cross<-0

lambda_smooth_vechh<-0.1*0:9
lambda_smooth_vech<-c(lambda_smooth_vechh,lambda_smooth_vechh[length(lambda_smooth_vechh)]+0.01*0:9,0.995,0.999,1)
lambda_smooth_vec<-lambda_smooth_vech[c(1+4*c(0:(length(lambda_smooth_vech)/4)),length(lambda_smooth_vech)-1,
                                        length(lambda_smooth_vech))]


b_mat<-matrix(nrow=L,ncol=(ncol(weight_func)-1)*length(lambda_smooth_vec))
for (i in 1:length(lambda_smooth_vec))
{
  lambda_smooth<-lambda_smooth_vec[i]

  mdfa_obj_decay<-MDFA_reg(L,weight_func,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj
  
  b_mat[,(i-1)*(ncol(weight_func)-1)+1:(ncol(weight_func)-1)]<-mdfa_obj_decay$b
}

par(mfrow=c(1,1))
mplot <- b_mat[,1+(0:(length(lambda_smooth_vec)-1)*(ncol(weight_func)-1))]
ax <- Lag + 0 : (L-1)
colo<-rainbow(ncol(mplot))
insamp<-1.e+90
plot_title <- "Series 1"
title_more<-paste("lambda_smooth=",lambda_smooth_vec,sep="")
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)

# Comments
# For increasing lambda_smooth the coefficients appear 'smoother': reduction of degrees of freedom (less overfitting)
# Smoothness (of the filter coefficients) is generally prefered (though not too smooth...)
# Extreme smoothness means that coefficients are linear: the number of degrees of freedom drops from L to 2
# Smothness does not control for decay of coefficients with increasing lag
# Strong smoothness can conflict with rapid/enforced decay


# Example 0.4 Cross-sectional tightness regularization 
# Main idea: if series in multivariate design are similar, then filter coefficients should be similar too

# We here add an explanatory series, EURJPY, to the design: we are thuis working with a bivariate design
data_cross<-na.exclude(diff(log_FX_mat[,c("EURUSD","EURUSD","EURJPY")]))
# Compute multivariate spectrum
for (i_series in 1:ncol(data_cross))
{
  x<-data_cross[,i_series]
  # Spectrum: white noise assumption
  
  if (i_series==1)
  {
    weight_func_mat<-per(x[paste("/",in_sample_span,sep="")],T)$DFT  
  } else
  {
    weight_func_mat<-cbind(weight_func_mat,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
  }
}
colnames(weight_func_mat)<-colnames(data_cross)

# First column is target; columns 2-3 are explanatory; columns 1 and 2 are identical because EURUSD is also explanatory
head(weight_func_mat)
# Set decay and smooth terms to zero
lambda_decay<-c(0,0)
lambda_smooth<-0
# Select a arnge of cross values in [0,1]
lambda_cross_vec<-c(0,0.05,0.1,0.15,0.2,0.25,0.4,1)
b_mat<-matrix(nrow=L,ncol=(ncol(weight_func_mat)-1)*length(lambda_cross_vec))
# Compute solutions for all cross-values
for (i in 1:length(lambda_cross_vec))#i<-1
{
  lambda_cross<-lambda_cross_vec[i]
  
  mdfa_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj
  
  b_mat[,(i-1)*(ncol(weight_func_mat)-1)+1:(ncol(weight_func_mat)-1)]<-mdfa_obj$b
}

tail(b_mat)

# Plot
par(mfrow=c(1,2))
mplot <- b_mat[,1+(0:(length(lambda_cross_vec)-1)*(ncol(weight_func_mat)-1))]
ax <- Lag + 0 : (L-1)
colo<-rainbow(ncol(mplot))
insamp<-1.e+90
plot_title <- "Series 1"
title_more<-paste("lambda_cross=",lambda_cross_vec,sep="")
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
mplot <- b_mat[,2+(0:(length(lambda_cross_vec)-1)*(ncol(weight_func_mat)-1))]
ax <- Lag + 0 : (L-1)
plot_title <- "Series 2"
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)

# Comments
#   For lambda_cross=0 we see that most weight is attributed to EURUSD (left panel corresponds to EURUSD in the above plot)
#   For increasing lambda_cross both filter coefficients (panels left/right) get closer together
#   For lambda_cross=1 both coefficients are identical: the same filter is applied to both series
#     In that case, MDFA uses information from both data sets to find a common filter
#     Possibly meaningful if both explanatory series correlate positively (otherwise this design would be a clear misspecification)
#     Number of freely determined parameters drops from 2*L to L: less overfitting
#-------------------------------------------------
# More examples: 
#  Trading applications are provided at the end of the file
#  The corresponding function mdfa_reg_trade_func is part of the file mdfa_trade_func.r
#  This function allows more flexibility than mdfa_mse_reg_trade_func used in the previous MSE-tutorial (the latter applies MSE without regularization only)


#-------------------------------------------------
# Example 1
# This example illustrates replication of original (unconstrained) MSE by MDFA_reg (regularization wrapper)
# General setting
# DFA-lowpass: applied to log-returns
# MSE
# Univariate
# Periodogram
#   Regularization controls parameters 


in_sample_span<-"2017-01-01"
asset<-"EURUSD"

x<-na.exclude(diff(log_FX_mat[,asset]))
# Spectrum: periodogram
weight_func<-cbind(per(x[paste("/",in_sample_span,sep="")],T)$DFT,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
K<-nrow(weight_func)-1
#weight_func[1,]<-1000
periodicity<-10
# Cutoff frequency
cutoff<-pi/periodicity
# Target 
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# Real-time estimate
Lag<-0
# Filter length
L<-100
L<-min(2*K,L)
# Regularization: replicates unconstrained MSE 
lambda_cross<-0
lambda_decay<-c(0,0)# try lambda_decay<-c(0,0.1)
lambda_smooth<-0
# MSE-criterion
lambda<-eta<-0

# MDFA_reg: wrapper for working with regularization
mdfa_reg_obj<-MDFA_reg(L,weight_func,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj


names(mdfa_reg_obj)
par(mfrow=c(1,1))
# Filter coefficients: they decay slowly and their appearance is noisy
ts.plot(mdfa_reg_obj$b)
# Amplitude: very noisy (overfitting)
ts.plot(abs(mdfa_reg_obj$trffkt))
# Shift: noisy too
ts.plot(Arg(mdfa_reg_obj$trffkt)/((0:K)*pi/K))


# Same as above but using unconstrained MSE-wrapper

mdfa_obj<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 

names(mdfa_obj)
# Filter coefficients
ts.plot(mdfa_obj$b)
# Amplitude
ts.plot(abs(mdfa_obj$trffkt))
# Shift
ts.plot(Arg(mdfa_obj$trffkt)/((0:K)*pi/K))

# Comparison of filter coefficients: MSE vs. regularization with zero weights
#  Both solutions should be identical i.e. the original MSE-solution is obtained when reg-weights vanish
cbind(mdfa_reg_obj$b,mdfa_obj$b)

#---------------------------------------------------------------------------------------------
# Example 2
# Playing with decay regularization
#   Most important type of regularization
#   Fundamental idea: if target is mid/short-term, then remote past of data is not useful
# Two parameters
#   1. Strength of decay: lambda_decay[2] in [0,1]
#   2. Rate of decay: lambda_decay[1] in [0,1]
# Decay regularization has an effect on scaling: coefficients are atrracted towards zero
#   -If turning-points are of interest, then this effect could be ignored (level-effect, zero-crossings are not affected)
#   -Otherwise: scale back to original levels, see below
#   -Or use constraints (i1<-T)

# Set all remaining parameters as in example above
source("Common functions/parameter_set.r")

# Example 2.0: no decay
lambda_decay<-c(0,0)
# Example 2.1
# Strong but slow decay
lambda_decay<-c(0.1,0.9)
# Example 2.2
# Strong and fast decay
lambda_decay<-c(0.7,0.9)
# Example 2.3
# very strong and fast decay
lambda_decay<-c(0.7,0.999)
# Example 2.4
# very strong and slow decay
lambda_decay<-c(0.01,0.999)

# MDFA_reg: wrapper for working with regularization
mdfa_reg_obj<-MDFA_reg(L,weight_func,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj


# Filter coefficients
ts.plot(mdfa_reg_obj$b)
# Amplitude
ts.plot(abs(mdfa_reg_obj$trffkt))
# Shift
ts.plot(Arg(mdfa_reg_obj$trffkt)/((0:K)*pi/K))

# Scale back
scaler<-sum(mdfa_reg_obj$b)
# Filter coefficients
ts.plot(mdfa_reg_obj$b/scaler)
# Amplitude
ts.plot(abs(mdfa_reg_obj$trffkt/scaler))
# Shift
ts.plot(Arg(mdfa_reg_obj$trffkt/scaler)/((0:K)*pi/K))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Example 3
# Playing with smoothness regularization: lambda_smoothness in [0,1]
#   Fundamental idea: filter coefficients should be smoothly changing over time (exception: seasonality)

# Set all remaining parameters as in example above
source("Common functions/parameter_set.r")

# Set lambda_decay back to zero
lambda_decay<-c(0,0)

# Example 3.0
# No smoothness
lambda_smooth<-0
# Example 3.1
# Mild smoothness
lambda_smooth<-0.5
# Strong smoothness
lambda_smooth<-0.99
# Extremely strong (maximum) smoothness: coefficients are 'linear'
lambda_smooth<-1


# MDFA_reg: wrapper for working with regularization
mdfa_reg_obj<-MDFA_reg(L,weight_func,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj


# Filter coefficients
ts.plot(mdfa_reg_obj$b)
# Amplitude
ts.plot(abs(mdfa_reg_obj$trffkt))
# Shift
ts.plot(Arg(mdfa_reg_obj$trffkt)/((0:K)*pi/K))



#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Example 4
# Cross-regularization
#   Fundamenatl idea: if sries in a multivariate design are 'similar', then 'similar' filters should be applied to them 
# Effect of lambda_cross: 
#   -lambda_cross imposes similar coefficients across series
#   -Need multivariate design: use all FX-series with EUR as explanatory data: target is EURUSD
data<-na.exclude(diff(log_FX_mat[,c("EURUSD","EURUSD","EURJPY","EURGBP")]))
# Compute multivariate spectrum
for (i_series in 1:ncol(data))
{
  x<-data[,i_series]
  # Spectrum: white noise assumption
  
  if (i_series==1)
  {
    weight_func_mat<-per(x[paste("/",in_sample_span,sep="")],T)$DFT  
  } else
  {
    weight_func_mat<-cbind(weight_func_mat,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
  }
}

# First column is target; columns 2-4 are explanatory; columns 1 and 2 are identical because EURUSD is also explanatory
head(weight_func_mat)


# Set all remaining parameters as in example above
source("Common functions/parameter_set.r")

# Set lambda_decay and lambda_smooth back to zero
lambda_decay<-c(0,0)
lambda_smooth<-0

# Example 4.0
# No cross-sectional constraint
lambda_cross<-0
# Example 4.1
# Mild cross-sectional constraint
lambda_cross<-0.5
# Example 4.2
# Strong cross-sectional constraint
lambda_cross<-0.9
# Example 4.3
# Very strong cross-sectional constraint
lambda_cross<-1


# MDFA_reg: wrapper for working with regularization
mdfa_reg_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj

par(mfrow=c(1,1))
# Filter coefficients
ts.plot(mdfa_reg_obj$b,col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])
# Amplitude
ts.plot(abs(mdfa_reg_obj$trffkt),col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])
# Shift
ts.plot(Arg(mdfa_reg_obj$trffkt)/((0:K)*pi/K),col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])

# Comments
#   -Imposing cross-sectional regularization is useful when data is positively cross-correlated
#     Counterexample: EURUSD and USDJPY are likely to be negatively cross-correlated
#     In the above example all pairs start in EUR: therefore assumption of positive cross-correlation is OK
#   -Imposing full constraint (all coefficients are the same) mixes information across all series
#     Degrees of freedom correspond to length L of filter


#------------------------------------------------------------------------------------
# Example 5: we here duplicate example 3 (simulation experiment) of tutorial 3 but we use regularization (instead of unconstrained design)
# Specifically
#   -We benchmark empirical DFA based on DFT, unconstrained (as in tutorial 3) and regularized,
#     against best possible MSE-design (assuming knowledge of true model) out-of-sample
# Expectation: in short samples (sparse information) performances of suitably regularized designs improve out-of-sample 

# Example 5.1: MA(1)
a1<-0.
b1<-0.7
true_model_order<-c(0,0,1)
# Example 5.2: ARMA with positive acf
a1<-0.6
b1<-0.7
true_model_order<-c(1,0,1)
# Example 5.3: AR with negative acf
a1<--0.9
b1<-0
true_model_order<-c(1,0,0)
# Example 5.4: close to noise (typical for log-returns of FX-data)
a1<--0.08
b1<-0.0
true_model_order<-c(1,0,0)
# Add any other processes...



# We generate anzsim realizations of length len of the arma-process
set.seed(1)
len<-1000
mse_true_arma<-mse_dfa<-NULL
# Number of simulations
anzsim<-500
# Length of in-sample span
#   For illustration we here use a short sample
#   In fact regularization is deemed useful when working with small samples (otherwise unconstrained designs work fine)
in_sample<-100
# Frequency grid for DFA based on true model
K_true<-600
# Lowpass target  
periodicity<-10
# Deliberately exagerate overfitting: we use filters of length 4*periodicity 
#   -This filter length is used for all designs
#   -Expected effect: when imposing regularization in short samples the performances improve out-of-sample
L<-4*periodicity
#   MSE design (no customization)
lambda<-eta<-0
# Use the most important (decay) regularization: mid-strength
#   See example 7 below for a formal derivation of these regularization settings
lambda_decay<-c(0.6,0.8)
# Smoothness is most useful when estimating long-term trends (large periodicity)
#   But it's usefulness might apply 'in general' too
lambda_smooth<-0.99
# Cross sectional regularization is not activated here since design is univariate  
lambda_cross<-0 
# Nowcast
Lag<-0
# Length of ideal filter
M<-100
mse_true<-mse_dft<-mse_reg_dft<-NULL
pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
# Loop through all simulations and collect out-of-sample forecast performances
for (i in 1:anzsim)#i<-1
{
  # Distinguish white noise  
  if (abs(a1)+abs(b1)>0)
  {
    # Generate series  
    x<-as.vector(arima.sim(n=len,list(ar=a1,ma=b1)))
  } else
  {
    x<-rnorm(len)
  }
  # Use in-sample span for model-estimation and for dft  
  x_insample<-x[1:in_sample]
  # True model: estimate model-parameters by relying on classic arima-function
  arima_true_obj<-arima(x_insample,order=true_model_order,include.mean=F)
  # Spectrum based on true model  
  spec<-arma_spectrum_func(ifelse(!is.na(arima_true_obj$coef["ar1"]),arima_true_obj$coef["ar1"],0),ifelse(!is.na(arima_true_obj$coef["ma1"]),arima_true_obj$coef["ma1"],0),K_true,F)$arma_spec
  weight_func<-cbind(spec,spec)
  colnames(weight_func)<-c("spectrum target","spectrum explanatory")
  weight_func_true<-weight_func
  cutoff<-pi/periodicity
  # target true model: frequnecy grid is not the same as for dft below i.e. K is different
  Gamma_true<-(0:(K_true))<=K_true*cutoff/pi+1.e-9
  mdfa_true_obj<-MDFA_mse(L,weight_func_true,Lag,Gamma_true)$mdfa_obj 
  b_true<-mdfa_true_obj$b
  # Use in-sample span for dft  
  weight_func_dft<-cbind(per(x_insample,F)$DFT,per(x_insample,F)$DFT)
  colnames(weight_func_dft)<-c("spectrum target","spectrum explanatory")
  K_dft<-nrow(weight_func_dft)-1
  Gamma_dft<-(0:(K_dft))<=K_dft*cutoff/pi+1.e-9
  # Compute unconstrained MSE-filter (same as in example 3 of tutorial 3)
  mdfa_dft_obj<-MDFA_mse(L,weight_func_dft,Lag,Gamma_dft)$mdfa_obj 
  b_dft<-mdfa_dft_obj$b
# New: use regularized design
  mdfa_reg_obj<-MDFA_reg(L,weight_func_dft,Lag,Gamma_dft,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj
  b_dft_reg<-mdfa_reg_obj$b
  
  # Filter data
  # 1. ideal filter  
  id_obj<-ideal_filter_func(periodicity,M,x)
  output_ideal<-id_obj$y
  # 2. DFA true
  filt_true_obj<-filt_func(x,b_true)
  output_dfa_true<-filt_true_obj$yhat
  # 3. DFA dft
  filt_dft_obj<-filt_func(x,b_dft)
  output_dfa_dft<-filt_dft_obj$yhat
  # 4. DFA dft regularized
  filt_dft_reg_obj<-filt_func(x,b_dft_reg)
  output_dfa_dft_reg<-filt_dft_reg_obj$yhat
  # Mean-square out-of-sample filter error  
  mse_true<-c(mse_true,mean((output_ideal-output_dfa_true)[(in_sample+1):len-M]^2,na.rm=T))
  mse_dft<-c(mse_dft,mean((output_ideal-output_dfa_dft)[(in_sample+1):len-M]^2,na.rm=T))
  mse_reg_dft<-c(mse_reg_dft,mean((output_ideal-output_dfa_dft_reg)[(in_sample+1):len-M]^2,na.rm=T))
  setTxtProgressBar(pb, i)
}

# Compute the ratio of root mean-square forecast errors:
#   The ratio cannot be larger than 1 asymptotically because our particular design distinguishes arma as the universally best possible design
result_mat<-matrix(c(sqrt(mean(mse_true)/mean(mse_dft)),sqrt(mean(mse_true)/mean(mse_reg_dft))),nrow=1,ncol=2)
colnames(result_mat)<-c("Unconstrained","Regularized")
result_mat
# Results: 
#   -Overfitting could be sucessively avoided by the regularized design
#   -Efficiency lies at 93% when compared to best possible design, assuming knowledge of true process (quite remarkable given the short estimation span)

# Note that we stuffed the above code-lines into a more convenient function-call
#   This function replicates the above simulation experiment
#   It computes additional statistics which are described in example 7 below
tic_simulation_experiment_func(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)
  

#----------------------------------------------------------------------------------------------
# Example 6: advanced regularization
# -We here exploit additional (computationally slightly more expensive) statistics which hint at 
#   'optimal regularization' settings
# -What is an 'optimal regularization' setting?
#   1. The idealized purpose of regularization is to contain overfitting without affecting (too much...) flexibility
#     -Ideally, an optimally regularized design is able to match all relevant features of a particular estimation problem 
#     -Ideally, an optimally regularized design performs as well out-of-sample as it does in-sample
#   2. More practically (less idealized) we expect that hinting at an 'optimally regularized design' will affect 
#       positively out-of-sample performances
#     -Ideally, the 'optimally regularized design' outperforms all other contenders out-of-sample
# -We here propose a novel statistic (see book for formal background) which (in a way) generalizes classic information
#   -We derive the effective degrees of freedom of a troika-regularized design (decay, smoothness and cross sectional similarity)
#   -From this we can derive the effective dimension of the shrinkage space
#   -From this we can derive the spurious decrease of the DFA-(optimization) criterion in the case of overparametrization (too richly parametrized filters)
#   -Form this we can design an 'information criterion' by compensating the spurious decrease (of the criterion) by 
#     a 'penalty'-term which effectively penalizes too-richly parametrized designs

# Let's do so!
#   Example 6 plays with the new edof statistic (effective degrees of freedom)
#   Example 7 plays with the new tic statistic (troika information criterion) which relies on edof
#   Both statistics are 'computationally expensive': therefore we can turn them off via a Boolean called troikaner
#     If trokaner==F then the ststistics are skipped and the calculation (optimization) is faster
#     If troikaner==T then the statistics are computed; the result of the optimization is not affected though


in_sample_span<-"2017-01-01"
# Compute spectrum
data<-na.exclude(diff(log_FX_mat[,c("EURUSD","EURUSD","EURJPY","EURGBP")]))
# Compute multivariate spectrum
for (i_series in 1:ncol(data))
{
  x<-data[,i_series]
  # Spectrum: white noise assumption
  
  if (i_series==1)
  {
    weight_func_mat<-per(x[paste("/",in_sample_span,sep="")],T)$DFT  
  } else
  {
    weight_func_mat<-cbind(weight_func_mat,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
  }
}


# First column is target; columns 2-4 are explanatory; columns 1 and 2 are identical because EURUSD is also explanatory
head(weight_func_mat)



# Example 6.1: replication of wrapper by generic estimation function
L<-20
# The regularization wrapper preselects additionally the following (hyper-) parameter values
K<-nrow(weight_func_mat)-1
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# nowcast
Lag<-0
# Periodicity
periodicity<-10
cutoff<-pi/periodicity
# MSE
lambda<-eta<-0
# No regularization
lambda_smooth<-lambda_cross<-0
lambda_decay<-c(0,0)

# Classic call (as above)
mdfa_reg_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj
# New call: with additional Boolean troikaner
#   If T then additional statistics will be computed (however computation time increases) but filter solution is not affected
#   Troika here means: the three supports of regularization, namely decay, smoothness, cross-similarity
troikaner<-T
mdfa_reg_T_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth,troikaner)$mdfa_obj


# Check: both filters are identical
cbind(mdfa_reg_T_obj$b,mdfa_reg_obj$b)

# But there is more output (statistics) available when troikaner==T
names(mdfa_reg_obj)
names(mdfa_reg_T_obj)


#------------------------
# Example 6.2 working with the degrees of freedom


# We estimated L*3 (trivariate filter) coefficients without regularization in the above example i.e. 20*3=60 
mdfa_reg_T_obj$edof

#-----------------------------
# Example 6.3 Imposing super strong decay
#   We here observe how this affects the degrees of freedom

# Super strong decay
lambda_decay<-c(0,1)
lambda_smooth<-0.
lambda_cross<-0.

troikaner<-T
mdfa_reg_T_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth,troikaner)$mdfa_obj

# We estimate L*3 coefficients and impose super-strong decay 
#   The parametgers vanish and the degrees of freedom are 0
#   Meaningful though not useful...
round(mdfa_reg_T_obj$edof,0)

#---------------------------
# Example 6.4 Imposing super strong smoothness
# We now impose superstrong smoothness: all other reg-terms vanish
#   In this case all coefficients are linear
#   A linear function is determined by 2 degrees of freedom
#   Trivariate design times two degrees of freedom implies a total of 6 degrees of freedom
# Let's check
lambda_decay<-c(0,0)
# Super strong smoothness
lambda_smooth<-1
lambda_cross<-0

troikaner<-T
mdfa_reg_T_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth,troikaner)$mdfa_obj

# We estimate L*3 coefficients and impose super-strong smoothness 
#   super strong smoothness means: coefficients are linear
#   a linear function requires two degrees of freedom (slope/intercept)
#   therefore we expect the dgrees of freedom to be 3*2=6
round(mdfa_reg_T_obj$edof,0)


#---------------------------
# Example 6.5 Imposing super strong cross-sectional
#  The coefficients must be identical across series
#  Therefore the degrees of freedom must by 1*L 

# Let's check
lambda_decay<-c(0,0)
lambda_smooth<-0
lambda_cross<-1

troikaner<-T
mdfa_reg_T_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth,troikaner)$mdfa_obj

# We estimate L*3 coefficients and impose super-strong cross sectional regularization 
#   Coefficients accross series must be identical (but are otherwise unconstrained)
#   Therefore we xpect the degrees of freedom to be 1*L (instead of 3*L)
round(mdfa_reg_T_obj$edof,0)

#-----------------------------
# Example 6.6 Imposing arbitrary regularization
#   We here observe how this affects the degrees of freedom

# Strong and fast decay
lambda_decay<-c(0.7,0.9)
lambda_smooth<-0.9
lambda_cross<-0.9

troikaner<-T
mdfa_reg_T_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth,troikaner)$mdfa_obj

# We estimate L*3 coefficients no regularization: in the above example L=20 i.e. degrees of freedom is 3*20=60 (without regularization)
#   The degrees of freedom are substantially smaller: around 5
round(mdfa_reg_T_obj$edof,0)

#-----------------------------------------------------------------------------------------------------
# Example 7: working with the generalized DFA information criterion (tic)
# Question: how do we have to choose the regularization parameters (lambda_decay, lambda_smooth, lambda_cross)?
#   -First possibility: cross-validation (determine settings such that out-of-sample MSE-performances are optimized, see example 5 above)
#   -Besides this purely empirical approach there is a formal alternative based on the dimension of the shrinkage-space (the degrees of freedom used above)
# Select lambda_decay, lambda_smooth, lambda_cross so as to minimize the DFA information criterion tic
#   -This criterion generalizes the classic Akaike criterion AIC to the more complex regularized signal extraction problem 
#   -It is not based on the idea of 'maximizing a likelihood' (as is AIC)
#   -Instead it derives the dimension of the regularized shrinkage space from which the spurious decrease of the 
#     (log of the) MSE-criterion can be obtained in the case of overparametrization (L too large)
#   -Compensating the spurious decrease of the (log of the) MSE-criterion by a penalty-term larger than this decrease then results in the new criterion

# Ideally we expect that
#   -Selecting lambda_decay, lambda_smooth, lambda_cross according to minimum tic should lead to best out-of-sample performances
# Let's check that idealized expectation
#   -The following code is the same as in example 5 above excpet that we set troikaner=T when calling the regularization wrapper in the loop

# Number of simulations
anzsim<-100

# Example 7.1: MA(1)
a1<-0.
b1<-0.7
true_model_order<-c(0,0,1)

# Example 7.1.1
# This regularization setting is a good choice in the case of 'nearly white noise'
#   -Degrees of freedom is around 10
#   -The tic value is small and out-of-sample performances (of the regularized design) are good (92% efficiency vs. 80% for unconstrained designs) 
lambda_smooth<-0.
lambda_decay<-c(0.6,0.8)
tic_simulation_experiment_func(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)

# Example 7.1.2
# Weaker lambda_decay
#   -Degrees of freedom is a bit larger: 13
#   -tic is sligthly larger than first setting and out-of-sample MSE is slightly worse (as desired)
lambda_decay<-c(0.3,0.8)
tic_simulation_experiment_func(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)


# Example 7.1.3
# Substantially weaker lambda_decay
#   -Degrees of freedom is larger: 23
#   -tic is sligthly larger than first setting and out-of-sample MSE is worse (as desired)
#   -Overfitting is starting to compromise out-of-sample performances
lambda_decay<-c(0.3,0.2)
tic_simulation_experiment_func(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)


# Example 7.1.4
# Stronger lambda_decay
#   -Degrees of freedom is smaller (than first design): 7
#   -tic is sligthly larger than first setting and out-of-sample MSE is worse (as desired)
#   -in this case our design is 'too simple': it cannot fit the (complex) structure of the estimation problem as well as first design
lambda_decay<-c(0.6,0.95)
tic_simulation_experiment_func(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)


# Example 7.1.5
# We add smoothness regularization to first design
lambda_smooth<-0.3
#   -Degrees of freedom is around 7
#   -The tic value is smaller than first design but out-of-sample performances are marginally worse (might be attributed to sampling error because we know formally that tic makes sense) 
lambda_decay<-c(0.6,0.8)
tic_simulation_experiment_func(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)


# Example 7.1.6
# Strong smoothness 
lambda_smooth<-0.99
#   -Degrees of freedom is around 3
#   -The tic value is smaller than first design and out-of-sample performances are marginally better 
#   -in this case the additional regularization shrinkage allows to fit the data well, without overfitting
#     tic suggests that 3 degrees of freedom (optimally shared across the 40 coefficients) is all we need to match the salient features of the estimation problem
#     out-of-sample MSE performances confirm that claim
lambda_decay<-c(0.6,0.8)
tic_simulation_experiment_func(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)

# Example 7.1.7
# Too much smoothing
lambda_smooth<-0.9999
#   -Degrees of freedom is 1 (now that's small...)
#   -The tic value is substantially larger than the previous design and out-of-sample performances are worse 
#   -Shrinkage is too strong...
lambda_decay<-c(0.6,0.8)
tic_simulation_experiment_func(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)


# Example 7.1.8
# Unconstrained design
lambda_smooth<-0.
#   -Degrees of freedom is 40 (which corresponds to L)
#   -The tic value is substantially larger than the above 'good' designs and out-of-sample performances are substantially worse 
#   -overparametrization (too many degrees of freedom) generates and exacerbates overfitting
lambda_decay<-c(0.,0.)
tic_simulation_experiment_func(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)

# Conclusion
#   -Best out-of-sample performances are obtained by designs which strongly shrink the original (unconstrained) estimation space
#   -the new information criterion --tic-- and out-of-sample MSE-performances correlate strongly (efficiency correlates negatively i.e. smaller tic implies higher efficiency)
#   -Designs with degrees of freedom between 13 (example 7.1.2), 10 (example 7.1.1) or 3 (example 7.1.6) perform quite similarly out-of-sample
#     -Therefore there is no urge of determining the uniquely absolute best possible design since MSE-performances are not 
#       too 'touchy' with respect to departures from the optimum
#     -Strategy: find a regularization setting with a fairly small tic (which is automatically greeted with fairly small degrees of freedom)
#       imposing an appealing mix of decay, smoothness and cross-sectional (the latter only for multivariate designs) regularization
#     -Forget the delusional concept of 'the single best outcome'
#   -Unnecessarily strong regularization may lead to unflexible designs with poorer out-of-sample performances

#------------------
# Example 7.2: Same as above examples 7.1.1-7.1.8 but use a 'nearly white' process (typical for log-returns of financial data)
#   The big picture is the same (nearly identical results as in example 1 up to sampling errors)
a1<--0.08
b1<-0.0
true_model_order<-c(1,0,0)


#-----------------------------------------------------------------------------------------------------------
# What did we learn
#   -Overfitting can seriously affect out-of-sample MSE-performances
#   -Multivariate designs suffer more heavily from overfitting
#   -The regularization troika based on decay, smoothness and similarity, imposes universally desirable filter-features
#     -Irrespective of the particular application, all three requirements are 'desirable'
#     -Therefore we may impose any combination of them in any kind of application (without fearing some 'loss of structural integrity' of the resulting filter)
#   -However, we don't know a priori how to weight the trems of the troika: 
#     -the strength (dimension of shrinkage space) is dependent on the structure of the particular estimation problem at hand
#   -We proposed a novel statistic, the troika information criterion (tic) which enables to find suitable regularization
#     weighting 
#     -Suitable here means: the tic correlates strongly (though negatively) with out-of-sample MSE-performances 
#   -In general the effective degrees of freedom (edof) of such optimized designs are markedly smaller than the original unconstrained designs
#     -The Troika is flexible enough (universal) to account for the relevant features of an estimation problem
#       with humbly-sized degrees of freedom


# Other/additional issues
#   1.The proposed regularization features are effective in terms of MSE- (level) performances, mainly
#     -Alternative applications for which MSE (or level) is not of utmost priority are not adressed explicitly 
#     -Example: trading or recession-calling
#       -In trading the analyst may decide to be positioned (long/short) depending on the sign of the filter output (positive/negative)
#       -Considering the filter output in terms of 'sign' (instead of level) is a scale-invariant application
#       -On the other hand, zero-shrinkage (by imposing a large decay-regularization for example) adresses 
#         MSE/level issues, explicitly (shrinking towards zero improves MSE under noisy conditions)
#       -Therefore zero-shrinkage cannot be claimed to address explicitly scale-invariant filter applications  
#     -Though indirectly most applications will still be addressed, at least indirectly/partially, because 
#         the troika-features may claim 'universality' (to some extent) 
#       -Filters with smooth and fast-decaying coefficients are felt to be useful irrespective of the 
#         particular application at hand (level/scale-dependent or scale-invariant)

#   2.Regularization and spectrum estimates
#     -The degrees of freedom statistic computed by the MDFA-package (edof) measures the dimension of the 
#       projection-space engrossed by the MDFA-estimate
#     -Shrinking that space to a smaller structurally commensurate ('smart') sub-space by imposing universally valid 
#       features (smoothness, decay, cross-similarity) addresses overfitting and out-of-sample performances (while still maintaining sufficient flexibility to tackle complex/generic estimation problems)
#     -However, overfitting is mainly (in fact: uniquely) due to our choice of the weighting-function (the dft) being 'noisy' 
#       -Overfitting could be avoided entirely by substituting the true spectrum (if known...) to the dft, as we did in tutorial 3 (nowcasting univariate MSE)
#         -The true spectrum has zero degrees of freedom: we do not fit the data (we just use the true model)
#       -In contrast, the dft has as many degrees of freedom as the original data: we can recover the original data from the dft by inverse fourier transform (idft)
#     -Therefore, MDFA inherits the degrees of freedom from its spectrum estimate: 
#       -if the spectrum is the dft then 
#         -the degrees of freedom of the dft corresponds to the number of observations 
#         -the degrees of freedom of the (unconstrained) MDFA corresponds to L (univariate) or L times #explanatory series (multivariate)
#         -overfitting (poor out-of-sample performances) by the (unconstrained) MDFA will be unavoidable, as L increases. 
#       -if the spectrum is a model-based spectrum (for example a parsimonious ARMA-spectrum, see tutorials 2 and 3) then 
#         -the degrees of freedom of MDFA cannot exceed the degrees of freedom of the fitted model
#           -An ARMA(1,1)-model would result in maximally 2 degrees of freedom for the MDFA, irrespective of L, 
#             as long as L/K is small 
#             -If L/K is 'large', then overfitting of the target at the discrete frequency-ordinates could result, recall tutorial 2
#             -Recall that for a model-based estimate K could be selected arbitrarily large (a large K 'consumes' more computation time, though)
#         -therefore, in this case (assuming the spectrum is based on a parsimonious model) the effectiveness of the proposed regularization troika 
#           would be mitigated: the troika is not necessarily/strictly required in such a setting
#       -if the spectrum is obtained from Burg's maximum entropy estimate then
#         -The degrees of freedom of MDFA cannot exceed the order of the autoregressive model fitted to the data
#           -Assuming, once again, that L/K is small...
#         -For large autoregressive model orders, the troika can be applied to contain (unavoidable) overfitting and to improve out-of-sample performances
#     -Now... given that overfitting might be contained/addressed/tackled by model-based spectra: 
#       -why don't we use more often (at all) model-based spectra in our applications?
#       -why do we use instead the noisy dft together with a fancy dimension-crusher (regularization troika)
#     -Answers
#       -In practice, models are invariably (always) misspecified; the dft never is (the dft IS the data: trivially, it cannot be misspecified)
#         -By using a model instead of the dft we trade bias (misspecification) against volatility (noise)
#       -While models for univariate designs perform well (as confirmed by a historical record of forecast competitions), models
#         for multivariate problems are tricky (there are numerous identification/causality issues and the number of parameters to be fitted increaes rapidly i.e. overfitting is an issue)
#       -Controlling the degres of freedom of the MDFA (based on dft) by our fancy dimension-crusher (the troika) addresses
#         bias (universality of the troika-requirements implies mitigation of misspecification) as well as volatility 
#         in a single-stroke, within an estimation paradigm that adresses the structure of the relevant problem (signal extraction/nowcasting) directly (in contrast to model-based approaches which address one-step forecasting, exclusively)
#     -Note that the troika could be applied to MDFA as based on arbitrary spectra (for example model-based spectra with few degrees of freedom)
#       -But applying 'universal' constraints to a design which is inherently misspecified misses somehow the point...

#   Wrap-up
#     -The extent of overfitting of the data by MDFA cannot exceed the extent of overfitting of the data-generating process by the spectral estimate
#       If the spectrum is not overfitted then the MDFA cannot overfit the data (assuming L/K is small): in-sample and out-of-sample performances will be congruent
#       However, a strongly misspecified spectrum (bias) will result in poor (in-sample as well as out-of-sample) performances
#       The dft in combination with the troika aims at good in-sample as well as congruent (and therefore equally good) out-of-sample performances
#     -The regularization troika is required (most effective) when using the (very noisy but 'not misspecified') dft as spectrum
#     -The combination of troika and dft mitigates bias (misspecification) and volatility (dimension-crusher) in a coherent 'single-stroke' approach 