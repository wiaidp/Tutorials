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
#devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)


#-----------------------------------------------------------------------------------------------
# Source common functions

source("Common functions/plot_func.r")
source("Common functions/mdfa_trade_func.r")
source("Common functions/data_load_functions.R")


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


#----------------------------------------------------------------------------------------------
# Example 5: advanced regularization
#  We here skip the regularization wrapper and use the fulle-fleshed mdfa_analytic function

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

# Set all remaining parameters as in example above
source("Common functions/parameter_set.r")


# The regularization wrapper preselects additionally the following (hyper-) parameter values
lin_eta<-F
weight_constraint<-rep(1/(ncol(weight_func_mat)-1),ncol(weight_func_mat)-1)
lin_expweight<-F
shift_constraint<-rep(0,ncol(weight_func_mat)-1)
grand_mean<-F
b0_H0<-NULL
c_eta<-F
weights_only<-F
weight_structure<-c(0,0)
white_noise<-F
synchronicity<-F
lag_mat<-matrix(rep(0:(L-1),ncol(weight_func_mat)),nrow=L)
troikaner<-F
i1<-i2<-F


# The above selection is 'just fine' for typical economic data (details will be provided in the book)
#  Once (hyper-)parameters are set, the principal estimation function is called in the wrapper

mdfa_obj<-mdfa_analytic(L,lambda,weight_func_mat,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                        lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                        b0_H0,c_eta,weight_structure,white_noise,
                        synchronicity,lag_mat,troikaner)


# This is the same as the reg-wrapper used in the previous examples

mdfa_reg_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj

# Check: both filters are identical
cbind(mdfa_obj$b,mdfa_reg_obj$b)


