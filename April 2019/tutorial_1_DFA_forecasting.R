# This tutorial is made of 6 exercises

# Purpose of tutorial: 
# DFA can replicate classic (MSE) one- and multi-step ahead forecasting by 
#   -specifying a corresponding target (allpass) and forecast horizon (Lag) 
#   -specifying a corresponding spectral estimate
# We illustrate, in particular, how classic ARIMA-based forecasting can be replicated in the DFA framework. 
#   -Once rooted into well-known territory, we can deploy the additional flexibility of DFA in the following tutorials 
# At the end we also illustrate pertinence of a non-parametric DFA approach based on the discrete fourier transform (dft) 
#     as well as Burg's max-entropy spectral estimate
#   -We show that 
#     1. Out-of-sample forecast performances of DFA based on dft are nearly as good as the universally best approach (which assumes knowledge of the true data-generating process)
#     2. Performances of DFA based on dft and of DFA based on max-entropy spectrum are indistinguishable  
#   -These results suggest pertinence of the non-parametric approach --DFA based on dft-- in a wide range of applications 
#     and in particular when models are likely to be misspecified (which is always the case for real-world data)

rm(list=ls())

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)


# Brief overview of wrappers and main MDFA function
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
source("Common functions/arma_spectrum.r")


#-------------------------------------------------------------------------------------------------
# Example 1: one-step ahead forecasting
#   DFA requires the user to supply a 'target' and a 'spectrum'; DFA then returns an optimal estimate 
#     -Optimality: Mean-square error (MSE) or beyond MSE (ATS-trilemma)
# In this example: we illustrate how to set-up one-step ahead forecasting in DFA
# General design decisions:
#   MSE (mean-square error criterion)
#   Univariate design
#   White noise spectrum

# DFA is specified in the frequency-domain, see McElroy/Wildi
#   -K defines the frequency-grid for estimation: all frequencies omega_j=j*pi/K, j=0,1,2,...,K are considered
#   -Tradeoff: larger K (denser grid) means better MSE (if true process is white noise) but longer computation-time
#   -We here consider a partition of the interval [0,pi] in K=600 equally-distanced frequency supports
K<-600
# Spectrum: white noise assumption
#   -First column is spectrum of target, second column is spectrum of explanatory variable
#   -In a univariate design target and explanatory data are the same
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
colnames(weight_func)<-c("spectrum target","spectrum explanatory")

# White noise: flat spectrum (all frequencies are loaded equally by the process)
plot(weight_func[,1],type="l",main=paste("White noise spectrum, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext(colnames(weight_func)[1],line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Target forecasting: 
#   -When forecastin a series we are interested in all frequencies equally
#   -Target (Gamma) is a forward-looking allpass filter in the frequency-domain
Gamma<-rep(1,K+1)
plot(Gamma,type="l",main=paste("Allpass forecast target, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()
# One step ahead: Lag=-1
Lag<--1
# Filter length: number of weights/coefficients of forecast filter
L<-10

# Estimation based on MDFA-MSE wrapper
mdfa_obj<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 

# We first plot the filter coefficients: 
#   -The true coefficients are zero (white noise spectrum)
#   -Our MSE-solution is close to (but not identical to) zero
#     -This is because the filter-grid is of finite resolution (K=600 above)
#     -Selecting a larger K will shrink estimates towards zero (but computation time increases accordingly)
#   -For any finite K we thus observe some overfitting
#     -As illustrated below the problem is negligible in the context of this tutorial
par(mfrow=c(1,1))
b<-mdfa_obj$b
colo<-rainbow(ncol(b))
plot(b[,1],type="l",main=paste("Filter coefficients",sep=""),
     axes=F,xlab="Lag",ylab="Coef",ylim=c(min(b),max(b)),col="black")
mtext(colnames(weight_func)[2],line=-1,col="black")
if (ncol(b)>1)
{
  for (i in 2:ncol(b))
  {
    lines(b[,i],col=colo[i])
    # We take the i+1 colname from weight_func because the first column is the target        
    mtext(colnames(weight_func)[i+1],line=-i,col=colo[i])
    
  }
}
axis(1,at=1:L,labels=1:L)
axis(2)
box()    

# The following function plots 
#   -the coefficients
#   -the time-shift and
#   -the amplitude of the filter
plot_estimate_func(mdfa_obj,weight_func,Gamma)

# Comments
#   -The amplitude function of the MSE-filter is close to zero, as desired (ideally, it would equate to zero)
#   -It is not exactly zero because K is finite
#   -For K/L>10 mismatch is generally negligible  


#-----------------------------------
# Example 2: classic one- and multi-step ahead forecasting ARMA-process (replicates SOTA ARMA-models)

# Frequency grid
K<-600
# Target forecasting: 
#   -We are interested in all frequencies equally
#   -Target Gamma is a forward-looking allpass filter in the frequency-domain
Gamma<-rep(1,K+1)

# Spectrum ARMA: try various processes
a1<-0.9
b1<-NULL
a1<-NULL
b1<-0.6

plot_T<-T

# This function computes the spectrum of the ARMA-process
spec<-arma_spectrum_func(a1,b1,K,plot_T)$arma_spec

# Fill into weight_func: target (first column) and explanatory (second column); both are identical for univariate problems
weight_func<-cbind(spec,spec)
colnames(weight_func)<-c("spectrum target","spectrum explanatory")

# Specify Lag
#   1. k-step ahead forecasting: Lag<--k (negative integer or negative real number: in the latter case one forecasts between two consecutive future time points)
#   2. Nowcast: Lag<-0
#   3. Backcast: Lag<-k (positive integer or positive real: in the latter case one interpolates between two consecutive past time points)
Lag<--1
# Filter length: number of weights/coefficients of forecast filter
L<-10

# Estimation based on MDFA-MSE wrapper
mdfa_obj<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 

names(mdfa_obj)
mdfa_obj$b

# The following function plots 
#   -the coefficients
#   -the time-shift and
#   -the amplitude of the filter
plot_estimate_func(mdfa_obj,weight_func,Gamma)



# Comments
# 1. Forecasting: Lag<0 integer-valued
#   -For an AR(1)-process and Lag=-1 (one-step ahead) the optimal forecast filter has weights (a1,0,0,...,0)
#   -For an AR(1)-process and Lag=-2 (two-steps ahead) the optimal forecast filter has weights (a1^2,0,0,...,0)
#   -For an AR(1)-process and Lag=-k (k-steps ahead) the optimal forecast filter has weights (a1^k,0,0,...,0)
#   -Our solution gets arbitrarily close for increasing K (denser grid)
#   -For K/L>10, mismatch is negligible  

#   -For an MA(1)-process and Lag=-1 (one-step ahead) the optimal forecast filter has weights (b1,-(b1^2),b1^3,-(b1^4),...)
#   -For an MA(1)-process and Lag=-2 (two-steps ahead) the optimal forecast filter has weights 0
#   -Our solution gets aribtrarily close for increasing K (denser grid)
#   -For K/L>10 mismatch is negligible 

# 2. Nowcasting: Lag=0
#   -For any ARMA-process the optimal filter has weights 1,0,0,....
#   -This is a trivial estimation problem because our target is an allpass: for lowpass/bandpass/highpass the estimation problem is nomore trivial

# 3. Backcasting: Lag>0
#   -For any process the optimal filter has weight 1 at 'lag' specified by Lag and otherwise 0
#   -This is a trivial estimation problem because our target is an allpass: for lowpass/bandpass/highpass the estimation problem is nomore trivial

# 4. For real-valued Lag the filter interpolates the data between consecutive (future or past) time-points



#-----------------------------------
# Example 3: compare DFA with classic arima R-code

# Specify forecast horizon

h<-5
len<-300
# Data: simulate ARMA (select any coefficients)

a1<-0.6
b1<-0.7
set.seed(0)
x<-arima.sim(n=len,list(ar=a1,ma=b1))

# Estimate model-parameters by relying on classic arima-function
arima_obj<-arima(x,order=c(1,0,1),include.mean=F)

# Diagnostics: model is OK
tsdiag(arima_obj)

# Compute forecasts
arima_pred<-predict(arima_obj,n.ahead=h)$pred

K<-600
# Target forecasting: 
#   -We are interested in all frequencies equally
#   -Target Gamma is a forward-looking allpass filter in the frequency-domain
Gamma<-rep(1,K+1)

# Spectrum: supply estimated filter coefficients
a1<-arima_obj$coef["ar1"]
b1<-arima_obj$coef["ma1"]
plot_T<-T

# This function computes the spectrum of the ARMA-process
spec<-arma_spectrum_func(a1,b1,K,plot_T)$arma_spec

# Fill spec into weight_func: target (first column) and explanatory (second column); both are identical for univariate problems
weight_func<-cbind(spec,spec)
colnames(weight_func)<-c("spectrum target","spectrum explanatory")

# Filter length: number of weights/coefficients of forecast filter
#   Typical economic data (close to random-walk or noise (after differencing)) has 'short' memory
#   Filter-length 10 is suitable for many applications 
#     Exception: seasonal data may require L to go back to yearly lags of the data
L<-10


# Compute one to h-steps ahead forecasts
#   Note that h-step ahead forecasts are 'direct': effective multi-step ahead filters (not iterated one-step rules)
#   This can make a difference when using non-ARMA-based spectra
dfa_forecast<-rep(NA,h)
for (i in 1:h)#i<-2
{
  Lag<--i
# Estimation based on MDFA-MSE wrapper
  mdfa_obj<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 
  
  b<-mdfa_obj$b
# Filter data (apply forecast filter)
  dfa_forecast[i]<-t(b)%*%x[length(x):(length(x)-L+1)]
  
}

# Compare DFA with above predict: both forecasts are virtually indistinguishable
#   For large K the difference vanishes
forecast_comparison<-cbind(c(rep(NA,21),arima_pred),c(rep(NA,21),dfa_forecast),c(x[(len-20):len],rep(NA,h)))
ts.plot(forecast_comparison[,1],ylim=c(min(forecast_comparison,na.rm=T),max(forecast_comparison,na.rm=T)),main="Data (black); classic ARIMA forecast (red) and DFA (green)",col="red")
lines(forecast_comparison[,2],col="green",lty=2)
lines(c(x[(len-20):len],rep(NA,h)),col="black")
abline(v=21)



#-----------------------------------
# Example 4: use non-parametric spectrum (discrete fourier transform: dft) instead of model-based spectrum 


# Use dft
weight_func<-cbind(per(x,T)$DFT,per(x,T)$DFT)
colnames(weight_func)<-c("spectrum target","spectrum explanatory")


K<-nrow(weight_func)-1
# Target forecasting: 
#   -We are interested in all frequencies equally
#   -Target Gamma is a forward-looking allpass filter in the frequency-domain
Gamma<-rep(1,K+1)


# Filter length: number of weights/coefficients of forecast filter
#   -The memory of most economic data is short: thus large L is not meaningful in most cases
#   -Examples: price series are close to random-walks, log-returns are close to noise
L<-10


# Compute one to h-steps ahead forecasts
dfa_forecast<-rep(NA,h)
for (i in 1:h)#i<-2
{
  Lag<--i
  # Estimation based on MDFA-MSE wrapper
  mdfa_obj<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 
  
  b<-mdfa_obj$b
  # Filter data (apply forecast filter)
  dfa_forecast[i]<-t(b)%*%x[length(x):(length(x)-L+1)]
  
}

par(mfrow=c(1,1))
# Compare DFA with above predict: both forecasts are virtually indistinguishable
#   For large K the difference vanishes
forecast_comparison<-cbind(c(rep(NA,21),arima_pred),c(rep(NA,21),dfa_forecast),c(x[(len-20):len],rep(NA,h)))
ts.plot(forecast_comparison[,1],ylim=c(min(forecast_comparison,na.rm=T),max(forecast_comparison,na.rm=T)),main="Data (black); classic ARIMA forecast (red) and DFA (green)",col="red")
lines(forecast_comparison[,2],col="green",lty=2)
lines(c(x[(len-20):len],rep(NA,h)),col="black")
abline(v=21)

# Comment
#   -Forecasts of 'non-parametric' DFA are close to (but visually slightly different than) model-based forecasts
#   -Overfitting (in particular for large L)
#   -Overfitting will be addressed in a later tutorial about regularization

#-----------------------------------------------------------
# Example 5: we compare the non-parametric DFA, relying on the dft, to the best possible forecast approach in an out-of-sample experiment
#   Specifically, for each realization we compute out-of-sample forecasts of non-parametric DFA and of best approach and take the mean of the squared forecast error
# Note that 
#   1. DFA based on non-parametric dft does not assume any 'a priori' knowledge  
#   2. This framework favors the arma-based approach 
#     -we assume knowledge of the true data-generating process (true model-orders) for the arma-forecast
#     -clearly, under these conditions, arma is the best possible approach (inherited from maximum likelihood concept)
#     -So we have a nice benchmark for the non-parametric DFA: if it performs well, then the dft is a meaningful statistic

# Example 5.1: MA(1)
a1<-0.
b1<-0.7
# Example 5.2: ARMA with positive acf
a1<-0.6
b1<-0.7
# Example 5.3: AR with negative acf
a1<--0.9
b1<-0
# Add any other processes...
# We generate anzsim realizations of length len of the arma-process
set.seed(1)
len<-300
mse_true_arma<-mse_dfa<-NULL
# Number of simulations
anzsim<-500

pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
# Loop through all simulations and collect out-of-sample forecast performances
for (i in 1:anzsim)
{

  x<-arima.sim(n=len,list(ar=a1,ma=b1))
# Use in-sample span for model-estimation and for dft  
  x_insample<-x[1:(len-1)]
# True model: estimate model-parameters by relying on classic arima-function
  arima_true_obj<-arima(x_insample,order=c(1,0,1),include.mean=F)
  # Compute forecasts
  arima_true_pred<-predict(arima_true_obj,n.ahead=1)$pred
# Use in-sample span for dft  
  weight_func<-cbind(per(x_insample,F)$DFT,per(x_insample,F)$DFT)
  colnames(weight_func)<-c("spectrum target","spectrum explanatory")
  K<-nrow(weight_func)-1
# Allpass target  
  Gamma<-rep(1,K+1)
# Default filter length for forecasting
  L<-10
# One-step ahead
  Lag<--1
# Compute MSE-filter
  mdfa_obj<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 
  b<-mdfa_obj$b
# Compute out-of-sample forecast  
  dfa_forecast<-t(b)%*%x_insample[length(x_insample):(length(x_insample)-L+1)]
# Mean-square out-of-sample forecast errors  
  mse_true_arma<-c(mse_true_arma,(x[length(x)]-arima_true_pred)^2)
  mse_dfa<-c(mse_dfa,(x[length(x)]-dfa_forecast)^2)
  setTxtProgressBar(pb, i)
}

# Compute the ratio of root mean-square forecast errors:
#   The ratio cannot be larger than 1 asymptotically because our particular design distinguishes arma as the universally best possible design
# Results: 
#   -for L=10, the ratio is typically 97% or larger: in the mean the non-parametric DFA performs as well (by all practical means) as the best possible forecast approach
#     Note that L=10 is fine when forecasting most 'typical' economic data (at least when data is seasonally adjusted)
#   -for L=100 the ratio drops to about 80% 
#     Quantification of (massive) overfitting: similar to fitting an AR(100)-model to the data (Burg max-entropy spectral estimate)
sqrt(mean(mse_true_arma)/mean(mse_dfa))



#-----------------------------------------------------------
# Example 6: we compare two spectral estimates when plugged as weighting functions into DFA
#   1. non-parametric (dft) and 
#   2. Burg's max-entropy (AR-) spectral estimate
#     The idea is that we do not need to identify a model for the data
#     Instead we just fit an AR(p)-model where p is 'sufficiently large' (in the code below we set p=L)
#     Based on the AR(p) model we can derive a (AR-based) spectrum which is used for DFA (instead of dft)


# Example 6.1: MA(1)
a1<-0.
b1<-0.7
# Example 6.2: ARMA with positive acf
a1<-0.6
b1<-0.7
# Example 6.3: AR with negative acf
a1<--0.9
b1<-0
# Add any other processes...
# We generate anzsim realizations of length len of the arma-process
set.seed(1)
len<-300
mse_burg<-mse_dfa<-NULL
# Number of simulations
anzsim<-500

pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
# Loop through all simulations and collect out-of-sample forecast performances
for (i in 1:anzsim)
{
  
  x<-arima.sim(n=len,list(ar=a1,ma=b1))
  # Use in-sample span for model-estimation and for dft  
  x_insample<-x[1:(len-1)]
  # Use in-sample span for dft  
  weight_func<-cbind(per(x_insample,F)$DFT,per(x_insample,F)$DFT)
  colnames(weight_func)<-c("spectrum target","spectrum explanatory")
  K<-nrow(weight_func)-1
  # Allpass target  
  Gamma<-rep(1,K+1)
  # Default filter length for forecasting
  L<-10
  # One-step ahead
  Lag<--1
  # Compute MSE-filter
  mdfa_obj<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 
  b<-mdfa_obj$b
  # Compute out-of-sample forecast  
  dfa_forecast<-t(b)%*%x_insample[length(x_insample):(length(x_insample)-L+1)]
# Burg max-entropy spectrum
#   We fit an AR(L)-model to the data and derive the spectrum from this AR(L)-model
# Fit  
  arima_burg_obj<-arima(x_insample,order=c(L,0,0),include.mean=F)
  ar_burg<-arima_burg_obj$coef
  ma_burg<-NULL
# Note that we can set K_burg to any (positive) value since we use a model-based (AR(L)-) spectrum  
  K_burg<-600
  plot_T<-F
# Derive spectrum from AR(L) model
  weight_func_burg<-arma_spectrum_func(ar_burg,ma_burg,K_burg,plot_T)$arma_spec
  weight_func_burg<-cbind(weight_func_burg,weight_func_burg)
  colnames(weight_func_burg)<-c("spectrum target","spectrum explanatory")
  # Allpass target  
  Gamma_burg<-rep(1,K_burg+1)
  # Filter length for Burg's estimate
  # In principle this could be larger than for dft (above) because AR-spectrum is smoother (less noisy than dft)
  L_burg<-10
  # One-step ahead
  Lag<--1
  # Compute MSE-filter
  mdfa_burg_obj<-MDFA_mse(L_burg,weight_func_burg,Lag,Gamma_burg)$mdfa_obj 
  b_burg<-mdfa_burg_obj$b
  # Compute out-of-sample forecast  
  dfa_burg_forecast<-t(b_burg)%*%x_insample[length(x_insample):(length(x_insample)-L+1)]
  # Mean-square out-of-sample forecast errors  
  mse_burg<-c(mse_burg,(x[length(x)]-dfa_burg_forecast)^2)
  mse_dfa<-c(mse_dfa,(x[length(x)]-dfa_forecast)^2)
  setTxtProgressBar(pb, i)
}

# Compute the ratio of root mean-square forecast errors:
#   The ratio cannot be larger than 1 asymptotically because our particular design distinguishes arma as the universally best possible design
# Results: 
#   -for L=10 the ratio is virtually 1 i.e. performances of both approaches are indistinguishable (up to random sampling error) 
#   -we conclude that dft (as a spectral estimate in DFA) is as good as maximum entropy spectral estimate
sqrt(mean(mse_burg)/mean(mse_dfa))



#-------------------------------------------------------------------------
# Wrap-up: what did we learn
# DFA can replicate classic (MSE) one- and multi-step ahead forecasting by 
#   -specifying a corresponding target (allpass) and forecast horizon (Lag) 
#   -specifying a corresponding spectral estimate
# A non-parametric DFA based on the dft performs nearly as well as the true model
