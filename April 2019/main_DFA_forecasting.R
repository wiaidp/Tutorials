# Purpose of tutorial: 
# DFA can replicate classic (MSE) one- and multi-step ahead forecasting by 
#   -specifying a corresponding target (allpass) and forecast horizon (Lag) 
#   -specifying a corresponding spectral estimate


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
colnames(weight_func)<-c("target","explanatory")

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
#   -The amplitude function of the MSE-filter is close to zero, as desired 
#   -It is not exactly zero because K is finite
#   -For K/L>10 mismatch is negligible by all practical means 


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
colnames(weight_func)<-c("target","explanatory")

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
#   -Our solution gets aribtrarily close for increasing K (denser grid)
#   -For K/L>10, mismatch is negligible by all practical means 

#   -For an MA(1)-process and Lag=-1 (one-step ahead) the optimal forecast filter has weights (b1,-(b1^2),b1^3,-(b1^4),...)
#   -For an MA(1)-process and Lag=-2 (two-steps ahead) the optimal forecast filter has weights 0
#   -Our solution gets aribtrarily close for increasing K (denser grid)
#   -For K/L>10 mismatch is negligible by all practical means 

# 2. Nowcasting: Lag=0
#   -For any process the optimal filter has weights 1,0,0,....
#   -This is a trivial estimation problem because our target is an allpass: for lowpass/bandpass/highpass the estimation problem is nomore trivial

# 3. Backcasting: Lag>0
#   -For any process the optimal filter has weight 1 at 'lag' specified by Lag and otherwise 0
#   -This is a trivial estimation problem because our target is an allpass: for lowpass/bandpass/highpass the estimation problem is nomore trivial

# 4. For real-valued Lag the filter interpolates the data between consecutive (future or past) time-points



#-----------------------------------
# Example 3: compare DFA with classic arima R-code

# Specify forecast horizon

h<-5

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
colnames(weight_func)<-c("target","explanatory")

# Filter length: number of weights/coefficients of forecast filter
L<-10


# Compute one to h-steps ahead forecasts
#   Note that h-step ahead forecasts are effective multi-step ahead filters (not iterated one-step)
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
colnames(weight_func)<-c("target","explanatory")


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


#-------------------------------------------------------------------------
# Wrap-up: what did we learn
# DFA can replicate classic (MSE) one- and multi-step ahead forecasting by 
#   -specifying a corresponding target (allpass) and forecast horizon (Lag) 
#   -specifying a corresponding spectral estimate
 
