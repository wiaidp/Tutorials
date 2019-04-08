# Todos
# Specify target in frequency-domain and time-domain
# Use different spectrum: flat, AR (discuss fit), customized/tweaked




# In previous tutorial (DFA and forecasting) we emphasized a particular target: Gamma was an allpass filter
# Here we propose generic targets, including lowpass, bandpass, highpass, Hodrick-Prescott, arbitrary,....

# Purpose of tutorial: illustrate (M)DFA user-interface
# Design: univariate, MSE (no customization), unconstrained (no regularization)
#   -Customization and regularization will be tackled in separate tutorials



# Start from scratch
rm(list=ls())



library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)

# Gamma and weight_func in the head of MDFA_mse allow for powerful interaction of the user with the estimation algorithm
head(MDFA_mse)

#-----------------------------------------------------------------------------------------------
# Source common functions

source("Common functions/plot_func.r")
source("Common functions/arma_spectrum.r")
source("Common functions/ideal_filter.r")


#-------------------------------------------------------------------------------------------------
# Play with target of DFA
#   DFA requires the user to supply a 'target' and a 'spectrum'; DFA then returns an optimal estimate 
#     -Optimality: Mean-square error (MSE) or beyond MSE (ATS-trilemma)
# In the following examples: illustration of various targets
#   This is more general than state-of-the-art forecast approaches
# Design:
#   MSE
#   Univariate
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

#-------------------
# Example 1: classic one step ahead forecasting (replicates SOTA ARMA-models)
# Target forecasting: 
#   -We are interested in all frequencies equally
#   -Target Gamma is a forward-looking allpass filter in the frequency-domain
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

# The following function plots 
#   -the coefficients
#   -the time-shift and
#   -the amplitude of the filter
plot_estimate_func(mdfa_obj,weight_func,Gamma)

# Comments
#   -The amplitude function of the MSE-filter is close to zero, as desired 
#   -It is not exactly zero because K is finite
#   -For K/L>10 mismatch is negligible by all practical means 

#--------------------------------------------------------------------------------------
# Example 2: lowpass
#   In contrast to forecasting we are not interested in all frequency-components
#   When specifying a lowpass design, the user generally wants to fade-out undesirable (high-frequency) 'noise'
#   Alternatively, the user might want to emphasize interesting low-frequency components: for example trend and/or cycle

# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
plot(Gamma,type="l",main=paste("Ideal lowpass, periodicity=",periodicity,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Comments
#   The target filter Gamma lets pass through (does not alter) components with durations larger than 2*periodicity
#   The target filter Gamma stops (eliminates) components with durations smaller than 2*periodicity
#   Filter is called 'ideal' lowpass
#   This filter is very easy to specify in the frequency-domain and it is intuitively appealing

# Since the transferfunction Gamma is real-valued, the filter is symmetric (imaginary parts of symmetric filters cancel)
# We now apply the ideal filter to white noise
len<-1000
set.seed(1)
x<-rnorm(len)
# Length of filter (true ideal filter is bi-infinite)
M<-100
# This function computes the filter coefficients and applies the filter to x
id_obj<-ideal_filter_func(periodicity,M,x)

# The coefficients are symmetric: the filter is not causal (it needs future data)  
ts.plot(id_obj$gamma)
# The filtered series (red line in plot below) is smooth: high-frequency noise has been damped
ts.plot(x)
lines(id_obj$y,col="red")
# Mean duration between consecutive zero-crossings of the filtered data (mean holding-time of an asset in trading with this filter)
#   Higher periodicities imply stronger smoothing and thus longer holding-times (less frequent trades)
id_obj$mean_holding_time

# Idea: the user can specify a periodicity which matches his particular interests (more flexible than classic time series approaches)
# Problem: the ideal filter cannot be applied in 'real-time' because it is not causal
# Solution: approximate Gamma by a real-time (causal) filter i.e. use DFA (see examples below)

#--------------------------------------------------------------------------------------
# Example 3: bandpass
#   In contrast to forecasting we are not interested in all frequency-components
#   When specifying a bandpass design, the user generally wants to fade-out undesirable (high- and low-frequency) components
#   Alternatively, the user might want to emphasize interesting cycle components
#   This type of filter is typically in business-cycle analysis, for example

# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity_low<-5
cutoff_low<-pi/periodicity_low
periodicity_high<-10
cutoff_high<-pi/periodicity_high
Gamma<-((0:(K))<=K*cutoff_low/pi)-((0:(K))<=K*cutoff_high/pi)+1.e-9


plot(Gamma,type="l",main=paste("Ideal lowpass, periodicities=",periodicity_low,",",periodicity_high,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Comments
#   The target filter Gamma lets pass through (does not alter) components with durations between 2*periodicity_low and 2*periodicity_high
#   The target filter Gamma stops (eliminates) all other components
#   Filter is called 'ideal' bandpass
#   This filter is very easy to specify in the frequency-domain and it is intuitively appealing

# The coefficients of the ideal bandpass can be calculated by using the previous function
#   The coefficients of the bandpass are the difference of the two lowpass filters
# We skip this calculation since it is not of immediate interest

# As for the ideal lowpass, the ideal bandpass is not causal and can be approximated by DFA



#--------------------------------------------------------------------------------------
# Example 4: more complex targets


# Target: specify cutoff=pi/periodicity of lowpass ideal target

Gamma<-rep(1,K+1)

Gamma[3:20]<-0.4
Gamma[37:98]<-0
Gamma[167:208]<-pi
Gamma[398:476]<-0.1

plot(Gamma,type="l",main="Arbitrary target",
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Comments
# As for the previous ideal lowpass and bandpass the above generic filter is not causal and can be approximated by DFA
# The user can specify any target that fits his research interest: more general than classic time series approaches


#-------------------------------------------------------------------------------------------
# Example 5: approximate non-causal ideal lowpass by causal filter
#   Interpret amplitude and time-shift functions

K<-600
# Spectrum: white noise assumption
#  First column is target, second column is explanatory variable: in a univariate design target and explanatory are the same
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
colnames(weight_func)<-c("target","explanatory")

periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
plot(Gamma,type="l",main=paste("Ideal lowpass, periodicity=",periodicity,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-200

# Estimation based on MDFA-MSE wrapper
mdfa_obj_mse<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 

b<-mdfa_obj_mse$b
plot_estimate_func(mdfa_obj_mse,weight_func,Gamma)



# Comments
#   1. The filter coefficients in the first plot are one-sided (causal filter)
#   2. Since the filter does not look into the future, the output will be delayed (with respect to the non-causal target)
#     The corresponding delay/lag can be seen in the second plot (time-shift)
#   3. The fit of the target (violet line last plot) by the (DFA-) amplitude function can be 
#     seen in the last plot
#   4. Amplitude and shift differ from target


len<-1000
set.seed(1)
x<-rnorm(len)
# Length of filter (true ideal filter is bi-infinite)
M<-100
# This function computes the filter coefficients and applies the filter to x
id_obj<-ideal_filter_func(periodicity,M,x)

# The coefficients are symmetric: the filter is not causal (it needs future data)  
ts.plot(id_obj$gamma)
# The filtered series (red line in plot below) is smooth: high-frequency noise has been damped
ts.plot(x)
lines(id_obj$y,col="red")
# Mean duration between consecutive zero-crossings of the filtered data (mean holding-time of an asset in trading with this filter)
#   Higher periodicities imply stronger smoothing and thus longer holding-times (less frequent trades)
id_obj$mean_holding_time

# Idea: the user can specify a periodicity which matches his particular interests (more flexible than classic time series approaches)





#-----------------------------------
# Example 1.2: classic one- and multi-step ahead forecasting ARMA-process (replicates SOTA ARMA-models)
# Target forecasting: 
#   -We are interested in all frequencies equally
#   -Target Gamma is a forward-looking allpass filter in the frequency-domain
Gamma<-rep(1,K+1)
plot(Gamma,type="l",main=paste("Allpass forecast target, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

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

# One step ahead: Lag=-1
Lag<-5
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
# 1. Forecasting: Lag<0
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


#-----------------------------------
# Example 1.3: compare DFA with classic arima R-code

# Specify forecast horizon

h<-10

# Data: simulate ARMA (select any coefficients)

a1<-0.6
b1<-0.7
set.seed(1)
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
lines(forecast_comparison[,2],col="green")
lines(c(x[(len-20):len],rep(NA,h)),col="black")
abline(v=21)



#-------------------------------------------------------------------------
# Wrap-up: what did we learn
# DFA can replicate classic (MSE) one- and multi-step ahead forecasting by 
#   -specifying a corresponding target (allpass) and forecast horizon (Lag) 
#   -specifying a corresponding spectral estimate
 
