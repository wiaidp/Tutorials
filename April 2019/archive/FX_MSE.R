rm(list=ls())

# Data
# Daily EURUSD has been included in MDFA-package (18-02-2017)
# But one can refresh the data from Quandl by selecting refresh_data<-T
#   Just plug-in the API key in the corresponding code line

refresh_data<-F


# Libraries: Quandl is required for refreshing the data
if (refresh_data)
  library(Quandl)

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)



#------------------------------------------------------------------
# Example 1: forecast based on white noise spectrum (data is assumed to be zero-mean white noise)
#   Forecast: Gamma=1 for all frequencies
#   Try various Lag (Lag=0,1,2 and Lag=-1,-2,...)
#   Try various K (denseness of frequency grid)
#   Look at the DFA criterion value

# Denseness of discret frequency grid: 
# This is an arbitrary integer: it means that omega_k=k*pi/K for k=0,1,...,K in the MDFA-criterion
#   The larger K, the better the approximation of the discrete sum (MDFA-criterion) to the integral (theoretical MSE)
#   The larger K, the longer the computation time
K<-240
# Filter length
#   It is assumed that L<= 2*K (otherwise the problem is ill-conditioned i.e. matrices cannot be inverted anymore)
L<-22
L<-min(2*K,L)
# One-step ahead forecast
Lag<--1
# cutoff: allpass target (forecasting means that we target everything i.e. Gamma is an allpass filter)
cutoff<-pi
# Target 
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9

plot(Gamma,type="l",main="Target",
     axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                     "4pi/6","5pi/6","pi"))
axis(2)
box()

# White noise spectrum: constant in [0,pi] (the value of the constant is irrelevant)
weight_func_white_noise<-matrix(rep(1,2*(K+1)),ncol=2)
tail(weight_func_white_noise)

# Estimation based on MDFA-MSE wrapper
mdfa_obj_mse<-MDFA_mse(L,weight_func_white_noise,Lag,Gamma)$mdfa_obj 

# Coefficients
tail(mdfa_obj_mse$b)
# Criterion: should be close to one of the most mythic numbers in mathematics (approximation improves when increasing K above)
mdfa_obj_mse$MS_error  
# Plot of amplitude
plot(abs(mdfa_obj_mse$trffkt),type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()


# This piece replicates the above outcome with the generic mdfa_analytic (it requires installation of the file control_default.r for setting all hyper-parameters)
if (F)
{
  path.pgm <- paste(getwd(),"/",sep="")
# Call default MSE-parameters  
  source(file=paste(path.pgm,"control_default.r",sep=""))
# By default Lag=0
  # One-step ahead forecast
  Lag<--1
# Estimate filter coefficients:
  mdfa_obj<-mdfa_analytic(L, lambda, weight_func_white_noise, Lag, Gamma, eta, cutoff, i1,i2, weight_constraint, lambda_cross, lambda_decay, lambda_smooth,lin_eta, shift_constraint, grand_mean, b0_H0, c_eta, weight_structure,white_noise, synchronicity, lag_mat, troikaner) 
# Filter coefficients: compare MDFA and previous DFA
  mdfa_obj$b
  
  plot(abs(mdfa_obj$trffkt),type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
       axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
  axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
  axis(2)
  box()
}

#-----------------------------------------------------------------------------------------------------------------
# Example 2: same as above but data is not assumed to be centered (white noise + non-vanishing mean)
#   Non-wanishing mean can be tackled by specifying a large value of weight_func_white_noise in frequency zero
#   Alternative: in the more general MSE-wrapper MDFA_mse_constraint one could set i1<-T 

# White noise spectrum
weight_func_white_noise_new<-weight_func_white_noise
# Frequency zero receives a large weight
#   Too large a value results in near singularity i.e. numerical error (be careful and stay realist)
#   Trick: 1000*max(weight_func_white_noise) 
weight_func_white_noise_new[1,]<-1000*max(abs(weight_func_white_noise[1,]))

# Estimation
mdfa_obj_mse<-MDFA_mse(L,weight_func_white_noise_new,Lag,Gamma)$mdfa_obj 

# Coefficients
b<-mdfa_obj_mse$b
ts.plot(b)

# Criterion
mdfa_obj_mse$MS_error  
# Plot of amplitude
plot(abs(mdfa_obj_mse$trffkt),type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

#-----------------------------------------------------------------------------------------------------------------
# Filter series: apply filter to log-returns of FX 

returns<-diff(log(FX))

par(mfrow=c(2,1))
plot(FX,main=paste(FX_pair,": ",index(FX)[1],"/",index(FX[length(FX)]),sep=" "))
plot(returns,main=paste("Log-retruns",FX_pair,": ",index(FX)[1],"/",index(FX[length(FX)]),sep=" "))

# We need ordinary time series (not xts) for filtering because of lag-structure when computing output
x<-as.double(returns)
names(x)<-index(FX)

# Filtering
yhat<-rep(NA,length(x))
for (i in L:length(x))
  yhat[i]<-sum(b*x[i:(i-L+1)])
names(yhat)<-names(x)

# Plot
par(mfrow=c(1,1))
plot(returns,main="Returns (black) and filtered series (red)")
lines(as.xts(yhat),col="red")


#--------------------------------------------------------------------------------------------------
# Compute performances: sign-rule and signal weighting
par(mfrow=c(2,1))
# Trading performance: sign-rule (don't forget to lag the signal by one day)
cum_log_perf<-cumsum(na.omit(lag(sign(as.xts(yhat)),1)*returns))
plot(cum_log_perf,main="Signum rule")

# Signal weigthing
#   We weight the trade by the filter-output (proxi for strength of signal (plausible indicator for confidence))
#   Normalize trading rule by scale of signal: divide by sqrt(var(yhat,na.rm=T)) (otherwise scale of performance is not meaningful)
cum_log_perf<-cumsum(na.omit(lag(as.xts(yhat),1)*returns))/sqrt(var(yhat,na.rm=T))
plot(cum_log_perf,main="Signal weighting")

#--------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
# Example 3: we look at a random-walk spectrum

# Random-walk
spec_random_walk<-1/(1-exp(1.i*(0:K)*pi/K))
# Remove singularity in frequency zero (otherwise computations stop by an error)
#   Just plug in a large number (1000*max)
spec_random_walk[1]<-1000*max(abs(spec_random_walk[-1]))
# Define the spectrum for target and for explanatory variables which are identical in a univariate setting 
#   We use absolute values of spec_random_walk for simplicity (we could leave complex values)
#   Note that this is not the pseudo-spectral density (which would correspond to the square)
#     In dfa_mse we would plug-in abs(spec_random_walk)^2 (which corresponds to the periodogram)
#     In MDFA we plug-in spec_random_walk (which corresponds to the DFT)
weight_func_random_walk<-abs(cbind(spec_random_walk,spec_random_walk))
head(weight_func_random_walk)

# Estimation
mdfa_obj_mse<-MDFA_mse(L,weight_func_random_walk,Lag,Gamma)$mdfa_obj 

# Coefficients
b<-mdfa_obj_mse$b
head(b)

# Criterion
mdfa_obj_mse$MS_error  
# Plot of amplitude
plot(abs(mdfa_obj_mse$trffkt),type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()


#------------------------------------------------------------------------------------------
# Example 4: General setting or linear combination of random-walk (example 3) and uncentered noise (example 2)

# Take any positive weight in [0,1]
alpha<-0.8
# Ensure that the weight is in [0,1]
alpha<-max(0,alpha)
alpha<-min(1,alpha)
# New general weight_func is alpha*noise+random-walk
#   For alpha=0 the random-walk is obtained
#   For alpha=1 the white noise is obtained
weight_func_general<-(alpha*weight_func_white_noise_new+(1-alpha)*weight_func_random_walk)

# Estimation
mdfa_obj_mse_general<-MDFA_mse(L,weight_func_general,Lag,Gamma)$mdfa_obj 

# Coefficients
b<-mdfa_obj_mse_general$b
ts.plot(b,main=paste("alpha=",alpha))

# Criterion
mdfa_obj_mse_general$MS_error  
# Plot of amplitude
plot(abs(mdfa_obj_mse_general$trffkt),type="l",main=paste("Amplitude concurrent, denseness=",K,", alpha=",alpha,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

#--------------------------------------------------------------------------------------
# That's all folks!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#--------------------------------------------------------------------------------------







