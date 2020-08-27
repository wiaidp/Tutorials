# To do: link to optimal filter paper and to companion paper



# Purpose of tutorial: 
#   -In first tutorial we applied DFA to forecasting and we replicated classic SOTA time series approaches
#   -In second tutorial We learned how to specify target (Gamma) and weighting-function (spectrum) 
#     (play with interface, see main_DFA_interface.r); we learned how to interpret amplitude and time-shift functions; 
#     we analyzed overfitting and proposed some easy diagnostics (rippled amplitude, shifts, coefficients, degradation of out-of-sample performances)
#   -Here we apply DFA to signal extraction (specifically: nowcasting of ideal lowpass) 
#     The examples could be straightforwardly extended to arbitrary targets and/or forecasting/backcasting of signals
#   We compare MSE-performances of best possible one-sided filter (assuming knowledge of the true model) vs. target filter
#   We compare best possible one-sided DFA with empirical DFA based on dft as well as DFA based on on Burg's maximum entropy 
#     estimate, both in-sample and out-of-sample
#   We show that lowpass nowcasting with DFA based on discrete fourier transform (dft) performs nearly as well (in terms of MSE performances) as best possible approach (assuming knowledge of true data generating process)
#     -Assuming some elementary care is taken in order to avoid overfitting
#     -Heavily overparametrized designs still perform honorably when compared to best possible approach (10%loss of efficiency 'only')
#   We also show that DFA based on dft and DFA based on Burg's max-entropy spectral estimate perform equally
#     These results suggest that dft or Burg-spectrum can be used for real-time signal extraction (lowpass nowcasting) in practice

# For illustration we here analyze univariat unconstrained MSE designs only
#   -Multivariate examples will be considered in next tutorial
#   -Customization and regularization will be tackled in separate tutorials

# Functions and parameters: throughout this tutorial we rely on the same function MDFA_mse as used previously and we learn how to (better) understand 
# its parameters in a univariate signal extraction (nowcasting) framework. 



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
source("Common functions/ideal_filter.r")
source("Common functions/mdfa_trade_func.r")

#---------------------------------------------
# Example 1 ARMA-process: best possible MSE estimate (assuming knowledge of true model)
#   MSE
#   Univariate

# K defines the frequency-grid for estimation: all frequencies omega_j=j*pi/K, j=0,1,2,...,K are considered
K<-600
# Spectrum ARMA(1,1): try various processes
#   Log-returns of typical economic data close to noise i.e. a1~0
a1<-0.1
b1<-NULL
plot_T<-T
# This function computes the spectrum of the ARMA-process
spec<-arma_spectrum_func(a1,b1,K,plot_T)$arma_spec
# Fill into weight_func: target (first column) and explanatory (second column); both are identical for univariate problems
weight_func<-cbind(spec,spec)
colnames(weight_func)<-c("spectrum target","spectrum explanatory")
# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity<-10
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length: L/K should be 'small'
L<-100

# Estimation based on MDFA-MSE wrapper
mdfa_obj_mse<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 

plot_estimate_func(mdfa_obj_mse,weight_func,Gamma)

# Apply filter to data
set.seed(1)
len<-1000
x<-as.vector(arima.sim(n=len,list(ar=a1,ma=b1)))
# Length of filter (true ideal filter is bi-infinite)
M<-100
# This function computes the filter coefficients and applies the filter to x
id_obj<-ideal_filter_func(periodicity,M,x)
output_ideal<-id_obj$y


# DFA output
b<-mdfa_obj_mse$b
filt_obj<-filt_func(x,b)

output_dfa<-filt_obj$yhat
output_dfa[1:(L-1)]<-NA

# The filtered series (red line in plot below) is smooth: high-frequency noise has been damped
ts.plot(output_ideal,col="blue",main="Output of ideal lowpass (blue) vs DFA (red)")
lines(output_dfa,col="red")

# Comments: output of ideal filter (blue) is not available towards sample-end: DFA computes a 'best' (MSE) one-sided filter which runs till sample-end

# Let's now zoom into the plot and have a closer look at both series
anf<-400
enf<-500
ts.plot(output_ideal[anf:enf],col="blue",main="Output of ideal lowpass (blue) vs DFA (red)")
lines(output_dfa[anf:enf],col="red")
abline(h=0)


#---------------------------------------------
# Example 2 ARMA-process: as above but rely on dft for estimating spectrum 
# Try various L and observe in-sample and out-of-sample performances in the table at the end of the example
#   1. L<-periodicity*2 (minimal length for removing a component with that periodicity)
#   2. L<-periodicity*4  (a1=0.9 or a1=0: amplitude still OK, out-of-sample OK too)
#   3. L<-periodicity*8 (a1=0.9: maybe we can observe start of overfitting... Look at amplitude for instance)
#   4. L<-periodicity*16 (overfitting: more severe when a1=0.9 (good news since log-returns typical economic data close to noise...))



in_sample<-300
# Use in-sample span for dft  
x_insample<-x[1:in_sample]
# Use in-sample span for dft  
weight_func_dft<-cbind(per(x_insample,F)$DFT,per(x_insample,F)$DFT)
colnames(weight_func_dft)<-c("spectrum target","spectrum explanatory")
K_dft<-nrow(weight_func_dft)-1
# Allpass target  
Gamma_dft<-(0:(K_dft))<=K_dft*cutoff/pi+1.e-9
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
#  L<-periodicity*2 is required for damping a component with that periodicity
L<-periodicity*2

# Estimation based on MDFA-MSE wrapper
mdfa_obj_dft_mse<-MDFA_mse(L,weight_func_dft,Lag,Gamma_dft)$mdfa_obj 

plot_estimate_func(mdfa_obj_dft_mse,weight_func_dft,Gamma_dft)



# Compare output of (non-causal) ideal filter, best MSE and DFA based on dft 

# DFA output
b_dft<-mdfa_obj_dft_mse$b
filt_dft_obj<-filt_func(x,b_dft)

output_dft_dfa<-filt_dft_obj$yhat
output_dft_dfa[1:(L-1)]<-NA

# The filtered series (red line in plot below) is smooth: high-frequency noise has been damped
ts.plot(output_ideal,col="blue",main="Output of ideal lowpass (blue) vs best possible (red) and DFA-dft (green)")
lines(output_dfa,col="red")
lines(output_dft_dfa,col="green")

# Comments: output of ideal filter (blue) is not available towards sample-end: DFA computes a 'best' (MSE) one-sided filter which runs till sample-end

# Let's now zoom into the plot and have a closer look at both series
anf<-400
enf<-500
ts.plot(output_ideal[anf:enf],col="blue",main="Output of ideal lowpass (blue) vs DFA (red) and DFA-dft (green)")
lines(output_dfa[anf:enf],col="red")
lines(output_dft_dfa[anf:enf],col="green")
abline(h=0)

# Compute Root-MSE-performances
dat_mat<-cbind(output_ideal,output_dfa,output_dft_dfa)
dat_mat_in_sample<-dat_mat[1:in_sample,]
dat_mat_out_sample<-dat_mat[(in_sample+1):nrow(dat_mat),]

# Remove all time points for which ideal filter is not defined (start and end of sample)
dat_mat_in_sample<-na.exclude(dat_mat_in_sample)
dat_mat_out_sample<-na.exclude(dat_mat_out_sample)
colnames(dat_mat_in_sample)<-colnames(dat_mat_out_sample)<-c("ideal","best MSE","dft")
mat_mse_result<-rbind(c(sqrt(mean((dat_mat_in_sample[,"ideal"]-dat_mat_in_sample[,"best MSE"])^2)),
sqrt(mean((dat_mat_in_sample[,"ideal"]-dat_mat_in_sample[,"dft"])^2))),
c(sqrt(mean((dat_mat_out_sample[,"ideal"]-dat_mat_out_sample[,"best MSE"])^2)),
sqrt(mean((dat_mat_out_sample[,"ideal"]-dat_mat_out_sample[,"dft"])^2))))

rownames(mat_mse_result)<-c("in sample","out sample")
colnames(mat_mse_result)<-c("Theoretically best (true model)","DFA based on dft")

mat_mse_result

#----------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# Example 3: Same as previous exercise but we conduct a simulation study based on multiple 
# realizations of an arma-process. This example corresponds to the simulation 
# studies in the previous tutorial (forecasting) but Gamma is now a lowpass (not an allpass) 
# and we emphasize nowcasting (Lag=0) 

# Example 3.1: MA(1)
a1<-0.
b1<-0.7
true_model_order<-c(0,0,1)
# Example 3.2: ARMA with positive acf
a1<-0.6
b1<-0.7
true_model_order<-c(1,0,1)
# Example 3.3: AR with negative acf
a1<--0.9
b1<-0
true_model_order<-c(1,0,0)
# Example 3.4: close to noise (typical for log-returns of FX-data)
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
in_sample<-300
# Frequency grid for DFA based on true model
K_true<-600
# Lowpass target  
periodicity<-10
# Default (reasonable) filter length for nowcasting
L<-2*periodicity
# Nowcast
Lag<-0
# Length of ideal filter
M<-100
mse_true<-mse_dft<-NULL
pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
# Loop through all simulations and collect out-of-sample forecast performances
for (i in 1:anzsim)
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
  # Compute MSE-filter
  mdfa_dft_obj<-MDFA_mse(L,weight_func_dft,Lag,Gamma_dft)$mdfa_obj 
  b_dft<-mdfa_dft_obj$b
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
  # Mean-square out-of-sample filter error  
  mse_true<-c(mse_true,mean((output_ideal-output_dfa_true)[(in_sample+1):len-M]^2,na.rm=T))
  mse_dft<-c(mse_dft,mean((output_ideal-output_dfa_dft)[(in_sample+1):len-M]^2,na.rm=T))
  setTxtProgressBar(pb, i)
}

# Compute the ratio of root mean-square forecast errors:
#   The ratio cannot be larger than 1 asymptotically because our particular design distinguishes arma as the universally best possible design
sqrt(mean(mse_true)/mean(mse_dft))

# Results: 
#   -for L=2*periodicity and in_sample=300, the ratio is typically around 97%: in the mean the non-parametric DFA performs as well (by all practical means) as the best possible forecast approach
#   -for L=2*periodicity and in_sample=100, the ratio is typically around 89%: in the mean the non-parametric DFA performs nearly as well as the best possible forecast approach
#     -Note that L=2*periodicity is fine for damping all components with durations shorter/equal periodicity
#     -But fitting L=2*periodicity=20 parameters for a time series of length in_sample=100 is 
#       not extremely smart (overfitting). See tutorial on regularization...



#---------------------------------------------------------------------------------------
# Example 4: same as above but we compare DFA-dft (as above) and DFA based on Burg's maximum entropy spectral estimate

# Example 4.1: MA(1)
a1<-0.
b1<-0.7
# Example 4.2: ARMA with positive acf
a1<-0.6
b1<-0.7
# Example 4.3: AR with negative acf
a1<--0.9
b1<-0
# Example 4.4: nearly noise
a1<--0.08
b1<-0
# Add any other processes...
# We generate anzsim realizations of length len of the arma-process
set.seed(1)
len<-1000
mse_true_arma<-mse_dfa<-NULL
# Number of simulations
anzsim<-500
# Length of in-sample span
in_sample<-300
# Frequency grid for DFA based on Burg's estimate
K_burg<-600
# Lowpass target  
periodicity<-10
# Default (reasonable) filter length for nowcasting: this is also used for the estimation of Burg's max-entropy spectrum
L<-2*periodicity
# Nowcast
Lag<-0
# Length of ideal filter
M<-100
mse_burg<-mse_dft<-NULL
pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
# Loop through all simulations and collect out-of-sample forecast performances
for (i in 1:anzsim)
{
  # Distinguish white noise  
  if (abs(a1)+abs(b1)>0)
  {
    # Generate series  
    x<-as.vector(arima.sim(n=len,list(ar=a1,ma=b1)))
  } else
  {
    x<-rnorm(len)
    spec<-rep(1,K_burg+1)
  }
  # Use in-sample span for model-estimation and for dft  
  x_insample<-x[1:in_sample]
  # Burg spectral estimate: use AR(L)
  # Restrict length of AR (otherwise numerical optimization fails)
  arima_burg_obj<-arima(x_insample,order=c(min(L,10),0,0),include.mean=F)
  # Spectrum based on burg model  
  spec<-arma_spectrum_func(arima_burg_obj$coef,NULL,K_burg,F)$arma_spec
  weight_func<-cbind(spec,spec)
  colnames(weight_func)<-c("spectrum target","spectrum explanatory")
  weight_func_burg<-weight_func
  cutoff<-pi/periodicity
  # target burg model: frequency grid is not the same as for dft below i.e. K is different
  Gamma_burg<-(0:(K_burg))<=K_burg*cutoff/pi+1.e-9
  mdfa_burg_obj<-MDFA_mse(L,weight_func_burg,Lag,Gamma_burg)$mdfa_obj 
  b_burg<-mdfa_burg_obj$b
  # Use in-sample span for dft  
  weight_func_dft<-cbind(per(x_insample,F)$DFT,per(x_insample,F)$DFT)
  colnames(weight_func_dft)<-c("spectrum target","spectrum explanatory")
  K_dft<-nrow(weight_func_dft)-1
  Gamma_dft<-(0:(K_dft))<=K_dft*cutoff/pi+1.e-9
  # Compute MSE-filter
  mdfa_dft_obj<-MDFA_mse(L,weight_func_dft,Lag,Gamma_dft)$mdfa_obj 
  b_dft<-mdfa_dft_obj$b
  # Filter data
  # 1. ideal filter  
  id_obj<-ideal_filter_func(periodicity,M,x)
  output_ideal<-id_obj$y
  # 2. DFA burg
  filt_burg_obj<-filt_func(x,b_burg)
  output_dfa_burg<-filt_burg_obj$yhat
  # 3. DFA dft
  filt_dft_obj<-filt_func(x,b_dft)
  output_dfa_dft<-filt_dft_obj$yhat
  # Mean-square out-of-sample filter error  
  mse_burg<-c(mse_burg,mean((output_ideal-output_dfa_burg)[(in_sample+1):len-M]^2,na.rm=T))
  mse_dft<-c(mse_dft,mean((output_ideal-output_dfa_dft)[(in_sample+1):len-M]^2,na.rm=T))
  setTxtProgressBar(pb, i)
}

# Compute the ratio of root mean-square forecast errors:

sqrt(mean(mse_burg)/mean(mse_dft))

# Results: 
#   -for L=2*periodicity and in_sample=300, the ratio is 98% i.e. both approaches are indistinguishable by all practical means

