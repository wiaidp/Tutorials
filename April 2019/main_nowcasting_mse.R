# Todos
#   -Add/discuss constraints


# Purpose of tutorial: 
#   In first tutorial we applied DFA to forecasting and we replicated classic SOTA time series approaches
#   Here we apply DFA to signal extraction (specifically: nowcasting of ideal lowpass) 
#     The examples could be straightforwardly extended to arbitrary targets and/or forecasting/backcasting of signals
#   We know how to specify target (Gamma) and weighting-function (spectrum) from previous tutorial (main_DFA_interface.r)
#   We here compare MSE-performances of best possible one-sided filter (assuming knowledge of the true model) vs. target filter
#   We compare best possible one-sided DFA with empirical DFA based on dft, in- and out-of-sample
#   In the latter case, we analyze overfitting (use various filter-lengths L) by
#     1. Comparing in- and out-of-sample performances
#     2. Analyzing important filter features: coefficients, amplitude, time-shift
# For illustration we here analyze univariate unconstrained MSE designs only
#   -Multivariate examples will be considered in next tutorial
#   -Customization and regularization will be tackled in separate tutorials


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
# Spectrum AR(1): try various processes
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
periodicity<-5
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
L<-periodicity*4

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
ts.plot(output_ideal[anf:enf],col="blue",main="Output of ideal lowpass (blue) vs DFA (red)")
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
