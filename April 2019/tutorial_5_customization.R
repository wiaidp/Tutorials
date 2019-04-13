# Todos
#   Link to trilemma paper




# Purpose of tutorial: 
# -Illustrate ATS-trilemma (see McElroy/Wildi)
# -Analyze filter characteristics (amplitude/shift) when emphasizing Timeliness (the T in the ATS-trilemma)
# -Analyze filter characteristics (amplitude/shift) when emphasizing Smoothness (the S in the ATS-trilemma)
# -Analyze filter characteristics (amplitude/shift) when emphasizing both T and S: power of trilemma (when compared to classic MSE-dilemma)
# -Compare univariate customized DFA to bivariate leading indicator (MSE-) MDFA of previous tutorial
# -Compare all designs to customized bivariate MDFA
# -Proceed to extensive (out-of-sample) simulation studies

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
source("Common functions/play_with_bivariate.r")
source("Common functions/functions_trilemma.r")
source("Common functions/compute_customized_designs.r")
#---------------------------------------------------------------------------------------------
# Example 1: emphasizing Timeliness
# Univariate

len<-300
a1<-0.0
b1<-NULL
set.seed(1)
x<-as.vector(arima.sim(n=len,list(ar=a1,ma=b1)))
# Spectrum
K<-600
plot_T<-T
weight_func<-cbind(arma_spectrum_func(a1,b1,K,plot_T)$arma_spec,arma_spectrum_func(a1,b1,K,plot_T)$arma_spec)
# Target: periodicity
periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# Specify filter length
L<-50
# Nowcast
Lag<-0
# Specify the lambdas: emphasize Timeliness
lambda_vec<-c(0,2^(0:7))
# Specify the fixed eta: do not emphasize Smoothness
eta_vec<-rep(0,length(lambda_vec))

for (i in 1:length(lambda_vec))#i<-1
{
  lambda<-lambda_vec[i]
  eta<-eta_vec[i]
  
  mdfa_obj<-MDFA_cust(L,weight_func,Lag,Gamma,cutoff,lambda,eta)$mdfa_obj
# Keep track of transferfunctions, filter coefficients and of filtered outputs 
  if (i==1)
  {
    trffkt_mat<-mdfa_obj$trffkt
    b_mat<-mdfa_obj$b
    yhat_mat<-filt_func(x,mdfa_obj$b)$yhat
      
  } else
  {
    trffkt_mat<-cbind(trffkt_mat,mdfa_obj$trffkt)
    b_mat<-cbind(b_mat,mdfa_obj$b)
    yhat_mat<-cbind(yhat_mat,filt_func(x,mdfa_obj$b)$yhat)
  } 
}

# Compare amplitude and time shifts
par(mfrow=c(2,2))
mplot<-abs(trffkt_mat)
dimnames(mplot)[[2]]<-paste("Amplitude (",lambda_vec,",",eta_vec,")",sep="")
ax<-rep(NA,ncol(mplot))
ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
plot_title<-"Amplitude functions"
insamp<-1.e+90
title_more<-dimnames(mplot)[[2]]
colo<-rainbow(ncol(mplot))
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
mplot<-abs(trffkt_mat)
# Scale amplitudes for easier visual inspection
for (i in 1:ncol(mplot))
  mplot[,i]<-mplot[,i]/max(mplot[,i])
dimnames(mplot)[[2]]<-paste("Scaled amplitude (",lambda_vec,",",eta_vec,")",sep="")
ax<-rep(NA,ncol(mplot))
ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
plot_title<-"Scaled amplitude functions"
insamp<-1.e+90
title_more<-dimnames(mplot)[[2]]
colo<-rainbow(ncol(mplot))
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
mplot<-Arg(trffkt_mat)/(pi*(0:K)/K)
dimnames(mplot)[[2]]<-paste("Time-shifts (",lambda_vec,",",eta_vec,")",sep="")
ax<-rep(NA,ncol(mplot))
ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
plot_title<-"Time-shifts"
insamp<-1.e+90
title_more<-dimnames(mplot)[[2]]
colo<-rainbow(ncol(mplot))
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)


# Compare filter outputs
anf<-250
enf<-300
par(mfrow=c(1,2))
mplot<-yhat_mat[anf:enf,]
dimnames(mplot)[[2]]<-paste("Filter outputs (",lambda_vec,",",eta_vec,")",sep="")
ax<-rep(NA,ncol(mplot))
ax[1+(0:6)*((nrow(mplot)-1)/6)]<-as.integer(1+(0:6)*((nrow(mplot)-1)/6))
plot_title<-"Filter outputs"
insamp<-1.e+90
title_more<-dimnames(mplot)[[2]]
colo<-rainbow(ncol(mplot))
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)


# Compare scaled filter outputs
mplot<-yhat_mat[anf:enf,]
# Scale outputs for easier visual inspection
for (i in 1:ncol(mplot))
  mplot[,i]<-mplot[,i]/max(abs(trffkt_mat)[,i])
dimnames(mplot)[[2]]<-paste("Filter outputs (",lambda_vec,",",eta_vec,")",sep="")
ax<-rep(NA,ncol(mplot))
ax[1+(0:6)*((nrow(mplot)-1)/6)]<-as.integer(1+(0:6)*((nrow(mplot)-1)/6))
plot_title<-"Scaled filter outputs"
insamp<-1.e+90
title_more<-dimnames(mplot)[[2]]
colo<-rainbow(ncol(mplot))
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)



#-------------------------------------------------------------------------
# Example 2: emphasizing Smoothness
# Univariate

# Specify the lambdas: do not emphasize Timeliness
lambda_vec<-rep(0,length(eta_vec))
# Specify the fixed eta: emphasize Smoothness
eta_vec<-0.3*0:6

# The following function replicates the above lengthy code of example 1
compute_customized_designs_func(lambda_vec,eta_vec,L,weight_func,Lag,Gamma,cutoff)
  

#-----------------------------------------------------------------------------------------------
# Example 3: emphasizing both Timeliness and Smoothness
# Univariate

# MSE vs. customized emphasizing S&T

eta_vec<-c(0,0.5)
lambda_vec<-c(0,30)

# The following function replicates the above lengthy code of example 1
compute_customized_designs_func(lambda_vec,eta_vec,L,weight_func,Lag,Gamma,cutoff)

#--------------------------------------------------------------------------------
# Example 4: compare bivariate leading indicator and univariate customized
#   See previous tutorial for a background to the play_bivariate_func function below
#     This function sets-up a MDFA-experiment based on a bivariate noisy leading-indicator design
# Weak autocorrelation (typical for log-returns of (positive) economic data)
a1<-0.08
# Noisy but not too much so
scale_idiosyncratic<-0.4
# Fairly long in-sample span
len<-300
# Modest filter length (periodicity is 6 so L=12 is needed to damp components in the stopband)
L<-12

play_obj<-play_bivariate_func(a1,scale_idiosyncratic,len,L)
  
play_obj$b__bivariate
play_obj$perf_mse





