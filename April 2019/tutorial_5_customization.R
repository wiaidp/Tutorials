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

eta_vec<-c(0,0.4)
lambda_vec<-c(0,50)

# The following function replicates the above lengthy code of example 1
compute_customized_designs_func(lambda_vec,eta_vec,L,weight_func,Lag,Gamma,cutoff)

#------------------------------------------------------------------------------------------------
# Example 4: the following function was used in the trilemma paper (McElroy/Wildi) 
#   It computes the full set of ATS-statistics and it performs complete simulation runs with detailed in-sample/out-of-sample results
#   It relies on DFA-code (which is replicated one-to-one by MDFA-univariate but which is much simpler/shorter)

eta_vec<-c(0,1.8)
# Specify the fixed lambda
lambda_vec<-c(0,128)

# Specify the processes: ar(1) with coefficients -0.9,0.1 and 0.9
#   One can specify different processes at once and the code will loop through 
#   For simplicity we here use positive and nearly zero autocorrelation
a_vec<-c(0.9,0.08,-0.9)
# Ordinary ATS-components
scaled_ATS<-F
# Generate a single realization of the processes
anzsim<-1
# Use periodogram
mba<-F
estim_MBA<-T
L_sym<-1000
# Length of long data (for computing the target)
len1<-3000
# Length of in-sample span
len<-120
# Frequency grid: length of periodogram i.e. len/2
K<-len/2
# Periodicity
periodicity<-12
cutoff<-pi/periodicity
# Specify filter length
L<-2*periodicity
# Real-time design
Lag<-0
# no constraints
i1<-i2<-F
# Use original (not differenced) data
dif<-F


# Proceed to estimation
for_sim_obj<-for_sim_out(a_vec,len1,len,cutoff,L,mba,estim_MBA,L_sym,
                         Lag,i1,i2,scaled_ATS,lambda_vec,eta_vec,anzsim,K,dif)



# ATS-components
#   -Decomposition of MSE into Accuracy (A), Timeliness (T) and Smoothness (S), see trilemma paper Wildi/Mcelroy
#   -Customization here emphasizes T (lambda>0) and S (eta>0) simultaneously at cost of A
#   -As can be seen in table below: A (of customized design) degrades massiveley in favour of substantial improvements by T and S
#   -MSE (of customized) degrades obviously, since customization is a deliberate departure of MSE-performances
#     A degrades overproportionally when improving S/T
ats_sym_ST<-for_sim_obj$ats_sym
# 2 Curvature, Peak Correlation, ...
amp_shift_mat_sim<-for_sim_obj$amp_shift_mat_sim
# 3. Amplitude and time-shifts
amp_sim_per<-for_sim_obj$amp_sim_per
shift_sim_per<-for_sim_obj$shift_sim_per
# 4. Output series
xff_sim<-for_sim_obj$xff_sim
# 5. Peak correlation and Curvature
amp_shift_mat_sim_ST<-for_sim_obj$amp_shift_mat_sim
dim_names<-for_sim_obj$dim_names
# Process 1 (positive acf) or 2 (zero acf) 
DGP<-1
ats_sym_ST[-1,,DGP,1]

# Amplitude and time shift
par(mfrow=c(1,2))
mplot<-amp_sim_per[,-1,DGP,1]
dimnames(mplot)[[2]]<-paste("Amplitude (",lambda_vec,",",eta_vec,")",sep="")
ax<-rep(NA,ncol(mplot))
ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
plot_title<-paste("Amplitude functions: a1=",a_vec[DGP],sep="")
insamp<-1.e+90
title_more<-dimnames(mplot)[[2]]
colo<-rainbow(ncol(mplot))
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
mplot<-shift_sim_per[,-1,DGP,1]
dimnames(mplot)[[2]]<-paste("Time-shifts (",lambda_vec,",",eta_vec,")",sep="")
ax<-rep(NA,ncol(mplot))
ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
plot_title<-paste("Time-shifts: a1=",a_vec[DGP],sep="")
insamp<-1.e+90
title_more<-dimnames(mplot)[[2]]
colo<-rainbow(ncol(mplot))
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)

par(mfrow=c(1,1))
series_vec<-c(2,3)
xf_per<-xff_sim[940:(940+len),,DGP,1]#dim(xff_sim)
dimnames(xf_per)[[2]]<-dim_names[[1]]
anf<-1
enf<-len
mplot<-scale(cbind(xf_per[,series_vec[1]],xf_per[,series_vec[2]])[anf:enf,])  #head(xf_per)
plot(as.ts(mplot[,1]),type="l",axes=F,col="red",ylim=c(min(na.exclude(mplot)),
                                                         max(na.exclude(mplot))),ylab="",xlab="",
main=paste("MSE (red) vs. Customized (cyan): a1=",a_vec[DGP],sep=""),lwd=1)
mtext("MSE", side = 3, line = -1,at=(enf-anf)/2,col=colo[1])
i<-2
lines(as.ts(mplot[,i]),col=colo[2],lwd=2)
mtext(paste("Customized: ",dimnames(xf_per)[[2]][series_vec[2]],sep=""), side = 3, line = -i,at=(enf-anf)/2,col=colo[2])
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*(enf-anf)/6),
     labels=as.integer(anf+(0:6)*(enf-anf)/6))
axis(2)
box()




#-----------------------------------------------------------------------------------------------
# Example 5: simulation exercise
#   -The ATS components are 'new' and therefore it's not immediately clear what they mean
#     -What does smaller S and T mean or: how has the practitioner to interpret these measures? 
#     -Explanations will be provided in new book but we here instead rely on well-known/established alternative statistics
#   -Alternative statistics: 
#     1. Instead of T we propose to look at peak-correlation
#       -Shift filter outputs of a specific one-sided design (for example classic MSE) with respect to target 
#         (output of symmetric filter) until the correlation between both series is maximized
#       -A smaller shift implies that the corresponding design is faster (smaller lag)
#       -Note that this concept (peak correlation) is scale invariant (scaling the filter output does not affect the measure)
#       -Expectation: emphasizing T (lambda>0) will result in faster filters (smaller shift at peak correlation)
#     2. Instead of S we propose to look at (relative) curvature
#       -'Smoothness' of a series (here: filter output) can be measured by looking at the (squared) second order differences
#       -If the noise-leakage is strong (poor stopband properties of the filter) then the filter-output will be 
#         noisy and squared second-order differences will be large
#       -In contrast: if the leakage is weak (strong suppression of noise) then the squared second order differences will be small
#       -The (relative) curvature is defined as follows: mean-squared second-order diffs divided by variance of series
#       -Note that this concept (relative curvature) is scale invariant (scaling the filter output does not affect the measure)
#       -Expectation: emphasizing S (eta>0) will result in smaller curvature

#  -Experimental design: in the following empirical experiment we
#    -compute data: multiple realizations of three differenet processes with positive/zero/negative autocorrelation 
#    -compute 
#      1. Symmetric target filters (in order to calculate peak correlation and (time-domain) MSEs)
#          These filter look into the future: therefore we expect that one-sided filters will be lagging (positive shift at peak correlation)
#      2. One-sided best possible MSE (assuming knowledge of the true model: no estimation): these filters are benchmarks
#        Ideally we would like empirical customized filters to outperform this benchmark in terms of lag/curvature out-of-sample....
#      3. Empirical (DFA) MSE and customized designs based on the periodogram (we assume 120 observations for the in-sample span)
#        -We have 3 customized designs
#          a. emphasize mainly T (specialized fast design) 
#          b. emphasize mainly S (specialized smooth design)  
#          c. 'best compromise' of T&S
#        -Our hope is: c. will outperform best MSE (assuming knowledge of true model) in terms of speed/curvature out-of-sample
#          for all processes considered
#    -Compute for each filter-realization of each process 
#      1. the peak-correlation (shift/lag)
#      2. the curvature
#      3. the (time-domain) MSE
#      In-sample and out-of-sample
#    -Plot the empirical distributions (box-plots) of all three measures: lag/curvature/MSE, in-sample and out-of-sample

#   -Expectations: 
#     1. MSE designs    
#      -Benchmark (best MSE assuming knowledge of true model) will outperform all other designs in terms of MSE out-of-sample
#       -DFA-MSE (based on periodogram) will outperform all customized designs in terms of MSE out-of-sample
#       -DFA-MSE (based on periodogram) will outperform benchmark in terms of MSE in-sample; but it will losse out-of-sample (overfitting)
#     2. Customized designs
#       -Emphasizing mainly T (fourth design below) will outperform all other filters in terms of 'peak correlation' (smallest delay with respect to target: ideally zero-shift)
#         But... this specialized design will be outperformed by classic MSE-design in terms of S (stronger leakage, noisy)
#       -Emphasizing mainly S (second design below) will outperform all other filters in terms of 'curvature' (smoothest output, strongest noise suppression) 
#         But... this specialized design will be outperformed by classic MSE-design in terms of T (larger lag)
#       -Best mix of S&T will outperform classic MSE in terms of peak-correlation (faster) AND curvature (stronger noise suppression)
#         This double score is not possible in a classic MSE-perspective
#         The ATS-trilemma allow to improve both S and T (peak-cor and curvature) at the expense of A (and MSE)
# Let's start

# Number of realizations for computing empirical distributions
anzsim<-100
# Specify the processes: ar(1) with coefficients -0.9,0.1 and 0.9
#   One can specify different processes at once and the code will loop through 
#   For simplicity we here use positive and nearly zero autocorrelation
a_vec<-c(0.9,0.08,-0.9)
# Specify the lambdas: first filter is MSE, second is specialized S (very smooth/high lag), 
#   third is 'best mix' (outperforms MSE in terms of S and T or in smaller lag and smaller curvature), 
#   last one is specialized T (smallest/vanishing lag but noisy) 
lambda_vec<-c(0,0,30,500)
# Specify the etas
eta_vec<-c(0,1.5,1,0.3)
# Ordinary ATS-components
scaled_ATS<-F
# Use periodogram
mba<-F
estim_MBA<-T
L_sym<-1000
# Length of long data (for computing the target)
len1<-3000
# Length of in-sample span
len<-120
# Frequency grid: length of periodogram i.e. len/2
K<-len/2
# Periodicity
periodicity<-12
cutoff<-pi/periodicity
# Specify filter length
L<-2*periodicity
# Real-time design
Lag<-0
# no constraints
i1<-i2<-F
# Use original (not differenced) data
dif<-F





# Proceed to simulation run
for_sim_obj<-for_sim_out(a_vec,len1,len,cutoff,L,mba,estim_MBA,L_sym,Lag,
                         i1,i2,scaled_ATS,lambda_vec,eta_vec,anzsim,M,dif)


# Extract sample performances
amp_shift_mat_sim<-for_sim_obj$amp_shift_mat_sim
amp_sim_per<-for_sim_obj$amp_sim_per
shift_sim_per<-for_sim_obj$shift_sim_per
xff_sim<-for_sim_obj$xff_sim
xff_sim_sym<-for_sim_obj$xff_sim_sym
ats_sym<-for_sim_obj$ats_sym
dim_names<-for_sim_obj$dim_names

dim(amp_shift_mat_sim)

# Plot: empirical distributions of MSEs, peak-correlation and curvature, in-sample and out-of-sample, for all processes specified in a_vec
colo<-c("red","orange","yellow","green","blue")#rainbow(length(lambda_vec)+1)

Perf_meas_sel<-c(3,7,4,8,5,9,6,10)
for (DGP in 1:length(a_vec))#DGP<-2
{
  par(mfrow=c(2,2))
  for (Perf_meas in Perf_meas_sel[1:4])
  {
    boxplot(list(amp_shift_mat_sim[1,Perf_meas,DGP,],amp_shift_mat_sim[2,Perf_meas,DGP,], amp_shift_mat_sim[3,Perf_meas,DGP,],amp_shift_mat_sim[4,Perf_meas,DGP,],amp_shift_mat_sim[5,Perf_meas,DGP,]),outline=T,names=c("Best MSE",paste("(",lambda_vec,",",eta_vec,")",sep="")),main=paste(dim_names[[2]][Perf_meas],", a1=",a_vec[DGP],sep=""),cex.axis=0.8,col=colo)
  }
  par(mfrow=c(1,2))
  for (Perf_meas in Perf_meas_sel[5:6])
  {
    boxplot(list(amp_shift_mat_sim[1,Perf_meas,DGP,],amp_shift_mat_sim[2,Perf_meas,DGP,],
                 amp_shift_mat_sim[3,Perf_meas,DGP,],amp_shift_mat_sim[4,Perf_meas,DGP,],amp_shift_mat_sim[5,Perf_meas,DGP,]),outline=T,
            names=c("Best MSE",paste("(",lambda_vec,",",eta_vec,")",sep="")),
            main=paste(dim_names[[2]][Perf_meas],", a1=",a_vec[DGP],sep=""),cex.axis=0.8,col=colo,notch=F)
  }
}



# Comparison: MSE vs. customized filter outputs
par(mfrow=c(1,1))
amp_shift_mat_sim<-for_sim_obj$amp_shift_mat_sim
amp_sim_per<-for_sim_obj$amp_sim_per
shift_sim_per<-for_sim_obj$shift_sim_per
xff_sim<-for_sim_obj$xff_sim
xff_sim_sym<-for_sim_obj$xff_sim_sym
ats_sym<-for_sim_obj$ats_sym
dim_names<-for_sim_obj$dim_names
DGP<-2
xf_per<-xff_sim[940:(940+2*len),,DGP,10]
dimnames(xf_per)[[2]]<-dim_names[[1]]
anf<-1
enf<-2*len
mplot<-scale(cbind(xf_per[,1],xf_per[,4])[anf:enf,])  #head(xf_per)
plot(as.ts(mplot[,1]),type="l",axes=F,col="red",ylim=c(min(na.exclude(mplot)),
                                                         max(na.exclude(mplot))),ylab="",xlab="",
main=paste("Benchmark MSE (red) vs. Customized balanced (green)",sep=""),lwd=2)
mtext("in sample",side = 3, line = -1,at=60,col="black")
mtext("out-of-sample",side = 3, line = -1,at=180,col="black")
mtext("Benchmark MSE", side = 3, line = -1,at=(enf-anf)/2,col="red")
i<-2
lines(as.ts(mplot[,i]),col=colo[4],lwd=2)
mtext(paste("Customized: ",dimnames(xf_per)[[2]][4],sep=""), side = 3, line = -i,at=(enf-anf)/2,col=colo[4])
abline(v=120)
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*(enf-anf)/6),
     labels=as.integer(anf+(0:6)*(enf-anf)/6))
axis(2)
box()




#---------------------------------------------------------------------------------------------------------
# Example 6: compare bivariate leading indicator and univariate customized
#   See previous tutorial for a background to the play_bivariate_func function below
#     This function sets-up a MDFA-experiment based on a bivariate noisy leading-indicator design

a1<-0.1
# Customization settings DFA
lambda_vec<-c(0,30)
eta_vec<-c(0,1)
# target
cutoff<-pi/12
len1<-2000
len<-120
L<-24
Lag<-0
i1<-i2<-F
# MDFA: MSE design
lambda_mdfa<-eta_mdfa<-0
troikaner<-F

# Run the competition: the new function handles the multivariate case
cust_leading_obj<-mdfa_mse_leading_indicator_vs_dfa_customized(anzsim,
                                                               a1,cutoff,L,lambda_vec,eta_vec,len1,len,i1,i2,Lag,
                                                               lambda_mdfa,eta_mdfa,troikaner)  
par(mfrow=c(1,2))
boxplot(list(cust_leading_obj$perf_in_sample[,1,1],cust_leading_obj$perf_in_sample[,1,2],cust_leading_obj$perf_in_sample[,1,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("Curvature in-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))
boxplot(list(cust_leading_obj$perf_out_sample[,1,1],cust_leading_obj$perf_out_sample[,1,2],cust_leading_obj$perf_out_sample[,1,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("Curvature out-of-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))

par(mfrow=c(1,2))
boxplot(list(cust_leading_obj$perf_in_sample[,2,1],cust_leading_obj$perf_in_sample[,2,2],cust_leading_obj$perf_in_sample[,2,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("Peak-Correlation in-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))
boxplot(list(cust_leading_obj$perf_out_sample[,2,1],cust_leading_obj$perf_out_sample[,2,2],cust_leading_obj$perf_out_sample[,2,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("Peak-Correlation out-of-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))

par(mfrow=c(1,2))
boxplot(list(cust_leading_obj$perf_in_sample[,3,1],cust_leading_obj$perf_in_sample[,3,2],cust_leading_obj$perf_in_sample[,3,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("MSE in-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))
boxplot(list(cust_leading_obj$perf_out_sample[,3,1],cust_leading_obj$perf_out_sample[,3,2],cust_leading_obj$perf_out_sample[,3,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("MSE out-of-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))


par(mfrow=c(1,1))

mplot<-scale(cust_leading_obj$filter_output_in_sample) 
dimnames(mplot)[[2]]<-dimnames(cust_leading_obj$filter_output_in_sample)[[2]]
colo_cust<-c("orange","green","brown")
plot(as.ts(mplot[,1]),type="l",axes=F,col=colo_cust[1],ylim=c(min(na.exclude(mplot)),max(na.exclude(mplot))),ylab="",xlab="",main=paste("Filter outputs: last realization",sep=""),lwd=1)
mtext(dimnames(mplot)[[2]][1], side = 3, line = -1,at=nrow(mplot)/2,col=colo_cust[1])
for (i in 2:(ncol(mplot)-1))
{
  lines(mplot[,i],col=colo_cust[i],lwd=1)
  mtext(dimnames(mplot)[[2]][i], side = 3, line = -i,at=nrow(mplot)/2,col=colo_cust[i])
}
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*nrow(mplot)/6),
     labels=c(1,rep(0,6))+as.integer((0:6)*nrow(mplot)/6))
axis(2)
box()


















