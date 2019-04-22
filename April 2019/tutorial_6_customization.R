# Todos
#   Link to trilemma paper





# Purpose of tutorial: 
# -Illustrate decomposition of classic MSE-norm into Accuracy, Timeliness and Smoothness components: ATS-trilemma (see McElroy/Wildi)
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
#      -DFA-MSE (based on periodogram) will outperform all customized designs in terms of MSE out-of-sample
#      -DFA-MSE (based on periodogram) will outperform benchmark in terms of MSE in-sample; but it will loose out-of-sample (overfitting)
#     2. Customized designs
#       -Emphasizing mainly T (fourth design below) will outperform all other filters in terms of 'peak correlation' (smallest delay with respect to target: ideally zero-shift)
#         But... this specialized design will be outperformed by classic MSE-design in terms of S (stronger leakage, noisy)
#       -Emphasizing mainly S (second design below) will outperform all other filters in terms of 'curvature' (smoothest output, strongest noise suppression) 
#         But... this specialized design will be outperformed by classic MSE-design in terms of T (larger lag)
#       -Best mix of S&T will outperform best possible MSE in terms of peak-correlation (faster) AND curvature (stronger noise suppression)
#         This double score is not possible in a classic MSE-perspective
#         The ATS-trilemma allow to improve both S and T (peak-cor and curvature) at the expense of A (and MSE): a trilemma is needed...

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
# Length symmetric filter
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

# Proceed to simulation
for_sim_obj<-for_sim_out(a_vec,len1,len,cutoff,L,mba,estim_MBA,L_sym,Lag,
                         i1,i2,scaled_ATS,lambda_vec,eta_vec,anzsim,K,dif)

# Extract sample performances
amp_shift_mat_sim<-for_sim_obj$amp_shift_mat_sim
amp_sim_per<-for_sim_obj$amp_sim_per
shift_sim_per<-for_sim_obj$shift_sim_per
xff_sim<-for_sim_obj$xff_sim
xff_sim_sym<-for_sim_obj$xff_sim_sym
ats_sym<-for_sim_obj$ats_sym
dim_names<-for_sim_obj$dim_names



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

DGP<-2
par(mfrow=c(1,1))
amp_shift_mat_sim<-for_sim_obj$amp_shift_mat_sim
amp_sim_per<-for_sim_obj$amp_sim_per
shift_sim_per<-for_sim_obj$shift_sim_per
xff_sim<-for_sim_obj$xff_sim
xff_sim_sym<-for_sim_obj$xff_sim_sym
ats_sym<-for_sim_obj$ats_sym
dim_names<-for_sim_obj$dim_names
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
#   -See previous tutorial for a background to the play_bivariate_func function below
#     The bivariate design relies on a leading indicator
#     It outperformed the univariate DFA in terms of MSE in-sample (of course) and also out-of-sample (expected)
#   -Given the above outcome (of simulation experiment) we can conclude that bivariate design will also outperform 
#     all customized designs in terms of MSE (in- and out-of-sample)
#   -Question: 
#     Can the 'best-mix' customized design outperform the bivariate MSE (with leading indicator) in terms of lead/curvature?
# Experimental design
#   -The function mdfa_mse_leading_indicator_vs_dfa_customized sets-up a corresponding experiment 
#   -It compares 
#     1. a univariate DFA-MSE (lambda=eta=0)
#     2  a univariate 'best-mix' customized DFA (see above experiment)
#     3. a bivariate MSE design based which adds a noisy leading-indicator to the set of explanatory series

# Select data generating process (ar(1)-coefficient)
a1<-0.08
# Number replications
anzsim<-500
# Customization settings DFA: MSE and 'best-mix'
lambda_vec<-c(0,30)
eta_vec<-c(0,1)
# target
periodicity<-12
cutoff<-pi/periodicity
# Full sample lengt (for applying symmetric ideal lowpass)
len1<-3000
# In-sample span (for estimation of spectrum: dft)
len<-240
L<-2*periodicity
# Nowcast
Lag<-0
# No restrictions
i1<-i2<-F
# MDFA: MSE design
lambda_mdfa<-eta_mdfa<-0
# Boolean for speeding up simulation (some statistics in MDFA are omitted: does not impact calculations (only computation time))
troikaner<-F

# Run the competition: univariate MSE and customized vs. bivariate leading indicator
cust_leading_obj<-mdfa_mse_leading_indicator_vs_dfa_customized(anzsim,a1,cutoff,L,lambda_vec,eta_vec,len1,len,i1,i2,Lag,lambda_mdfa,eta_mdfa,troikaner)  

# The following comments assume a1<-0.08 (almost white noise i.e. log-returns of typical (positive) economic time series)
# 1. Curvature
#     -Best-mix customized (gren) outperforms bivariate (brown) out-of-sample (stronger noise suppression)
#     -Bivariate (brown) marginally better than univariate MSE (orange) out-of-sample 
par(mfrow=c(1,2))
boxplot(list(cust_leading_obj$perf_in_sample[,1,1],cust_leading_obj$perf_in_sample[,1,2],cust_leading_obj$perf_in_sample[,1,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("Curvature in-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))
boxplot(list(cust_leading_obj$perf_out_sample[,1,1],cust_leading_obj$perf_out_sample[,1,2],cust_leading_obj$perf_out_sample[,1,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("Curvature out-of-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))

# 2. Lag at peak-correlation
#     -Best-mix customized (gren) outperforms bivariate (brown) out-of-sample (lead by one time-point)
#       This outcome is both unexpected and remarkable
#     -Bivariate (brown) outperforms univariate MSE (orange) out-of-sample by one time-point 
#       Expected outcome (because of leading indicator in bivariate design)
par(mfrow=c(1,2))
boxplot(list(cust_leading_obj$perf_in_sample[,2,1],cust_leading_obj$perf_in_sample[,2,2],cust_leading_obj$perf_in_sample[,2,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("Peak-Correlation in-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))
boxplot(list(cust_leading_obj$perf_out_sample[,2,1],cust_leading_obj$perf_out_sample[,2,2],cust_leading_obj$perf_out_sample[,2,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("Peak-Correlation out-of-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))

# 3. MSE
#     No surprise i.e. everything as expected
boxplot(list(cust_leading_obj$perf_in_sample[,3,1],cust_leading_obj$perf_in_sample[,3,2],cust_leading_obj$perf_in_sample[,3,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("MSE in-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))
boxplot(list(cust_leading_obj$perf_out_sample[,3,1],cust_leading_obj$perf_out_sample[,3,2],cust_leading_obj$perf_out_sample[,3,3]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),"MDFA-MSE Leading Indicator"),main=paste("MSE out-of-sample, a1=",a1,sep=""),cex.axis=0.8,col=c("orange","green","brown"))


# Compare filter outputs
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

#---------------------------------------------------------------------------------------------------------
# Example 7: same as example 6 but we now allow for customization of the bivariate design
# Contenders in this competition: uni and bivariate MSE as well as uni and bivariate customized (the latter is new)


# Use the same settings as above but add a 'best-mix' customization for the bivariate filter
#   You might have to load the settings in the previous example in order to run this piece (if not done yet)
lambda_mdfa<-c(0,30)
eta_mdfa<-c(0,1.)


cust_leading_obj<-mdfa_mse_leading_indicator_vs_dfa_customized(anzsim,a1,
                                                               cutoff,L,lambda_vec,eta_vec,len1,len,i1,i2,Lag,lambda_mdfa,eta_mdfa,troikaner)  


colo<-rainbow(length(lambda_mdfa)+length(lambda_vec))
par(mfrow=c(1,2))
# 1. Curvature
#   in-sample
boxplot(list(cust_leading_obj$perf_in_sample[,1,1],cust_leading_obj$perf_in_sample[,1,2],cust_leading_obj$perf_in_sample[,1,3],cust_leading_obj$perf_in_sample[,1,4]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),paste("MDFA(",lambda_mdfa,",",eta_mdfa,")",sep="")),main=paste("Curvature in-sample, a1=",a1,sep=""),cex.axis=0.8,col=colo)
#   out-of-sample
boxplot(list(cust_leading_obj$perf_out_sample[,1,1],cust_leading_obj$perf_out_sample[,1,2],cust_leading_obj$perf_out_sample[,1,3],cust_leading_obj$perf_out_sample[,1,4]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),paste("MDFA(",lambda_mdfa,",",eta_mdfa,")",sep="")),main=paste("Curvature out-of-sample, a1=",a1,sep=""),cex.axis=0.8,col=colo)



# 2. Peak correlation 
#   in-sample
boxplot(list(cust_leading_obj$perf_in_sample[,2,1],cust_leading_obj$perf_in_sample[,2,2],cust_leading_obj$perf_in_sample[,2,3],cust_leading_obj$perf_in_sample[,2,4]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),paste("MDFA(",lambda_mdfa,",",eta_mdfa,")",sep="")),main=paste("Peak Correlation in-sample, a1=",a1,sep=""),cex.axis=0.8,col=colo)
#   out-of-sample 
boxplot(list(cust_leading_obj$perf_out_sample[,2,1],cust_leading_obj$perf_out_sample[,2,2],cust_leading_obj$perf_out_sample[,2,3],cust_leading_obj$perf_out_sample[,2,4]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),paste("MDFA(",lambda_mdfa,",",eta_mdfa,")",sep="")),main=paste("Peak Correlation out-of-sample, a1=",a1,sep=""),cex.axis=0.8,col=colo)


# 3. MSE 
#   in sample
boxplot(list(cust_leading_obj$perf_in_sample[,3,1],cust_leading_obj$perf_in_sample[,3,2],cust_leading_obj$perf_in_sample[,3,3],cust_leading_obj$perf_in_sample[,3,4]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),paste("MDFA(",lambda_mdfa,",",eta_mdfa,")",sep="")),main=paste("MSE in-sample, a1=",a1,sep=""),cex.axis=0.8,col=colo)
#   out-of-sample 
boxplot(list(cust_leading_obj$perf_out_sample[,3,1],cust_leading_obj$perf_out_sample[,3,2],cust_leading_obj$perf_out_sample[,3,3],cust_leading_obj$perf_out_sample[,3,4]),outline=T,names=c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),paste("MDFA(",lambda_mdfa,",",eta_mdfa,")",sep="")),main=paste("MSE out-of-sample, a1=",a1,sep=""),cex.axis=0.8,col=colo)

# Multivariate customized 
#   outperforms all other contenders with respect to lead and curvature out-of-sample (though outperformance with respect to customized DFA is modest)
#   outperforms customized univariate in terms of MSE out-of-sample
#   is outperformed in terms of MSE by both MSE-designs out-of-sample
# Conclusion: 
#   1. Performances in terms of Smoothess (smaller curvature) and Timeliness (smaller lag at peak correlation) are
#     obtained mainly by design, i.e. by the art of optimizing the relevant features of the filter, rather than by adding 
#     a (cheating...) leading time series.
#   2. MSE-performances may degrade substantially as a side-effect of addressing S and T at the expense of A
#     but some users (me included) really don't care about that collateral damage to MSE

#---------------------------------------------------------------------------------------------------------------------------------
# Wrap-up
# -The MSE-norm can be split into four error components, Accuracy, Timeliness, Smoothness and Residual, which are weighted equally in the original MSE norm
#   -S measures the noise suppression by the one-sided filter in the stopband
#   -T measures the shift (delay) of the one-sided filter in the passband
#   -A measures the level-tracking of the one-sided filter in the passband
#   -R measures the contribution of the shift in the stopband: in our applications this component vanishes invariably (because the target vanishes in the stopband)
#   -The remaining A, T and S account for the dilemma of optimizing amplitude (A and S) as well as phase (T) functions, recall tutorial 2 example 5
#     -One cannot improve both amplitude and phase functions fits simultaneously and arbitrarily well across the full frequency-band
# -Playing with the trilemma
#   -Emphasizing any (A-, T- or S-) error component (assigning unilaterally a larger weight to that component uniquely) inflates the other ones as well as the sum (i.e. MSE)
#   -Emphasizing either combination of two components (assigning larger weights to these two components) inflates overprortionally the remaining one as well as MSE
#   -Emphasizing S and T (at costs of A) improves simultaneously the amplitude (noise suppression in stop-band) as well as the phase (shift in the passband)
#      at cost of amplitude in passband (A or Accuracy-component: accounts for level-tracking ability of the one-sided filter) 
#   -Therefore amplitude and phase function fits can be improved simultaneously on parts of the frequency band: the stop-band (amplitude) and the passband (phase)   
#   -Classic econometric approaches are immanently incapable of tackling that problem because the classic maximum likelihood approach emphasizes a 
#     (one-step ahead) forecast problem; the forecasting target is an allpass filter i.e. there is no split of the frequency-band into pass- and stopbands
#     or, stated otherwise, S (smoothness) does not exist.
# -A comparison with classic econometric/time series approaches
#   -Classic one-step ahead mean-square forecasting decomposes the MSE-norm into A and T only (there is no S): classic econometric approaches emphasize a 
#     dilemma; they are by design incapable of addressing the proposed trilemma
#   -By tracking arbitrary targets (see tutorial 2) the DFA is more general than classic econometric approaches
#     1. DFA can replicate classic (one-steap ahead allpass) approaches, see tutorial 1
#     2. By allowing more generic targets (lowpass, bandpass) a trilemma can be spanned upon pass- and stop-bands of the target
#     3. The trilemma enables to address amplitude AND phase fits simultaneously in relevant frequency-bands
# -Interpretation of S and T
#   -Smoothness is intimately related to the classic curvature statistic (mean-square second order differences) which measures ... well... the curvature (i.e. smoothness) of the filter output
#   -Timeliness is intimately related to the classic peak-correlation concept (shift outputs of one-sided and of target filters until correlation is maximized)
#   -Since S and T can be improved simulatenously in the ATS-trilemma, at costs of A and overall MSE, we conclude that curvature and lag (at peak-correlation) 
#     can be improved simultenously, too.
#   -Our simulation studies above confirm this claim, in-sample as well as out-of-sample
#     -Curvature as well as lag at peak-correlation can be improved both substantially i.e. improvements are not marginal
#     -But A(ccuracy) and MSE-performances degrade: 
#       -This loss is mainly due to shrinkage of the amplitude function in the passband 
#       -The shrinkage could be remedied easily, at least to some extent, by re-scaling of the filter output. 
#       -Stated otherwise: part of the loss in A- and MSE-performances could be overcome by a very simple transformation which does not affect 
#         the scale-invariant (relative) curvature and lag at peak-correlation measures
# -Summary of empirical studies
#   -MSE designs outperform customized designs in- and out-of-sample, as expected
#   -The bivariate leading indicator MSE-design outperforms all other contenders in-sample and out-of-sample in terms of MSE (assuming overfitting is not excessively heavy)
#   -A univariate suitably customized DFA filter can outperform all MSE-designs in terms of curvature and (lag at) peak-correlation simultaneously in-sample and out-of-sample
#     1. the best possible univariate MSE-approach assuming knowledge of the true data-generating process is outperformed
#     2. more surprisingly, perhaps, the bivariate MSE leading indicator design is outperformed: this last result confirms that gains can be substantial in both dimensions at once

# -Cautionary words: avoid confusions
#   -All reported measures are aggregates of stochastic events
#   -Improved S and T or, equivalently, improved curvature and peak correlation imply that the corresponding filter output improves 'in the mean'
#     -At some turning-points the lead can be larger or smaller (than indicated by the aggregate peak-correlation number): there is variation (no determinism)
#     -Sometimes the output is less (or more) smooth than assumed by the aggregate curvature number
#   -Applying a customized filter to white noise can (and does) not improve our inferentail ability about the (completely random) future
#     -Neither forecasting nor any derived statistic (for example trading performances based on the sign of the filter output) can be improved
#     -In the mean, the customized filter will cross the zero line earlier than the MSE-filter: this will indeed be observed
#     -But this feature would be completely 'useless' in the context of an iid process (independent identically distributed): no utility could be derived for the user
#   -But all is not lost... in real-world markets liquidity is finite (adjustments are not immediate) and a substantial share of traders are relying on classic ('slow') MA-filters. 
#     -In such a context, improved timing by customized designs could deliver 
#     -However, experience suggests that faster is not always better
#     -In any case, the ATS-trilemma enables to trigger whatever is felt judicious




