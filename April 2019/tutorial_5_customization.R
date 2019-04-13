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
#   It computes much more statistics and it can be used for simulations
#   It relies on DFA-code (which is the same as MDFA-univariate)

eta_vec<-c(0,1.8)
# Specify the fixed lambda
lambda_vec<-c(0,128)

# Specify the processes: ar(1) with coefficients -0.9,0.1 and 0.9
a_vec<-0.9
# Ordinary ATS-components
scaled_ATS<-F
# Generate a single realization of the processes
anzsim<-1
# Specify filter length
L<-24
# Use periodogram
mba<-F
estim_MBA<-T
M<-len/2
L_sym<-1000
# Length of long data (for computing the target)
len1<-3000
# Length of estimation sample
len<-120
# cutoff
cutoff<-pi/12
# Real-time design
Lag<-0
# no constraints
i1<-i2<-F
# difference data
dif<-F


# Proceed to estimation
for_sim_obj<-for_sim_out(a_vec,len1,len,cutoff,L,mba,estim_MBA,L_sym,
                         Lag,i1,i2,scaled_ATS,lambda_vec,eta_vec,anzsim,M,dif)



# 1 ATS
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
i_process<-1

ats_sym_ST[-1,,i_process,1]

DGP<-1
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
ki<-1
xf_per<-xff_sim[940:(940+len),,ki,1]#dim(xff_sim)
dimnames(xf_per)[[2]]<-dim_names[[1]]
anf<-1
enf<-len
mplot<-scale(cbind(xf_per[,series_vec[1]],xf_per[,series_vec[2]])[anf:enf,])  #head(xf_per)
plot(as.ts(mplot[,1]),type="l",axes=F,col="red",ylim=c(min(na.exclude(mplot)),
                                                         max(na.exclude(mplot))),ylab="",xlab="",
main=paste("MSE (red) vs. Customized (cyan): a1=",a_vec[ki],sep=""),lwd=1)
mtext("MSE", side = 3, line = -1,at=(enf-anf)/2,col=colo[1])
i<-2
lines(as.ts(mplot[,i]),col=colo[2],lwd=2)
mtext(paste("Customized: ",dimnames(xf_per)[[2]][series_vec[2]],sep=""), side = 3, line = -i,at=(enf-anf)/2,col=colo[2])
axis(1,at=c(1,rep(0,6))+as.integer((0:6)*(enf-anf)/6),
     labels=as.integer(anf+(0:6)*(enf-anf)/6))
axis(2)
box()




#-----------------------------------------------------------------------------------------------
# Example 6: simulation exercise
#   Empirical distributions of MSEs, peak-correlation and cusrvature, in-sample and out-of-sample, for all processes specified in a_vec
anzsim<-100
# Specify the processes: ar(1) with coefficients -0.9,0.1 and 0.9
a_vec<-c(0.9,0.1,-0.9)
# Ordinary ATS-components
scaled_ATS<-F
# Specify the lambdas
lambda_vec<-c(0,0,30,500)
# Specify the etas
eta_vec<-c(0,1.5,1,0.3)
# Specify filter length
L<-24
# Use periodogram
mba<-F
estim_MBA<-T
M<-len/2
# Length of symmetric target filter (for computing MSEs)
L_sym<-2*939
# Length of long data
len1<-2000
# Length of estimation sample
len<-120
# cutoff
cutoff<-pi/12
# Real-time design
Lag<-0
# no constraints
i1<-i2<-F
# difference data
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



# Plot: empirical distributions of MSEs, peak-correlation and cusrvature, in-sample and out-of-sample, for all processes specified in a_vec
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
ki<-2
xf_per<-xff_sim[940:(940+2*len),,ki,10]
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
# Example 5: compare bivariate leading indicator and univariate customized
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


















