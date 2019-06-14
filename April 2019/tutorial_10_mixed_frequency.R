



# Purpose of tutorial: 

rm(list=ls())

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA",force=TRUE)
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
source("Common functions/mdfa_trade_func.r")
source("Common functions/ideal_filter.r")
source("Common functions/MSE_perf.r")
source("Common functions/mixed_freq_simulation_embed_vs_fold.r")

# Specify experimental setting
# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100





set.seed(10)
x_low<-rep(NA,len)
eps_low<-sigma_low*arima.sim(n=len,list(ar = c(0.09), ma = c(0)))

eps_high<-arima.sim(n=period_high*len,list(ar = c(0.09), ma = c(0)))

# Low-frequency data is flow-data i.e. differences of high-frequency data sum up to build diff of low-frequency data
#   Then add an idiosyncratic error term eps_low
for (i in 1:len)
  x_low[i]<-sum(eps_high[(i-1)*period_high+1:period_high])+eps_low[i]

# Tucker's Idea
#   Embedding of high-frequency data in low-frequency sampling-scheme (period_high columns)
x_high_embed<-matrix(ncol=period_high,nrow=len)   
index_mat<-NULL
# Use high-freq differences (for example monthly) or low-freq diff (for example quarterly)
for (i in 1:len)
{
  for (j in 1:period_high)
  {
    if (!high_freq_diff)
    {
# 1. low-freq diff
      index<-(i-2)*period_high+(1:period_high)+j
    } else
    {
# 2. high-freq diff
      index<-(i-1)*period_high+j
    }
    if (min(index)<1)
      index<-rep(1,period_high)
    x_high_embed[i,j]<-sum(eps_high[index])
    index_mat<-cbind(index_mat,as.vector(index))
  }
}
tail(t(index_mat)) 
cor(x_high_embed[,c(1,3)])

data_mat<-cbind(x_low,x_high_embed)
colnames(data_mat)<-c("Low",paste("High lag ",(period_high-1):0,sep=""))
tail(data_mat)#dim(data_mat)
ts.plot(data_mat)

#acf(data_mat)

#Use target series a sexplanatory variable: not suitable for GDP since GDP releases are delayed and noisy (subject to revisions)
if (!target_as_explanatory)
{
  weight_func_embed_h<-spec_comp(nrow(data_mat[1:insamp,]), cbind(data_mat[1:insamp,1],data_mat[1:insamp,2:ncol(data_mat)]), 0)$weight_func
  colnames(weight_func_embed_h)<-c("target",colnames(data_mat)[2:ncol(data_mat)])
} else
{
  weight_func_embed_h<-spec_comp(nrow(data_mat[1:insamp,]), cbind(data_mat[1:insamp,1],data_mat[1:insamp,]), 0)$weight_func
  colnames(weight_func_embed_h)<-c("target",colnames(data_mat))
}
weight_func_embed<-weight_func_embed_h#dim(weight_func_embed)
# Resolution of frequency-grid
K<-nrow(weight_func_embed)-1

weight_func_embed[,1]<-exp(-lead*1.i*(0:K)*pi/K)*weight_func_embed_h[,1]

cutoff<-pi/periodicity
# Target (in frequency domain)
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length: L/K should be 'small'
L<-2*periodicity
# Estimate filter coefficients

if (sigma_low==0)
{
  print("System is singular: sigma_low must be larger than zero")
} else
{
  mdfa_obj_mixed<-MDFA_mse(L,weight_func_embed,Lag,Gamma)$mdfa_obj 
}


ts.plot(abs(mdfa_obj_mixed$trffkt))


b_mixed<-mdfa_obj_mixed$b
colnames(b_mixed)<-colnames(weight_func_embed)[2:ncol(weight_func_embed)]
rownames(b_mixed)<-paste("Lag ",0:(L-1),sep="")#dim(b_mat)
b_mixed

# Filtering: one-sided filter
yhat_mixed<-filt_func(data_mat[,colnames(b_mixed)],b_mixed)$yhat

# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
y<-ideal_filter_func(periodicity,M,data_mat[,1])$y

ts.plot(cbind(yhat_mixed,y),col=c("blue","red"))

# MSE performances

MSE_embed<-MSE_perf_func(insamp,y,yhat_mixed,len,mdfa_obj_mixed)$perf_mat
  


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Marc's 1. idea

weight_func_highh<-spec_comp(length(eps_high[1:(insamp*period_high)]),as.matrix(cbind(eps_high[1:(insamp*period_high)],eps_high[1:(insamp*period_high)])), 0)$weight_func[,1]
weight_func_high<-weight_func_highh#/sqrt(period_high)
weight_func_lowhh<-spec_comp(length(x_low[1:(insamp)]),as.matrix(cbind(x_low[1:insamp],x_low[1:insamp])), 0)$weight_func[,1]
weight_func_lowh<-c(weight_func_lowhh,rep(0.0,length(weight_func_high)-length(weight_func_lowhh)))
# Rescale low-freq data to high-freq scale 
weight_func_low<-weight_func_lowh*sqrt(period_high)
weight_func<-weight_funch<-cbind(weight_func_low,weight_func_low,weight_func_high)
if (!target_as_explanatory)
{
  weight_func<-weight_funch<-cbind(weight_func_low,weight_func_high)
} else
{
  weight_func<-weight_funch<-cbind(weight_func_low,weight_func_low,weight_func_high)
}

# abs(DFT) of target is larger because target is sum of high-freq-data + noise
ts.plot(abs(weight_func),col=c("blue","red"))

K<-nrow(weight_func)-1

weight_func[,1]<-exp(-lead*1.i*(0:K)*pi/K)*weight_funch[,1]#dim(weight_func)


# Target (in frequency domain)
#   -For low-frequency target the cutoff has to be adjusted with respect to high-freq scale
#   -This must correspond to periodicity selected when computing the ideal filter below (otherwise targets differ)
Gamma<-(0:(K))<=K*cutoff/(period_high*pi)+1.e-9#pi/18
ts.plot(Gamma)

# Filter length high-freq: account for narrowing high-freq time-scale
L<-2*period_high*periodicity

lin_eta <- F
lambda <- 0
eta <- 0
weight_constraint <- rep(1/(ncol(weight_func) - 1), ncol(weight_func) - 
                           1)
lambda_cross <- lambda_smooth <- 0
lambda_decay <- c(0, 0)
lin_expweight <- F
shift_constraint <- rep(0, ncol(weight_func) - 1)
grand_mean <- F
b0_H0 <- NULL
c_eta <- F
weights_only <- F
weight_structure <- c(0, 0)
white_noise <- F
synchronicity <- F
del_high<-1
# Lag matrix of low-frequency series must be expanded according to period_high 
if (!target_as_explanatory)
{
  lag_mat <- matrix(cbind(period_high*(0:(L - 1)),del_high*(0:(L - 1))), nrow = L)
} else
{
  lag_mat <- matrix(cbind(period_high*(0:(L - 1)),period_high*(0:(L - 1)),del_high*(0:(L - 1))), nrow = L)
}
troikaner <- F
i1 <- i2 <- F


mdfa_obj <- mdfa_analytic(L, lambda, weight_func, Lag, Gamma, 
                          eta, cutoff, i1, i2, weight_constraint, lambda_cross, 
                          lambda_decay, lambda_smooth, lin_eta, shift_constraint, 
                          grand_mean, b0_H0, c_eta, weight_structure, white_noise, 
                          synchronicity, lag_mat, troikaner)

ts.plot(abs(mdfa_obj$trffkt),col=rainbow(ncol(mdfa_obj$trffkt)))
b_mixed<-mdfa_obj$b
colnames(b_mixed)<-colnames(weight_func)[2:ncol(weight_func)]
rownames(b_mixed)<-paste("Lag ",0:(L-1),sep="")#dim(b_mat)
b_mixed

# Filtering: one-sided filter

if (!target_as_explanatory)
{
  b_high<-b_mixed[,1:ncol(b_mixed)]
  yhat_highh<-filt_func(as.matrix(eps_high),as.matrix(b_high))$yhat
  
  # Extract the monthly values at the right month in the quarter
  yhat_high<-yhat_highh[round(-lead*period_high)+(length(yhat_highh)-as.integer(length(yhat_highh)/period_high)*period_high)+period_high*(1:(as.integer(length(yhat_highh)/period_high)))]
  yhat_agg<-yhat_high
  
} else
{
  yhat_low<-filt_func(as.matrix(x_low),as.matrix(b_mixed[,1]))$yhat
  b_high<-NULL
  for (i in 1:length(b_mixed[,2]))
    b_high<-c(b_high,c(b_mixed[i,2],rep(0,del_high-1)))
  
  yhat_highh<-filt_func(as.matrix(eps_high),as.matrix(b_high))$yhat
  
  # Extract the monthly values at the right month in the quarter
  yhat_high<-yhat_highh[-lead*period_high+(length(yhat_highh)-as.integer(length(yhat_highh)/period_high)*period_high)+period_high*(1:(as.integer(length(yhat_highh)/period_high)))]
  
  
  yhat_agg<-apply(cbind(yhat_low,yhat_high),1,sum)
}
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
y<-ideal_filter_func(periodicity,M,x_low)$y

ts.plot(cbind(yhat_agg,y),col=c("blue","red"))

# MSE performances
MSE_fold<-MSE_perf_func(insamp,y,yhat_agg,len,mdfa_obj)$perf_mat

perf_mat<-rbind(MSE_embed,MSE_fold)
rownames(perf_mat)<-c("Embed","Fold")
perf_mat


#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 1. Same settings as above

anzsim<-100
# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0

# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 2. Same as 1 but target_as_explanatory==T

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-T
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 3. Same as 1 but period_high<-3 (quarterly monthly) and L<-4*

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-3
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0

# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100

# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))

#-------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 4. Same as 3 but L<-4*periodicity

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-3
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Filter length: L too large leads to overfitting
L<-4*periodicity
# No regularization
lambda_cross<-0

# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100

# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 5. Same as 1 but sigma_low<-1

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-1
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 6. Same as 5 but periodicity<-2 
# We examine three cases: 6.1, 6.2 and 6.3
# 6.1: L<-2*periodicity
#   In this case L=2*periodicity might be a limiting factor i.e. too few freely determined coefficients
#   Selecting larger L (4*periodicity or 6*periodicity) will improve out-of-sample perfs (and folding then outperforms embedding) 

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-1
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-2
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


# Here for the first time performances of fold are marginally worse than embed out-of-sample

#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 6.2 Same as 6.1 but L<-6*periodicity 
#   Selecting larger L (4*periodicity or 6*periodicity) will improve out-of-sample perfs (and folding then outperforms embedding) 

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-1
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-2
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-6*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))

#--------------------------------------------------------------------------------------------------------
# 6.3: L<-4*periodicity or L<-8*periodicity (both were computed), insamp<-2400, len<-3000
#   Because L=2*periodicity might be a limiting factor we here set L<-4*periodicity
#   In order to avoid overfitting we select a very large in-sample span
#     So we can better compare asymptotic performances of embedding vs. folding
#     Embeding still marginally better than folding (but t-values are around 0.02)

# Sample length
len<-3000
# In-sample length
insamp<-2400
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-1
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-T
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-2
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
#   Both settings were computed/stored
L<-4*periodicity
L<-8*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# 6.4: as 6.1 but target_as_explanatory<-T
#   Now folding outperforms embedding in-sample (MSE and criterion value) and out-of-sample

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-1
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-T
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-2
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# 6.5: as 6.1 but periodicity<-1
#   Here we can see the problem of folding!!!!!!!!!!!!!

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-1
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-1
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# 6.6: as 6.5 but sigma_low<-0.0001

# Regularization: absolute cross-sectional




#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 7. Same as 6 but target_as_explanatory<-T 

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-1
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-T
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-2
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))

# Here fold is again marginally (not significantly) better out-of-sample

#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 8. Same as 1 but with autocorrelated data ar_low<-0.4, ar_high<-0.8 (the latter ought to be larger because it's on high-freq scale)

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.4
# DGP of differenced data
ar_high<-0.8
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))

#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 9. Same as 8 but L<-4*periodicity

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.4
# DGP of differenced data
ar_high<-0.8
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-4*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 10. Same settings as 1 but short insamp<-120


# Sample length
len<-1200
# In-sample length
insamp<-120
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 11. Same settings as 10 but period_high<-3 and periodicity<-2: this implies smaller L (less overfitting)


# Sample length
len<-1200
# In-sample length
insamp<-120
# Folding rate
period_high<-3
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-2
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))

#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 12. Same as 1 but high_freq_diff<-F

# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-F
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# No regularization
lambda_cross<-0


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#--------------------------------------------------------------------------------------------------------
# Simulation experiment: compare embedding (Tucker) and folding (Marc)
# 13. Same as 1. but lambda_cross<-1 
#   -Model: low-freq data is stock data i.e. diffs of low-freq data are sum of high-freq diffs within low-freq span
#   -If this model assumption is true/pertinent, then all embedded high-freq series must be treated equally
#     i.e. impose full cross-sectional regularization
#   -The empirical results below confirm pertinence: embedding with regularization performs as well as folding 
#     and the former outperforms embedding (without regularization) in 6.1 out-of-sample

anzsim<-100
# Sample length
len<-1200
# In-sample length
insamp<-600
# Folding rate
period_high<-6
# strength of low-freq indiosyncratic component
sigma_low<-3
# Use differences of high-freq data on high-freq scale (for example monthly) or on low-freq scale (for example quarterly)
high_freq_diff<-T
# Use low-freq target data as explanatory too (not suitable when target is GDP: publication-lag and revisions)
target_as_explanatory<-F
# DGP of differenced idiosyncratic component (low-freq data)
#   Low-freq data is flow-data (sum of high-freq: for example sum in months of quarter) + idiosyncratic (independent) component
ar_low<-0.09
# DGP of differenced data
ar_high<-0.09
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
# Target: specify cutoff of ideal lowpass
periodicity<-6
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Filter length: L too large leads to overfitting
L<-2*periodicity
# Full cross-sectional regularization
lambda_cross<-1


# Computations need approx 5-10 mins: Results were previously stored
perform_computations<-F

if (perform_computations)
{
  set.seed(1)
  anzsim<-100
  perf_mat<-matrix(rep(0,2*3),nrow=2,ncol=3)
  perf_array<-array(dim=c(anzsim,dim(perf_mat)))
  var_perf<-perf_mat
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  
  for (i in 1:anzsim)
  {
    perf_math<-simulation_embed_vs_fold_reg(len,sigma_low,ar_low,ar_high,period_high,high_freq_diff,target_as_explanatory,lead,periodicity,M,L,lambda_cross)$perf_mat
    perf_array[i,,]<-perf_math
    perf_mat<-perf_mat+perf_math
    setTxtProgressBar(pb, i)
  }  
  for (i in 1:nrow(var_perf))
  {
    for (j in 1:ncol(var_perf))
    {
      var_perf[i,j]<-var(perf_array[,i,j])
    }
  }
  perf_mat<-perf_mat/anzsim
  save(perf_mat,file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  save(var_perf,file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
} else
{
  load(file=paste("output/perf_mat_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
  load(file=paste("output/var_perf_",len,"_",insamp,"_",period_high,"_",sigma_low,"_",high_freq_diff,"_",target_as_explanatory,"_",ar_low,"_",ar_high,"_",round(lead,3),"_",periodicity,"_",L,"_",lambda_cross,sep=""))
}

perf_mat
# Test for significance of differences: larger than 2 (in abs) means significance
sqrt(anzsim)*(perf_mat[1,]-perf_mat[2,])/(2*sqrt(apply(var_perf,2,sum)))


#------
# Shorter in-sample
# One month ahead