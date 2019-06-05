



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

insamp<-1200

period_high<-6

sigma_low<-2
len<-1200
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
high_freq_diff<-T
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
target_as_explanatory<-T
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
# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high

weight_func_embed[,1]<-exp(-lead*1.i*(0:K)*pi/K)*weight_func_embed_h[,1]

# Specify cutoff of ideal lowpass
periodicity<-6
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
mean((y-yhat_mixed)^2,na.rm=T)
if (!target_as_explanatory)
{
  print("The next MSE is cheating in this configuration because it uses low-freq data as estimate")
}  
mean((y-data_mat[,1])^2,na.rm=T)
mdfa_obj_mixed$MS_error
mean(diff(diff(yhat_mixed[L:length(yhat_mixed)]))^2)


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Marc's 1. idea
target_as_explanatory<-T

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

# Lead time: fractional 1/period_high corresponds to 1 time-unit on high-frequency scale
lead<-0/period_high
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
mean((y-yhat_agg)^2,na.rm=T)
if (!target_as_explanatory)
{
  print("The next MSE is cheating in this configuration because it uses low-freq data as estimate")
}  
mean((y-data_mat[,1])^2,na.rm=T)
mdfa_obj$MS_error


mean(diff(diff(yhat_agg[L:length(yhat_agg)]))^2)

