rm(list=ls())

# Load libraries
library(mFilter)
library(Quandl)
library(tis)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA",force=T)
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)

# Source functions
source("Common functions/hpFilt.r")
source("Common functions/plot_func.r")



# One can refresh GDP from Quandl (a Quandl-API key is needed, though) 
download_data_from_Quandl<-F

if (download_data_from_Quandl)
{  
  start_year<-1947 
  start_date=paste(start_year,"-01-01",sep="")
# Last data point
  end_date<-format(Sys.time(), "%Y-%m-%d")
  end_year<-as.double(substr(end_date,1,4))
# Load Real GDP
#Title:               Real Gross Domestic Product, 3 Decimal
#Series ID:           GDPC96
#Source:              US. Bureau of Economic Analysis
#Release:             Gross Domestic Product
#Seasonal Adjustment: Seasonally Adjusted Annual Rate
#Frequency:           Quarterly
#Units:               Billions of Chained 2009 Dollars
#Date Range:          1947-01-01 to 2014-07-01
#Last Updated:        2014-11-25 7:56 AM CST
#Notes:               A Guide to the National Income and 
#                     Product Accounts of the United States 
#                     (NIPA) - 

  Quandl.api_key("Provide your quandl API-key")

  mydata<-Quandl(c("FRED/GDPC96"),start_date=start_date,end_date=end_date,type="xts")
} else
{
# Load mydata (GDP-series)
  load("mydata")
}

# Prepare GDP-data
start_year<-1960
end_date<-format(Sys.time(), "%Y-%m-%d")
end_year<-as.double(substr(end_date,1,4))
start_date=paste(start_year,"-01-01",sep="")
# Select data between start_year and end_year
data_sample<-mydata[paste("/",end_date,sep="")]
data_sample<-data_sample[paste(start_date,"/",sep="")]
lgdp <- ts(100*log(data_sample),start=start_year,frequency=4)
nobs <- length(lgdp)


# The function hpFilt derives the MA-coefficients of the implicit ARIMA(0,2,2)-model underlying the HP-filter
#   -The HP-filter assumes a particular ARIMA(0,2,2)-model for the data-gerenating process
#   -The MA-coeffcients of the (twice differenced) stationary data are determined by lambda
#   -See 
head(hpFilt)

x<-lgdp
# Series length
len<-L_hp<-length(x)
# Select lambda
lambda_hp<-1600
q<-1/lambda_hp

hp_filt_obj<-hpFilt(q,L_hp)

tail(hpFilt,2)
hp_filt_obj$ma_model
ma_coeff<-hp_filt_obj$ma_model[2:3]



K<-2*len
K<-1200
# proceed to filtering
x_hp <- hpfilter(x,type="lambda", freq=lambda_hp)
# Extract the coefficients of the symmetric trend:
#   hpfilter generates coefficients of the HP-gap (see below):
#   we here transform back to trend filter
parm<-diag(rep(1,len))-x_hp$fmatrix

# Plots: filter coefficients and series
par(mfrow=c(2,1))
title_more<-NA
mplot<-cbind(parm[,len/2],parm[,1])
plot_title<-"HP lambda=1600: symmetric (red) and real-time (blue)"
axis_d<-1:len-1
insamp<-1.e+99
colo<-c("red","blue")
mplot_func(mplot,axis_d,plot_title,title_more,insamp,colo)
mplot<-cbind(rep(NA,len),rep(NA,len),rep(NA,len),x,x_hp$trend)
plot_title<-"Log US-GDP (blue) vs HP-Trend (red)"
plot(mplot[,4],col="blue",xlab="",ylab="",main=plot_title)
nberShade()
lines(mplot[,5],col="red")



# Compute pseudo-spectral density underlying Wiener-Kolmogorov derivation of HP 
#   (see McElroy (2008) or Maravall-Kaiser p.179)
# For lambda=1600 the MA coefficients are -1.77709 and 0.79944
# Note that MDFA is fed with the square-root of the spectrum
#   (this would correspond to the absolute value of the DFT)
weight_func_h<-abs((1+ma_coeff[1]*exp(-1.i*(0:(K))*pi/(K))+
                      ma_coeff[2]*exp(-1.i*2*(0:(K))*pi/(K)))/
                     (1-exp(-1.i*(0:(K))*pi/(K)))^2)
# Specify (square-root) spectra of target (first column) 
#   and of explanatory variable (second column): target and 
#   explanatory are the same here (univariate filter)
weight_func<-cbind(weight_func_h,weight_func_h)
# Compute target Gamma: HP-trend symmetric filter, see McElroy (2008)
Gamma<-0:(K)
for (k in 0:(K))
{
  omegak<-k*pi/(K)
  Gamma[k+1]<-(1/lambda_hp)/(1/lambda_hp+abs(1-exp(1.i*omegak))^4)
}


par(mfrow=c(2,1))
colo<-c("blue","red")
insamp<-1.e+99
mplot<-as.matrix(Gamma)
plot_title<-"Target Gamma, lambda=1600"
freq_axe<-rep(NA,K+1)
freq_axe[1]<-0
freq_axe[1+(1:6)*K/6]<-c(paste(c("",2:5),"pi/6",sep=""),"pi")
mplot_func(mplot,freq_axe,plot_title,title_more,insamp,colo)
# Plot log spectrum: weight_func must be squared
mplot<-as.matrix(c(NA,log(weight_func[2:(K+1),1]^2)))
plot_title<-"Log pseudo-spectrum, lambda=1600"
mplot_func(mplot,freq_axe,plot_title,title_more,insamp,colo)





weight_func_hp<-weight_func
K<-nrow(weight_func_hp)-1
# Frequency zero is infinity (unit root)
#   The singularity is removed by imposing first and second order 
#     restrictions
#   For numerical computations we set the spectrum arbitrarily 
#     to zero in freq. zero
weight_func_hp[1,]<-0
# Filter length is identified with sample length
L<-len
# Set default settings for MDFA (MSE, no regularization)
# First and second order constraints are imposed
#   (the level constraint weight_constraint=1 is set in the default settings)
i1<-T
i2<-T
# Cutoff: the frequency at which the target drops below 0.5
cutoff<-pi*which(Gamma<0.5)[1]/length(Gamma)
# Real-time (nowcast)
Lag<-0
# Alternative (identical) context-specific estimation: 
weight_constraint<-1
shift_constraint<-0

imdfa_hp<-MDFA_mse_constraint(L,weight_func_hp,Lag,Gamma,i1,i2,weight_constraint,shift_constraint)$mdfa_obj

par(mfrow=c(2,1))
colo<-c("blue","red")
insamp<-1.e+99
mplot<-cbind(imdfa_hp$b,parm[1:L,max(0,Lag)+1])
rownames(mplot)<-paste("Lag ",0:(nrow(mplot)-1))
colnames(mplot)<-c("Replication by DFA","HP-real-time")
plot_title<-"Replication HP-real-time by DFA: Lags 0-240"
freq_axe<-rownames(mplot)
title_more<-c("DFA","HP")
mplot_func(mplot,freq_axe,plot_title,title_more,insamp,colo)
mplot<-mplot[1:21,]
rownames(mplot)<-paste("Lag ",0:(nrow(mplot)-1))
colnames(mplot)<-c("Replication by DFA","HP-real-time")
plot_title<-"Replication HP-real-time by DFA: Lags 0-20"
freq_axe<-rownames(mplot)
title_more<-c("DFA","HP")
mplot_func(mplot,freq_axe,plot_title,title_more,insamp,colo)


