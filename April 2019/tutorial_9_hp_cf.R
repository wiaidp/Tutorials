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

#-------------
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

#---------------------------------------------------------------------------------
# Replication Hodrick-Prescott (HP-) filter

# The function hpFilt derives the MA-coefficients of the implicit ARIMA(0,2,2)-model underlying the HP-filter
#   -The HP-filter assumes a particular ARIMA(0,2,2)-model for the data-gerenating process
#   -The MA-coeffcients of the (twice differenced) stationary data are determined by lambda
#     -See paper by McElroy (in literature folder on github) for details
#   -The ARIMA(0,2,2) is necessary for deriving weight_func, the pseudo-spectrum, of MDFA
head(hpFilt)

x<-lgdp
# Series length
len<-L_hp<-length(x)
# Select lambda: 1600 is typically used for quarterly data 
lambda_hp<-1600
q<-1/lambda_hp

# This function derives the MA-coefficients of the ARIMA(0,2,2)-model: it is not related to MDFA
hp_filt_obj<-hpFilt(q,L_hp)

tail(hpFilt,2)
hp_filt_obj$ma_model
ma_coeff<-hp_filt_obj$ma_model[2:3]
# The MA-coefficients depend on lambda
ma_coeff


#----------------------
# Apply the HP-filter to GDP
#   We here rely on the R-package mFilter 
#   Below, we will replicate the HP-filter by MDFA

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

#-------------------------
# Replicate HP-filter by MDFA
#   -For that purpose we need weight_func (the pseudo-spectrum of the ARIMA(0,2,2)) as well as Gamma (the target symmetric filter)

# Select resolution of frequency-grid
K<-1200

# 1.Compute pseudo-spectral density underlying Wiener-Kolmogorov derivation of HP 
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

# 2.Compute target Gamma: HP-trend symmetric filter, see McElroy (2008)
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


#-------------------
# Set-up MDFA
#   -We here rely on the MSE-wrapper MDFA_mse_constraint 
#     -The constraints assume that
#       1. Shift of one-sided filter vanishes in frequency-zero: i2<-T, shift_constraint<-0
#       2. Amplitude of one-sided filter is one in frequency zero: i1<-T, weight_constraint<-1
#       It is necessary to impose these constraints because of the implicit ARIMA(0,2,2)-model assumed by HP-filter
#   -Note that the pseudo-spectrum of the ARIMA(0,2,2) is infinite in frequency zero which would shut-down the numerical optimization
#     -We skip this issue by assigning an arbitrary (zero-) value in frequency-zero
#     -This arbitrary value could be ignored because we impose constraints i1==T and i2==T in frequency zero
#     -Any other arbitrary number could be used without affecting the resulting estimate
#   -We want a nowcast: Lag=0
#   -cutoff is not relevant in a MSE-setting (it's relevant for customization): we could assign any value to cutoff

weight_func_hp<-weight_func
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
# Constraints: amplitude is one, shift is zero
weight_constraint<-1
shift_constraint<-0

imdfa_hp<-MDFA_mse_constraint(L,weight_func_hp,Lag,Gamma,i1,i2,weight_constraint,shift_constraint)$mdfa_obj

# Compare coefficients of one-sided HP-filter by R-package mFilter and by MDFA 
#   -Both series are virtually indistinguishable
#   -The minuscule differences could be reduced further (arbitrarily small) by selecting a higher resolution of the frequency-grid above (at cost of numerical speed)
# Thus the HP-filter could be replicated by MDFA
#   -In principle we could now apply customization (see MDFA-book for illustration)
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


# Check constraints
# first-order: should give 1
print(paste("Transfer function in frequency zero: ",
            round(sum(imdfa_hp$b),3),sep=""))
# second-order: time-shift should vanish
print(paste("Time-shift in frequency zero: ",
            round((1:(L-1))%*%imdfa_hp$b[2:L],10),sep=""))



#--------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Replication Christiano-Fitzgerald (CF-) filter
# CF-filter differs from HP
#   It is an ideal bandpass: this will be the target for MDFA
#   The implicit model assumption is a random-walk: this will be the spectrum for MDFA


x<-lgdp
# Series length
len<-length(x)
# Resolution of frequency-grid
#   Selecting a higher resolution would tighten the approximation of CF by DFA
K<-1200
# Upper and lower cutoffs of bandpass: lengths in quarters
len1<-8
len2<-40
cutoff1<-2*pi/len1
cutoff2<-2*pi/len2
# Specify target bandpass
Gamma_cf<-((0:K)>K*cutoff2/pi)&((0:K)<K*cutoff1/pi)
# Specify (square-root of) implicit pseudo-spectral density (of random-walk)
weight_func_cf<-matrix(rep(1/abs(1-exp(1.i*(0:K)*pi/K)),2),ncol=2)
K<-nrow(weight_func_cf)-1
# Remove singularity in frequency zero (one can assign an arbitrary value)
weight_func_cf[1,]<-0
# Filter length
L<-len
# Filter constraints: the CF-filter assumes an I(1)-model: i1<-T, i2<-F
i1<-T
i2<-F
# Constraints: amplitude is zero (CF-filter is a abndpass); 
# shift is irrelevant because i2<-F: we could assign any arbitrary value to shift_constraint
weight_constraint<-0
shift_constraint<-0

# Proceed to estimation
imdfa_cf<-MDFA_mse_constraint(L,weight_func_cf,Lag,Gamma_cf,i1,i2,weight_constraint,shift_constraint)$mdfa_obj

omega_k<-pi*(0:K)/K
amp_mse<-abs(imdfa_cf$trffkt)
mplot<-as.matrix(amp_mse)
mplot[1,]<-NA
colnames(mplot)<-NA
ax<-rep(NA,nrow(mplot))
ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6",
                                   "4pi/6","5pi/6","pi")
plot_title<-paste("Amplitude of one-sided Christiano Fitzgerald with 
                  cutoffs pi/",len2/2,", pi/",len1/2,sep="")
insamp<-1.e+90
title_more<-dimnames(mplot)[[2]]
colo<-c("blue","red")
mplot_func(mplot, ax, plot_title, title_more, insamp, colo)


print(paste("Transfer function in frequency zero: ",
            round(sum(imdfa_cf$b),3),sep=""))
# Check second-order: time-shift is unconstrained
print(paste("Time-shift in frequency zero: ",
            round((1:(L-1))%*%imdfa_cf$b[2:L],3),sep=""))

#------------------------------------
# We now verify that the above one-sided filter, obtained by MDFA, replicates the one-sided CF-filter
#   -For that purpose we would have relied on the mFilter-package: 
#     -Unfortunately mFilter is wrong, as can be seen below 
#   -Therefore we rely on the classic model-based solution for deriving the one-sided filter and we check that this corresponds to MDFA


# First attempt: based on mFilter
# We here rely on R-package mFilter for the CF-filter: comparisons with the MDFA-package are provided below
x_cf<-cffilter(x,pu=len2,pl=len1,root=F,drift=F, nfix=NULL,theta=1)
parm_cf<-x_cf$fmatrix

# Check first-order constraint: should give 0 i.e. amplitude in frequency zero should vanish
print(paste("Transfer function in frequency zero: ",
            round(sum(parm_cf[,1]),3),sep=""))
  x_cf_T<-cffilter(x,pu=len2,pl=len1,root=T,drift=F, nfix=NULL,theta=1)
parm_cf_T<-x_cf_T$fmatrix

# Check first-order: should give 0
print(paste("Transfer function in frequency zero: ",
            round(sum(parm_cf_T[,1]),3),sep=""))

#  Obviously, the code does not seem to work properly. 
x_mF_T<-mFilter(x,filter="CF",pu=len2,pl=len1,root=T,drift=F, 
                  nfix=NULL,theta=1)
parm_mF_T<-x_mF_T$fmatrix

# Check first-order: should give 0
print(paste("Transfer function in frequency zero: ",
            round(sum(parm_mF_T[,1]),3),sep=""))

#  Neither $mFilter$ nor $cffilter$  comply with the first-order constraint in frequency zero. Selecting $root=T$ results in a severly misspecified filter.
colo<-c("blue","red","green")
mplot<-cbind(imdfa_cf$b,parm_cf[,1],parm_mF_T[,1])
colnames(mplot)<-c("DFA","mFilter: plot=F","mFilter: plot=T")
plot_title<-"Real-time CF-filters: DFA (blue) vs. mFilter (red and green)"
freq_axe<-paste("Lag ",0:(len-1),sep="")
title_more<-colnames(mplot)
mplot_func(mplot,freq_axe,plot_title,title_more,insamp,colo)

# As can be seen, DFA (blue) and $mFilter$ or $cffilter$ based on $root=F$ (red) match closely up to the boundaries, lags 0 and $\Sexpr{len}$, where non-negligible discrepancies can be observed. 
# In contrast, the coefficients of $mFilter$ for $root=T$ (green line) are `off the mark'.

#----------------------
# Cross-check real-time DFA coefficients with the (time-domain) model-based solution proposed in section \ref{time_domain}, see fig.\ref{z_CF_us_real_log_gdp_fb}. Hint: our R-code relies on \ref{mba_coef_td}, whereby $a_1=1$ (random-walk).
#   Finally we can verify replication of CF-filter by MDFA
ord<-100000
b<-0:ord
b[1+1:ord]<-(sin((1:ord)*2*pi/len1)-sin((1:ord)*2*pi/len2))/(pi*(1:ord))
b[1]<-2/len1-2/len2
# Real-time filter based on for- and backcasts
b_finite<-b[1:len]
# The lag-0 coefficient is augmented by forecasts
b_finite[1]<-b_finite[1]+sum(b[2:ord])
# The lag-len coefficient is augmenetd by backcasts
b_finite[len]<-b_finite[len]+sum(b[(len+1):ord])
# Compare DFA and model-based coefficients
#   Both coefficients overlap almost perfectly: tighter approximations could be obtained 
#     by increasing K (resolution of frequency-grid)
colo<-c("red","blue")
mplot<-cbind(b_finite,imdfa_cf$b)
plot_title<-"Real-time CF-filters: Forecast/backcast (red) vs. DFA (blue)"
freq_axe<-rep(NA,len)
freq_axe[1]<-0
freq_axe[(1:6)*len/6]<-paste("Lag ",as.integer(1+(1:6)*len/6),sep="")
mplot_func(mplot,freq_axe,plot_title,title_more,insamp,colo)
