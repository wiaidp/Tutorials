# Todos
#   -Add a model-based example
#   -Show that reg replicates MSE when all lambdas are zero
#   -Amplitude, shift, coeff, spectrum, target


# Here we use MDFA_reg-function


rm(list=ls())

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)


#-----------------------------------------------------------------------------------------------
# Source common functions

source("Common functions/plot_func.r")
source("Common functions/mdfa_trade_func.r")
source("Common functions/data_load_functions.R")


#-----------------------------------------------------------------------------------------------
# Data: FX 

data_from_IB<-T
hour_of_day<-"16:00"
i_series_vec<-c(1,2,3,6,7,8)
reload_sp500<-F
path.dat<-paste(getwd(),"/",sep="")

data_load_obj<-data_load_gzd_trading_func(data_from_IB,hour_of_day,reload_sp500,path.dat)

mydata<-data_load_obj$mydata

FX_mat<-data_load_obj$xts_data_mat[,i_series_vec]
log_FX_mat<-log(FX_mat)

colnames(log_FX_mat)

plot_T<-T
anf_plot<-"2000-10-01/"



#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Here we apply to FX what we learned in regularization and customization tutorials
# Multivariate design: target EURUSD, all series are used as explanatory data 
# Fundamental idea
#   Unconstrained designs cannout account 'directly' for a priori knowledge 
#   Imposing constraints (regularization) might go hand-in-hand with applying a priori knowledge
#   For that purpose it is wise to format the data in view of complying with constraints 
# Example: use a priori knowledge for defining the data-matrix
#   -Sine we target EURUSD we here use only pairs with either EUR or USD (i.e. we skip GBPJPY)
#   -We change signs for 'inverted' explanatory (i.e. we change sign of USDJPY)
# In this context regularization (in particular cross-sectional regularization) could 
#   be meaningful beyond the straightforward 'shrinkage' of the parameter space 
in_sample_span<-"2018-01-01"
datah<-na.exclude(diff(log_FX_mat[,c("EURUSD","EURUSD","EURJPY","GBPUSD","USDJPY","EURGBP")]))
head(datah)
data<-datah
data[,"USDJPY"]<--datah[,"USDJPY"]
colnames(data)[5]<-"JPYUSD"
head(data)
# Compute multivariate spectrum
for (i_series in 1:ncol(data))
{
  x<-data[,i_series]
  # Spectrum: white noise assumption
  
  if (i_series==1)
  {
    weight_func_mat<-per(x[paste("/",in_sample_span,sep="")],T)$DFT  
  } else
  {
    weight_func_mat<-cbind(weight_func_mat,per(x[paste("/",in_sample_span,sep="")],T)$DFT)
  }
}

# First column is target; columns 2-6 are explanatory; columns 1 and 2 are identical because EURUSD is also explanatory;
#   column 5 corresponds to the DFT of -USDJPY (negative sign)
head(weight_func_mat)
weight_func<-weight_func_mat

# Set all remaining parameters as in previous tutorials above
source("Common functions/parameter_set.r")

# Imposing regularity: 
#-------------------------
# Example 1
#   We start with the most important one: decay
#   We can impose a strong regularization with a reasonably fast decay 
#     This will freeze a lots of unnecessary degrees of freedom (alternatively we may select a smaller filter length L)
#   Fundamental idea: if target is mid/short-term, then remote past of data is not useful
lambda_decay<-c(0.3,0.99)
# Cross and smooth are set to zero
lambda_cross<-0.00
lambda_smooth<-0.
# MDFA_reg: wrapper for working with regularization
mdfa_reg_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj

par(mfrow=c(1,1))
# Filter coefficients: largest weight is attributed to EURUSD, coefficients are not excessively noisy and the decay rapidly
ts.plot(mdfa_reg_obj$b,col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])
# Amplitude: amplitudes are intuitively appealing
ts.plot(abs(mdfa_reg_obj$trffkt),col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])
# Shift
ts.plot(Arg(mdfa_reg_obj$trffkt)/((0:K)*pi/K),col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])

# Comments
#   1. The coefficients look similar already (without imposing cross-sectional similarity)
#   2. The amplitude functions look fine (shape corresponds to targeted cutoff)
#   3. The shifts look OK, too

#----------
# Trading with the above designs
data_filter<-data[,2:ncol(data)]
# Compute performance on next day: 
#   -lag_fx must be >= 1
#   -if lag_fx>1 then execution will be delayed by lag_fx days
lag_fx<-1

# MSE: 
lambda<-0
eta<-0

mdfa_trade_obj<-mdfa_reg_trade_func(K,periodicity,L,Lag,lag_fx,data_filter,plot_T,weight_func,lambda_cross,lambda_decay,lambda_smooth,lambda,eta)


# Customization
lambda<-5
eta<-0.3

mdfa_trade_obj<-mdfa_reg_trade_func(K,periodicity,L,Lag,lag_fx,data_filter,plot_T,weight_func,lambda_cross,lambda_decay,lambda_smooth,lambda,eta)


#--------------------------------------------------------------------------------
# Example 2: as above with additional cross-sectional tightness

lambda_decay<-c(0.3,0.99)
# No smoothness
lambda_smooth<-0.
# Strong similarity
lambda_cross<-0.9
# MDFA_reg: wrapper for working with regularization
mdfa_reg_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj

par(mfrow=c(1,1))
# Filter coefficients: very similar across series, decay papidly, fairly smooth
ts.plot(mdfa_reg_obj$b,col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])
# Amplitude
ts.plot(abs(mdfa_reg_obj$trffkt),col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])
# Shift
ts.plot(Arg(mdfa_reg_obj$trffkt)/((0:K)*pi/K),col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])

# Comment: 
#   lambda_smooth is not critical: any value in [0.5,0.99] just seems fine
#   In case of doubt: larger might be marginally better (freezing unnecessary degrees of freedom)
#   We might select lambda_cross>0, too, but the coefficients look just fine 'as is'

# Trading with the above designs
data_filter<-data[,2:ncol(data)]
# Compute performance on next day: 
#   -lag_fx must be >= 1
#   -if lag_fx>1 then execution will be delayed by lag_fx days
lag_fx<-1

# MSE: 
lambda<-0
eta<-0

mdfa_trade_obj<-mdfa_reg_trade_func(K,periodicity,L,Lag,lag_fx,data_filter,plot_T,weight_func,lambda_cross,lambda_decay,lambda_smooth,lambda,eta)


# Customization
lambda<-5
eta<-0.3

mdfa_trade_obj<-mdfa_reg_trade_func(K,periodicity,L,Lag,lag_fx,data_filter,plot_T,weight_func,lambda_cross,lambda_decay,lambda_smooth,lambda,eta)


#------------
# Example 3: as above with additional strong smoothness

lambda_decay<-c(0.3,0.99)
# Cross and smooth are set to large values
#   Strong cross means that filter coefficients are similar: common multivariate filter
#   Extract information from all series to infer the common filter coefficients 
lambda_cross<-0.9
lambda_smooth<-0.99

# MDFA_reg: wrapper for working with regularization
mdfa_reg_obj<-MDFA_reg(L,weight_func_mat,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj

par(mfrow=c(1,1))
# Filter coefficients: decay rapidly, are smooth and similar
ts.plot(mdfa_reg_obj$b,col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])
# Amplitude
ts.plot(abs(mdfa_reg_obj$trffkt),col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])
# Shift
ts.plot(Arg(mdfa_reg_obj$trffkt)/((0:K)*pi/K),col=rainbow(ncol(data)-1))
for (i in 1:(ncol(data)-1))
  mtext(colnames(data)[i+1],side=3,line=-i,col=rainbow(ncol(data)-1)[i])

# Comment: 
#   lambda_smooth is not critical: any value in [0.5,0.99] just seems fine
#   In case of doubt: larger might be marginally better (freezing unnecessary degrees of freedom)
#   We might select lambda_cross>0, too, but the coefficients look just fine 'as is'


# Trading with the above designs
data_filter<-data[,2:ncol(data)]
# Compute performance on next day: 
#   -lag_fx must be >= 1
#   -if lag_fx>1 then execution will be delayed by lag_fx days
lag_fx<-1

# MSE: 
lambda<-0
eta<-0

mdfa_trade_obj<-mdfa_reg_trade_func(K,periodicity,L,Lag,lag_fx,data_filter,plot_T,weight_func,lambda_cross,lambda_decay,lambda_smooth,lambda,eta)


# Customization
lambda<-5
eta<-0.3

mdfa_trade_obj<-mdfa_reg_trade_func(K,periodicity,L,Lag,lag_fx,data_filter,plot_T,weight_func,lambda_cross,lambda_decay,lambda_smooth,lambda,eta)

