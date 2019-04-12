# Todos
#   -Example with leading indicator
#   -Add/discuss constraints




# Purpose of tutorial: 

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
# Example 1 bivariate leading indicator design
# -multivariate designs are particularly useful when incorporating a leading indicator into the set of explanatory variables 
#   (assuming that such a series exists...)
# -we here propose a simulation experiment based on a bivariate design where the leading indicator is a noisy version
#   of the target data, shifted one time point ahead
# -we then analyze the outcome of a bivariate MDFA vs. a univariate DFA (filter coefficients)
# -we compute out-of-sample (time-domain) MSEs as well as in-sample time-domain as well as frequency-domain MSEs 
#   (the latter frequency-domain MSEs are the (M)DFA criterion values i.e. optimal coefficients minimize these expressions)
# -we analyze various settings: in-sample length, L (2*L for bivariate design), data-generating process, noisiness (of leading indicator) 
# -The analysis is splitted into three 'toggle breakpoints' where specific outcomes are analyzed

# General design
#  -MSE
#  -Bivariate vs. univariate


# Total sample size (long in order to obtain reliable time-domain MSE-performances)
#   Also: the non-causal ideal lowpass needs more data (assumes knowledge of future observations)
lenh<-10000
# In-sample span: for increasing ratios L/len overfitting will be magnified: 
#   note that bivariate filter requires 2*L coefficients: more prone to overfitting than univariate (for identical L)
len<-120
# Specify AR-process: positive, none, negative autocorrelation (though the latter is rather uncommon in typical economic data: negative autocorrelation may arise from overdifferencing the data)
a1<-0.9
a1<-0.1
# Generate series for each AR(1)-process
set.seed(1)
# We generate lenh+1 data points because construction of the leading indicator will consume one observation
x_long<-arima.sim(list(ar=a1),n=lenh+1)
set.seed(2)
# Generate a leading indicator as second explanatory variable: leading indicator is x_long+scale_idiosyncratic*noise shifted
#   one time point in the future
# Scaling of the idiosyncratic noise
scale_idiosyncratic<-10.4
if (abs(scale_idiosyncratic)==0)
  print("Design is not uniquely specified since both explanatory series are identical (up to a shift)")
eps<-rnorm(lenh+1)
indicator<-x_long+sqrt(var(x_long))*scale_idiosyncratic*eps
# Data: first column=target, second column=x, third column=shifted (leading) indicator
data_matrix<-na.exclude(cbind(x_long[1:(lenh)],x_long[1:lenh],c(indicator[2:(lenh+1)])))
dimnames(data_matrix)[[2]]<-c("target","x","leading indicator")
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
# Extract in-sample span from the long sample: the first M observations are lost for initialization of ideal filter
#   If we want to compare in-sample time-domain and in-sample frequency domain MSEs then we need to skip the first M observations
data_matrix_in_sample<-data_matrix[M:(M+len),]
# One can see that the second explanatory (third column) leads the target data (but it is noisy)
tail(round(data_matrix_in_sample,4))

# Specify target
periodicity<-6
cutoff<-pi/periodicity
# Ideal filter: used for evaluating time-domain MSE
# Length of ideal lowpass (M is the half-length: effective length is 2*M-1 since filter is symmetric)
M<-100
y<-ideal_filter_func(periodicity,M,x_long)$y[1:lenh]

# Spectrum
# One can use either the function per
weight_func_bivariate<-cbind(per(data_matrix_in_sample[,1],T)$DFT,per(data_matrix_in_sample[,2],T)$DFT,per(data_matrix_in_sample[,3],T)$DFT)
# Or one can use spec_comp (the latter is slightly more general than per)
#   Here we can supply the full data matrix data_matrix_in_sample at once: spec_comp then returns the multivariate dft 
weight_func_bivariate<-spec_comp(nrow(data_matrix_in_sample), data_matrix_in_sample, 0)$weight_func
# Resolution of frequency-grid
K<-nrow(weight_func_bivariate)-1
# Target (in frequency domain)
Gamma<-(0:(K))<=K*cutoff/pi
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length: L/K should be 'small'
L<-2*periodicity
# Estimate filter coefficients
mdfa_obj_bivariate<-MDFA_mse(L,weight_func_bivariate,Lag,Gamma)$mdfa_obj 
# Filter coefficients
b__bivariate<-mdfa_obj_bivariate$b
dimnames(b__bivariate)[[2]]<-c("x","leading indicator")
dimnames(b__bivariate)[[1]]<-paste("Lag ",0:(L-1),sep="")#dim(b_mat)


# First 'toggle breakpoint': interpreting the filter coefficients
head(b__bivariate)


# Filtering bivariate: used for computing time domain MSEs
yhat_bivariate_leading_indicator<-filt_func(data_matrix[,2:ncol(data_matrix)],b__bivariate)$yhat

# Estimation univariate DFA: benchmark (we know from previous tutorial that this is a tough competitor)
weight_func_univariate<-weight_func_bivariate[,1:2]
mdfa_obj_univariate<-MDFA_mse(L,weight_func_univariate,Lag,Gamma)$mdfa_obj 
# Filter coefficients
b_univariate<-mdfa_obj_univariate$b
# Filtering univariate: this is used for computing time-domain MSEs
yhat_univariate<-filt_func(as.matrix(data_matrix[,2]),b_univariate)$yhat

# Compute MSE-performances
y_target_leading_indicator<-y
anf<-M
enf<-M+len
mse_in<-round(c(mean(na.exclude((yhat_bivariate_leading_indicator-y_target_leading_indicator)[anf:enf])^2),
           mean(na.exclude((yhat_univariate-y_target_leading_indicator)[anf:enf])^2)),3)
anf<-M+1+len
enf<-lenh
mse_out<-round(c(mean(na.exclude((yhat_bivariate_leading_indicator-y_target_leading_indicator)[anf:enf])^2),
                mean(na.exclude((yhat_univariate-y_target_leading_indicator)[anf:enf])^2)),3)
perf_mse<-rbind(mse_out,mse_in,c(round(mdfa_obj_bivariate$MS_error,3),round(mdfa_obj_univariate$MS_error,3)))
colnames(perf_mse)<-c("bivariate MDFA","DFA")
rownames(perf_mse)<-c("Sample MSE out-of-sample (time domain)","Sample MSE in-sample (time domain)","Criterion values in-sample (frequency-domain)")

# Second 'toggle breakpoint': analyze MSEs, in-sample and out-of-sample, time domain and frequency domain
round(perf_mse,3)


# Plot/compare filtered series: ideal lowpass vs. DFA (univariate) and MDFA (bivariate)
par(mfrow=c(1,1))
mplot<-mplot_all<-cbind(yhat_univariate,y,yhat_bivariate_leading_indicator)
ymin<-min(mplot,na.rm=T)
ymax<-max(mplot,na.rm=T)
ts.plot(mplot[,1],main=paste("Sample MSE MDFA: ",ylab="",
                            round(perf_mse[1],3),", DFA: ",round(perf_mse[2],3),sep=""),col="blue",
        ylim=c(ymin,ymax))
lines(mplot[,2],col="red")
lines(mplot[,3],col="green")
mtext("DFA", side = 3, line = -2,at=len/2,col="blue")
mtext("target", side = 3, line = -1,at=len/2,col="red")
mtext("MDFA", side = 3, line = -3,at=len/2,col="green")

# Third 'toggle breakpoint': analyze filtered series (smoothness, lead)
#   Zoom into above plot
anf<-300
enf<-400
mplot<-mplot_all[anf:enf,]
ymin<-min(mplot,na.rm=T)
ymax<-max(mplot,na.rm=T)
ts.plot(mplot[,1],main=paste("Sample MSE MDFA: ",ylab="",
                             round(perf_mse[1],3),", DFA: ",round(perf_mse[2],3),sep=""),col="blue",
        ylim=c(ymin,ymax))
lines(mplot[,2],col="red")
lines(mplot[,3],col="green")
mtext("DFA", side = 3, line = -2,at=len/2,col="blue")
mtext("target", side = 3, line = -1,at=len/2,col="red")
mtext("MDFA", side = 3, line = -3,at=len/2,col="green")

# Comments
#   -Bivariate leading indicator design (MDFA, green line) tends to anticipate turning-points (lies to the left of univariate DFA, blue line)
#   -Both outputs of causal filters (blue/green) are noisier than target (red) and are delayed
# Customization
#   -Reduce delay
#   -Reduce noise leakage
#   -Both?



scale_idiosyncratic<-10.4# weights of second series are zero, in-sample MDFA better, out-of-sample worse
scale_idiosyncratic<-0.4


len<-120
len<-300
len<-1000

L<-2*periodicity
L<-4*periodicity
L<-8*periodicity








