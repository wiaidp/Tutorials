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
source("Common functions/play_with_bivariate.r")

#---------------------------------------------
# Explaining the experiment: bivariate leading indicator design
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
len<-300
# Specify AR-process: positive, none, negative autocorrelation (though the latter is rather uncommon in typical economic data: negative autocorrelation may arise from overdifferencing the data)
a1<-0.9
# Generate series for each AR(1)-process
set.seed(1)
# We generate lenh+1 data points because construction of the leading indicator will consume one observation
x_long<-arima.sim(list(ar=a1),n=lenh+1)
set.seed(2)
# Generate a leading indicator as second explanatory variable: leading indicator is x_long+scale_idiosyncratic*noise shifted
#   one time point in the future
# Scaling of the idiosyncratic noise
scale_idiosyncratic<-.4
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
ts.plot(mplot[,1],main=paste("Out-of-sample MSE MDFA: ",ylab="",
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
ts.plot(mplot[,1],main=paste("Out-of-sample MSE MDFA: ",ylab="",
                             round(perf_mse[1],3),", DFA: ",round(perf_mse[2],3),sep=""),col="blue",
        ylim=c(ymin,ymax))
lines(mplot[,2],col="red")
lines(mplot[,3],col="green")
mtext("DFA", side = 3, line = -2,at=len/2,col="blue")
mtext("target", side = 3, line = -1,at=len/2,col="red")
mtext("MDFA", side = 3, line = -3,at=len/2,col="green")


#--------------------------------------------------------------------------------
# Example 1
#   This example replicates the above lengthy code
#   For that purpose we stiffed the above lengthy code into a function play_bivariate_func and we play with
#     -a1, scale_idiosyncratic (noisiness of leading indicator), len (in-sample span) and L (filter length: bivariate filter requires 2*L degrees of freedom) 

# Strong autocorrelation
a1<-0.9
# Noisy but not too much so
scale_idiosyncratic<-0.4
# Fairly long in-sample span
len<-300
# Modest filter length (periodicity is 6 so L=12 is needed to damp components in the stopband)
L<-12

play_obj<-play_bivariate_func(a1,scale_idiosyncratic,len,L)
  
play_obj$b__bivariate
play_obj$perf_mse

# Comments:
#  1. filter coefficients (toggle point 1 ):
#   The coefficients applied to the first series look like an 'ordinary' (one-sided) lowpass
#   The coefficients applied to the second series (leading indicator) assign most weight to last data point
#     This is meaningful because last data-point looks into the future (leading)
#     but since the series is noisy, the weight (coefficient) is not dominating
#  2. MSE performances (toggle point 2):
#   -Both in-sample measures (row 2: time-domain; row 3: frequency-domain criterion value) are close as desired
#     -A heavy mismatch is indicative of overfitting: the frequency-grid would not be dense enough (K/L too small, ill-conditioned)
#   -Out-of-sample performances of both filters are worse (overfitting)
#   -Bivariate filter (first column) outperforms univariate (sedond column) in-sample as well as out-of-sample
#     -This last result is expected at least in-sample; however, out-of-sample outperformance is less trivial because
#      of increased overfitting by the bivariate design (requiring 2*L degrees of freedom)
#  3. Plot of filter outputs (toggle point 3):
#   -The output of the bivariate filter is generally anticipating turning-points, as expected
#   -Both outputs (of one-sided filters) are slightly delayed with respect to target (though sometimes the one-sided filters seem to anticipate target)
#     -This is because the data is smooth already (strong autocorrelation) so that the filtering-task is not too difficult
#--------------------------------------------------------------------------------
# Example 2: same as above but we make the filtering task more difficult by specifying a 'close to white' process
# Almost white noise (close to fitting an ARMA-model to log-returns of EURUSD)
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

# Comments:
#  1. filter coefficients (toggle point 1 ):
#   The coefficients applied to the first series are decaying more slowly than in example 1: heavier smoothing of the noise
#   As previously, the coefficients applied to the second series (leading indicator) assign most weight to last data point
#     This is meaningful because last data-point looks into the future (leading)
#     but since the series is noisy, the weight (coefficient) is not dominating
#  2. MSE performances (toggle point 2):
#   -Both in-sample measures (row 2: time-domain; row 3: frequency-domain criterion value) are tightly matched, as desired
#     -A heavy mismatch is indicative of overfitting: the frequency-grid would not be dense enough (K/L too small, ill-conditioned)
#   -Out-of-sample performances of both filters are worse (overfitting) but less so (proportionally) than in example 1
#     Heavier smoothing mitigates (a bit) overfitting
#   -Efficiency gains by bivariate design are commensurate with example 1
#  3. Plot of filter outputs (toggle point 3):
#   -Filtered series of one-sided designs are less smooth than in example 1 (input data is much noisier) and slightly more delayed (due to heavier smoothing)
#   -bivariate filter generally leads univariate design by one time unit, as desired

#--------------------------------------------------------------------------------
# Example 3: same as example 2 (noisy data) but we make the leading indicator very noisy: thus we expect the univariate DFA to work as well as bivariate MDFA
# Almost white noise (close to fitting an ARMA-model to log-returns of EURUSD)
a1<-0.08
# Very noisy leading indicator
scale_idiosyncratic<-10
# Fairly long in-sample span
len<-300
# Modest filter length (periodicity is 6 so L=12 is needed to damp components in the stopband)
L<-12

play_obj<-play_bivariate_func(a1,scale_idiosyncratic,len,L)

play_obj$b__bivariate
play_obj$perf_mse

# Comments:
#  1. filter coefficients (toggle point 1 ):
#   The coefficients applied to the first series are decaying similarly to example 2
#   In contrast to previous example 2, the coefficients applied to the second series are close to zero
#     This is meaningful because the second explanatory series is independent noise: there is no useful information for improving estimation
#  2. MSE performances (toggle point 2):
#   -Both in-sample measures (row 2: time-domain; row 3: frequency-domain criterion value) are tightly matched, as desired
#     -A heavy mismatch is indicative of overfitting: the frequency-grid would not be dense enough (K/L too small, ill-conditioned)
#   -Out-of-sample performances of both filters are worse (overfitting) as in example 2 
#   -The bivariate filter marginally outperforms the univariate filter in-sample; but it is marginally outperformed out-of-sample
#     -This is 'as expected' due to overfitting
#     -However we note that the loss of out-of-sample performance is small which is fine (the more general bivariate design does not loose much when 'unnecessary' explanatory data is added) 
#  3. Plot of filter outputs (toggle point 3):
#   -Both one-sided filter outputs nearly overlap, as expected ideally

#--------------------------------------------------------------------------------
# Example 4: same as example 2 but we make the in-sample span much shorter in order to magnify overfitting
# Almost white noise (close to fitting an ARMA-model to log-returns of EURUSD)
a1<-0.08
# Noisy but not too much so
scale_idiosyncratic<-0.4
# Short sample: might be useful for non-stationary data where the generating process changes over time
len<-100
# Modest filter length (periodicity is 6 so L=12 is needed to damp components in the stopband)
#   But note that we fit 2*12=24 coefficients to 100 observations only (49 complex-valued frequency-domain equations + 2 real-valued at frequencies 0 and pi which gives 49*2+2=100 frequency-domain observations)
L<-12

play_obj<-play_bivariate_func(a1,scale_idiosyncratic,len,L)

play_obj$b__bivariate
play_obj$perf_mse

# Comments:
#  1. filter coefficients (toggle point 1 ):
#   -Remarkably close to results in example 2
#  2. MSE performances (toggle point 2):
#   -In-sample measures (row 2: time-domain; row 3: frequency-domain criterion value) differ: 
#     -the criteria values (third row) are too small because both filters overfit the coarse freqeuncy-grid (only 50 frequencies available)
#     -the bivariate filter (first column) 'overfits' more hevaily, as expected
#   -Out-of-sample performances of both filters are worse (overfitting) as in example 2 
#     -However, a comparison of out-of-sample MSEs (example 2 and here) reveals pretty similar performances: which is rather unexpected but welcome
#     -Efficiency gains of the bivariate filter are similar to example 2 above which is, once again, rather unexpected but welcome
#  3. Plot of filter outputs (toggle point 3):
#   -The bivariate design is faster out-of-sample, es desired (but not necessarily as expected given that we estimate 24 coefficients on a sample of length 100)

#--------------------------------------------------------------------------------
# Example 5: same as example 2 but 'extreme' overfitting: 2*L/len=0.5 i.e. the bivariate design fits 50 coefficients to 100 data-points 
#   What happens with out-of-sample performances?

# Almost white noise (close to fitting an ARMA-model to log-returns of EURUSD)
a1<-0.08
# Noisy but not too much so
scale_idiosyncratic<-0.4
# Short sample (a bit more than a quarter of daily data): 
#   might be interesting in the context of non-stationary data where the generating process changes rapidly over time
len<-100
# Too much...
L<-25

play_obj<-play_bivariate_func(a1,scale_idiosyncratic,len,L)

play_obj$b__bivariate
play_obj$perf_mse
yhat_mat_overfitting<-play_obj$mplot_all
# Comments:
#  1. filter coefficients (toggle point 1 ):
#   -Not too ugly (still interpretable: the largest weight of the second filter is still assigned to the last data point)
#  2. MSE performances (toggle point 2):
#   -Out-of-sample performances of bivariate still a smidge better than univariate 
#  3. Plot of filter outputs (toggle point 3):
#   -The bivariate design is still faster out-of-sample, es desired (but not necessarily as expected given that we estimate 50 coefficients on a sample of length 100)

#--------------------------------------------------------------------------------
# Example 6: to be contrasted with example 5 above
#   As example 5 (or 2) but we use a very long sample in order to reveal 'asymptotic' behaviour
#   Expectation: for each filter all MSE-numbers should be tightly matched

# Almost white noise (close to fitting an ARMA-model to log-returns of EURUSD)
a1<-0.08
# Noisy but not too much so
scale_idiosyncratic<-0.4
# Very long sample (more than 10 years of daily data)
len<-3000
# moderate filter length
L<-50

play_obj<-play_bivariate_func(a1,scale_idiosyncratic,len,L)

play_obj$b__bivariate
play_obj$perf_mse

# Comments:
#  1. filter coefficients (toggle point 1 ):
#   -as expected 
#  2. MSE performances (toggle point 2):
#   -as expected: all MSE numbers are tightly matched (per filter)
#  3. Plot of filter outputs (toggle point 3):
#   -as expected

#--------------------------------------------------------------------------------
# Example 7: same as example 6 (asymptotics) above but we illustrate that longer filter lengths do not improve performances much asymptotically
#   As a consequence our rule of selecting L<-2*periodicity is OK in particular in cases where data is sparse and/or the generating-process is non-stationary (changing rapidly over time)

# Almost white noise (close to fitting an ARMA-model to log-returns of EURUSD)
a1<-0.08
# Noisy but not too much so
scale_idiosyncratic<-0.4
# Very long sample (more than 10 years of daily data)
len<-3000
# moderate filter length
L<-50

play_obj<-play_bivariate_func(a1,scale_idiosyncratic,len,L)

play_obj$b__bivariate
play_obj$perf_mse
yhat_mat_asymptotic<-play_obj$mplot_all



# To conclude we here compare the filter outputs of the 'massively overfitted' MDFA (example 5 above) and of the 
#   best possible 'asymptotic'  MDFA and DFA (very large in-sample span, large L)
anf<-(nrow(yhat_mat_asymptotic)-100)
enf<-nrow(yhat_mat_asymptotic)
mplot<-cbind(yhat_mat_asymptotic[,3],yhat_mat_overfitting[,3],yhat_mat_asymptotic[,1])[anf:enf,]
ts.plot(mplot[,1],ylim=c(min(mplot),max(mplot)),main="Asymptotic MDFA (blue) and DFA (green) with large L vs. massively overfittef (red) bivariate MDFAs",col="blue")
lines(mplot[,2],col="red")
lines(mplot[,3],col="green")
mtext("Asymptotic MDFA",line=-1,col="blue")
mtext("Asymptotic DFA (univariate)",line=-2,col="green")
mtext("Massively overfitted MDFA",line=-3,col="red")
abline(h=0)

# Comment:
#   The overfitted design performs remarkably well... 





