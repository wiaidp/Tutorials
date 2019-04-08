# Todos
# Specify target in frequency-domain and time-domain
# Use different spectrum: flat, AR (discuss fit), customized/tweaked
# Play with Lag
# Try a smaller L



# In previous tutorial (DFA and forecasting) we emphasized a particular target: Gamma was an allpass filter
# Here we propose generic targets, including lowpass, bandpass, highpass, Hodrick-Prescott, arbitrary,....

# Purpose of tutorial: illustrate (M)DFA user-interface
# Design: univariate, MSE (no customization), unconstrained (no regularization)
#   -Customization and regularization will be tackled in separate tutorials



# Start from scratch
rm(list=ls())



library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
#devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)

# Gamma and weight_func in the head of MDFA_mse allow for powerful interaction of the user with the estimation algorithm
head(MDFA_mse)

#-----------------------------------------------------------------------------------------------
# Source common functions

source("Common functions/plot_func.r")
source("Common functions/arma_spectrum.r")
source("Common functions/ideal_filter.r")
source("Common functions/mdfa_trade_func.r")


#-------------------------------------------------------------------------------------------------
# Play with target of DFA
#   DFA requires the user to supply a 'target' and a 'spectrum'; DFA then returns an optimal estimate 
#     -Optimality: Mean-square error (MSE) or beyond MSE (ATS-trilemma)
# In the following examples: illustration of various targets
#   This is more general than state-of-the-art forecast approaches
# Design:
#   MSE
#   Univariate
#   White noise spectrum

# DFA is specified in the frequency-domain, see McElroy/Wildi
#   -K defines the frequency-grid for estimation: all frequencies omega_j=j*pi/K, j=0,1,2,...,K are considered
#   -Tradeoff: larger K (denser grid) means better MSE (if true process is white noise) but longer computation-time
#   -We here consider a partition of the interval [0,pi] in K=600 equally-distanced frequency supports
K<-600
# Spectrum: white noise assumption
#   -First column is spectrum of target, second column is spectrum of explanatory variable
#   -In a univariate design target and explanatory data are the same
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
colnames(weight_func)<-c("target","explanatory")

# White noise: flat spectrum (all frequencies are loaded equally by the process)
plot(weight_func[,1],type="l",main=paste("White noise spectrum, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext(colnames(weight_func)[1],line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

#-------------------
# Example 1: classic one step ahead forecasting (replicates SOTA ARMA-models)
# Target forecasting: 
#   -We are interested in all frequencies equally
#   -Target Gamma is a forward-looking allpass filter in the frequency-domain
Gamma<-rep(1,K+1)
plot(Gamma,type="l",main=paste("Allpass forecast target, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()
# One step ahead: Lag=-1
Lag<--1
# Filter length: number of weights/coefficients of forecast filter
L<-10

# Estimation based on MDFA-MSE wrapper
mdfa_obj<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 

# The following function plots 
#   -the coefficients
#   -the time-shift and
#   -the amplitude of the filter
plot_estimate_func(mdfa_obj,weight_func,Gamma)

# Comments
#   -The amplitude function of the MSE-filter is close to zero, as desired 
#   -It is not exactly zero because K is finite
#   -For K/L>10 mismatch is negligible by all practical means 

#--------------------------------------------------------------------------------------
# Example 2: lowpass
#   In contrast to forecasting we are not interested in all frequency-components
#   When specifying a lowpass design, the user generally wants to fade-out undesirable (high-frequency) 'noise'
#   Alternatively, the user might want to emphasize interesting low-frequency components: for example trend and/or cycle

# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
plot(Gamma,type="l",main=paste("Ideal lowpass, periodicity=",periodicity,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Comments
#   The target filter Gamma lets pass through (does not alter) components with durations larger than 2*periodicity
#   The target filter Gamma stops (eliminates) components with durations smaller than 2*periodicity
#   Filter is called 'ideal' lowpass
#   This filter is very easy to specify in the frequency-domain and it is intuitively appealing

# Since the transferfunction Gamma is real-valued, the filter is symmetric (imaginary parts of symmetric filters cancel)
# We now apply the ideal filter to white noise
len<-1000
set.seed(1)
x<-rnorm(len)
# Length of filter (true ideal filter is bi-infinite)
M<-100
# This function computes the filter coefficients and applies the filter to x
id_obj<-ideal_filter_func(periodicity,M,x)

# The coefficients are symmetric: the filter is not causal (it needs future data)  
ts.plot(id_obj$gamma)
# The filtered series (red line in plot below) is smooth: high-frequency noise has been damped
ts.plot(x)
lines(id_obj$y,col="red")
# Mean duration between consecutive zero-crossings of the filtered data (mean holding-time of an asset in trading with this filter)
#   Higher periodicities imply stronger smoothing and thus longer holding-times (less frequent trades)
id_obj$mean_holding_time

# Idea: the user can specify a periodicity which matches his particular interests (more flexible than classic time series approaches)
# Problem: the ideal filter cannot be applied in 'real-time' because it is not causal
# Solution: approximate Gamma by a real-time (causal) filter i.e. use DFA (see examples below)

#--------------------------------------------------------------------------------------
# Example 3: bandpass
#   In contrast to forecasting we are not interested in all frequency-components
#   When specifying a bandpass design, the user generally wants to fade-out undesirable (high- and low-frequency) components
#   Alternatively, the user might want to emphasize interesting cycle components
#   This type of filter is typically in business-cycle analysis, for example

# Target: specify cutoff=pi/periodicity of lowpass ideal target
periodicity_low<-5
cutoff_low<-pi/periodicity_low
periodicity_high<-10
cutoff_high<-pi/periodicity_high
Gamma<-((0:(K))<=K*cutoff_low/pi)-((0:(K))<=K*cutoff_high/pi)+1.e-9


plot(Gamma,type="l",main=paste("Ideal lowpass, periodicities=",periodicity_low,",",periodicity_high,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Comments
#   The target filter Gamma lets pass through (does not alter) components with durations between 2*periodicity_low and 2*periodicity_high
#   The target filter Gamma stops (eliminates) all other components
#   Filter is called 'ideal' bandpass
#   This filter is very easy to specify in the frequency-domain and it is intuitively appealing

# The coefficients of the ideal bandpass can be calculated by using the previous function
#   The coefficients of the bandpass are the difference of the two lowpass filters
# We skip this calculation since it is not of immediate interest

# As for the ideal lowpass, the ideal bandpass is not causal and can be approximated by DFA



#--------------------------------------------------------------------------------------
# Example 4: more complex targets


# Target: specify cutoff=pi/periodicity of lowpass ideal target

Gamma<-rep(1,K+1)

Gamma[3:20]<-0.4
Gamma[37:98]<-0
Gamma[167:208]<-pi
Gamma[398:476]<-0.1

plot(Gamma,type="l",main="Arbitrary target",
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Comments
# As for the previous ideal lowpass and bandpass the above generic filter is not causal and can be approximated by DFA
# The user can specify any target that fits his research interest: more general than classic time series approaches


#-------------------------------------------------------------------------------------------
# Example 5: approximate non-causal ideal lowpass by causal filter
#   Interpret amplitude and time-shift functions

K<-600
# Spectrum: white noise assumption
#  First column is target, second column is explanatory variable: in a univariate design target and explanatory are the same
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
colnames(weight_func)<-c("target","explanatory")
weight_func_noise<-weight_func
periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
plot(Gamma,type="l",main=paste("Ideal lowpass, periodicity=",periodicity,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-200

# Estimation based on MDFA-MSE wrapper
mdfa_obj_noise_mse<-MDFA_mse(L,weight_func_noise,Lag,Gamma)$mdfa_obj 

b<-mdfa_obj_noise_mse$b
plot_estimate_func(mdfa_obj_noise_mse,weight_func_noise,Gamma)



# Comments
#   1. The filter coefficients in the first plot are one-sided (causal filter)
#   2. Since the filter does not look into the future, the output will be delayed (with respect to the non-causal target)
#     The corresponding delay/lag can be seen in the second plot (time-shift)
#   3. The fit of the target (violet line last plot) by the (DFA-) amplitude function can be 
#     seen in the last plot
#   4. Amplitude and shift differ from target: 
#     MSE-criterion in DFA makes an optimal 'mixed-fit' of both functions


# Compare output of (non-causal) ideal filter and DFA real-time (MSE solution) 

len<-1000
set.seed(1)
x<-rnorm(len)
# Length of filter (true ideal filter is bi-infinite)
M<-100
# This function computes the filter coefficients and applies the filter to x
id_obj<-ideal_filter_func(periodicity,M,x)
output_ideal<-id_obj$y

# DFA output
filt_obj<-filt_func(x,b)
  
output_dfa<-filt_obj$yhat
output_dfa[1:(L-1)]<-NA

# The filtered series (red line in plot below) is smooth: high-frequency noise has been damped
ts.plot(output_ideal,col="blue",main="Output of ideal lowpass (blue) vs DFA (red)")
lines(output_dfa,col="red")

# Comments: output of ideal filter (blue) is not available towards sample-end: DFA computes a 'best' (MSE) one-sided filter which runs till sample-end

# Let's now zoom into the plot and have a closer look at both series
anf<-400
enf<-500
ts.plot(output_ideal[anf:enf],col="blue",main="Output of ideal lowpass (blue) vs DFA (red)")
lines(output_dfa[anf:enf],col="red")
abline(h=0)

# We can see that DFA-output (red) is noisier (amplitude does not vanish in stopband) and slightly shifted to the right (time-shift is not zero)
#   Both effects (noise leakage and delay) can be explained/understood by looking at amplitude and time-shift functions
#   Beyond MSE (customization): improve noise suppression and reduce lag

# Play with parameters
# 1.  one can change periodicity
# 2.  one can modify Lag: 
#       Lag=0: nowcast
#       Lag<0: forecast (this is not a forecast of the original data but of the filtered target output)
#       Lag>0: backcast (in contrast to forecasting, backcasting of ideal lowpass is not trivial)
#               Backcasting is used when revising historical data (typically macro-economic data)

# Idea: user can specify a target which matches his particular interests (more flexible than classic time series approaches)
#   DFA derives an 'optimal' (here: MSE) real-time estimate (output is available towards sample-end)
#   The estimate is 'as close as possible' (mean-square sense) to (output of) target

#---------------------------------------------------------------------------
# Example 6: same as above but AR(1) series instead of noise


K<-600

# Spectrum ARMA: try various processes
a1<-0.9
b1<-NULL
plot_T<-T
# This function computes the spectrum of the ARMA-process
spec<-arma_spectrum_func(a1,b1,K,plot_T)$arma_spec
# Fill into weight_func: target (first column) and explanatory (second column); both are identical for univariate problems
weight_func<-cbind(spec,spec)
colnames(weight_func)<-c("target","explanatory")
weight_func_ar1<-weight_func

periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
plot(Gamma,type="l",main=paste("Ideal lowpass, periodicity=",periodicity,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-200

# Estimation based on MDFA-MSE wrapper
mdfa_obj_ar1_mse<-MDFA_mse(L,weight_func_ar1,Lag,Gamma)$mdfa_obj 

b<-mdfa_obj_ar1_mse$b
plot_estimate_func(mdfa_obj_ar1_mse,weight_func,Gamma)

# Note that we altered only the spectrum: AR(1) instead of white noise
# We now compare the previous amplitude function (example 5,white noise) with the new one (ar(1)) and  and interpret the result

plot_compare_two_DFA_designs(mdfa_obj_mse,mdfa_obj_ar1_mse,weight_func_noise,weight_func_ar1,Gamma)
  
# Interpretation
#   In the top-plot (white noise) the spectrum is flat: each frequency is equally important in DFA-criterion (see Wildi/McElroy)
#   In the bottom-plot (ar(1)) the spectrum is larger towards the lower frequencies (assuming a1>0).
#     The lower frequencies dominate in an AR(1)-process with a1>0
#     Therefore the match of the target (violet) and DFA-amplitude (black) is better towards the lower ferquencies for the AR(1)-process (bottom plot).
#     In contrast, the fit towards the higher frequencies is poorer for the AR(1) (when compared to white noise)

# Idea: the spectrum modulates the quality of the fit:
#  DFA-amplitude (time-shift) are matching the target better at those frequencies which are more heavily loaded (larger spectrum)

# Let's illustrate the previous idea (DFA MSE-criterion) by altering deliberately the spectrum
#---------------------------------------------------------------------------------------
# Example 7: same as above but with altered (tweaked) spectrum


# Let's modify the spectrum in frequency pi/2

weight_func_tweaked<-weight_func_ar1
weight_func_tweaked[nrow(weight_func_tweaked)/2,]<-100

# AR(1) spectrum with a huge spike at pi/2
plot(abs(weight_func_tweaked[,1]),type="l",main=paste("Tweaked ar(1)-spectrum, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Comments
#   1. We artificially enlarge the spectrum at pi/2
#   2. We therefore expect the fit of the amplitude function to improve towards pi/2
# Let's verify this conjecture

periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
plot(Gamma,type="l",main=paste("Ideal lowpass, periodicity=",periodicity,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-200

# Estimation based on MDFA-MSE wrapper
mdfa_obj_tweaked_mse<-MDFA_mse(L,weight_func_tweaked,Lag,Gamma)$mdfa_obj 

plot_estimate_func(mdfa_obj_tweaked_mse,weight_func_tweaked,Gamma)

# Let's compare the original AR(1) with the tweaked one
plot_compare_two_DFA_designs(mdfa_obj_tweaked_mse,mdfa_obj_ar1_mse,weight_func_tweaked,weight_func_ar1,Gamma)


# What happened?
#   Both amplitude functions are quite similar except that the tweaked one approaches zero at frequenvy pi/2
#   This is because 
#     1. the target is zero at frequency pi (this component will be eliminated by the ideal lowpass)
#     2. the tweaked spectrum is very large: thus DFA tries to damp that component more effectively (the amplitude approaches zero)
#  The fit of the time-shift improves similarly, see next example
#---------------------------------------------------------------------------------------
# Example 8: same as example 6 (AR(1) but with altered (tweaked) spectrum


# Let's modify the spectrum at a (arbitrary) low frequency
weight_func_tweaked<-weight_func_ar1
weight_func_tweaked[10,]<-100

# AR(1) spectrum with a huge spike at pi*(10-1)/600
plot(abs(weight_func_tweaked[,1]),type="l",main=paste("Tweaked ar(1)-spectrum, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Comments
#   1. We artificially enlarge the spectrum at pi/2
#   2. We therefore expect the fit of the amplitude function to improve towards pi/2
# Let's verify this conjecture

periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
plot(Gamma,type="l",main=paste("Ideal lowpass, periodicity=",periodicity,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length
L<-200

# Estimation based on MDFA-MSE wrapper
mdfa_obj_tweaked_mse<-MDFA_mse(L,weight_func_tweaked,Lag,Gamma)$mdfa_obj 

plot_estimate_func(mdfa_obj_tweaked_mse,weight_func_tweaked,Gamma)

# More generally: the fit of target by amplitude and time-shift functions improves at frequencies more heavily loalded by spectrum
#   DFA: weighted approximation whereby weight is supplied by spectrum
# DFA-MSE criterion makes a best possible compromise of time-shift (delay) and of amplitude (noise leakage) fitting

# Target is arbitrary: 
#   one/multi-step ahead forecasting, 
#   lowpass: nowcast, forecast backcast
#   bandpass, highpass, aribtrary target

# Spectrum is arbitrary
#   model-based (ARMA)
#   non-parametric (dft)
#   other...

# Customization (see later tutorial)
#   Control noise-leakage (reliability of signals)
#   Control time-shift (faster trading execution)

#---------------------------------------------------------------------------------------
# Example 9: same as previous example 8 (AR(1) with altered/tweaked) spectrum but smaller L
# Purpose: introduction to overfitting


# Same tweaked AR(1) spectrum as in exercise 8
weight_func_tweaked<-weight_func_ar1
weight_func_tweaked[10,]<-100

# AR(1) spectrum with a huge spike at pi*(10-1)/600
plot(abs(weight_func_tweaked[,1]),type="l",main=paste("Tweaked ar(1)-spectrum, denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Comments
#   1. In the previous exercise we observed that the amplitude and the time-shift functions adapted to
#     nearly perfect fits of the target at the spike-frequency
#   2. In most real-world applications this 'perfect fit' would be termed as overfitting
#   3. We could overfit the target (as weighted by the spectrum) because our filter was sufficiently flexible (large L)
#   4. What happens if we restrict the flexibility (degrees of freedom) by selecting a smaller L (filter length)?
# Let's do and see

periodicity<-5
cutoff<-pi/periodicity
Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
plot(Gamma,type="l",main=paste("Ideal lowpass, periodicity=",periodicity,", denseness=",K,sep=""),
     axes=F,xlab="Frequency",ylab="Amplitude",col="black")
# We take 2-nd colname from weight_func because the first column is the target        
mtext("Target",line=-1,col="black")
axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                  "4pi/6","5pi/6","pi"))
axis(2)
box()

# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length: much smaller than previous exercises
L<-2*K+1

# Estimation based on MDFA-MSE wrapper
mdfa_obj_tweaked_mse<-MDFA_mse(L,weight_func_tweaked,Lag,Gamma)$mdfa_obj 

plot_estimate_func(mdfa_obj_tweaked_mse,weight_func_tweaked,Gamma)


# Comments: assuming L<-10 has been selected
#   In contrast to previous exercise 8 the fit of target by amplitude and by shift at spectral spike is 'less perfect'
#   The smaller L=10 (less degrees of freedom) does not allow the filter to match arbitrary frequencies arbitrary well
# Additional/alternative experiments
#   1. What happens for very large L (say L<-600 or L=2*K); what happens if L>2*K
#     Hint: if L=2*K then the fit is perfect at all frequencies
#       The amplitude matches the target perfectly (both curves are indistinguishable)
#       The shift is zero in the passband
#       This is a very bad example of extreme overfitting!!!!!
#       For L>2*K the system is singular (the numerical optimization breaks down): there are more parameters available than system-equations
#   2. For L=10: look what happens when selecting a much larger peak-value: set weight_func_tweaked[10,]<-10000
#     Hint: the degrees of freedom are modestly-sized so the filter cannot overfit
#       For increasing (arbitrarily large) peak-values of the spectrum the amplitude and the shift match the target arbitrarily well towards the peak-frequency
#       This 'perfect fit' at the peak frequency induces an increasing mismatch (of the DFA-MSE) at all other frequencies
#       The optimization has to make a trade-off i.e. it cannot improve the fit at the peak-frequency without loosing at other frequencies


#-------------------------------------------------------------------------
# Wrap-up: what did we learn? 
# 1. User-interface:
#   The user-interface of the MSE-wrapper is made of L,weight_func,Lag,Gamma
#   -Gamma is a generic target: we learned allpass (classic forecasting), lowpass, bandpass, arbitrary shapes
#     The user is supposed to provide a target that matches his research interests
#   -weight_func is the weighting function modulating the fit of the target
#     weight_func has typically the meaning of a spectrum
#     one could rely on parametric (model-based) or non-parametric (dft) spectrum estimates
#     the user can tweak weight_func in order to obtain a filter which outperforms at particular frequencies
#   -L is the filter length and as such it accounts for the degrees of freedom: 
#     if L/K is 'small', then overfitting is contained/mitigated
#     very strong overfitting is obtained for L=2*K: in this case the target is matched perfectly by the DFA-MSE solution at the discrete frequency-ordinates (but the fit is awful at all other frequencies in [0,pi])
#   -Lag: 
#     allows for forecasting (Lag<0), nowcasting (Lag=0) or backcasting (Lag>0) of the target  
#     if Gamma is an allpass, then Lag<0 means ordinary (MSE-) forecasting
#     if Gamma is an ideal lowpass, then 
#       -Lag<0 means a forecast of the output of the ideal lowpass (not a forecast of the original series)
#       -Lag=0 means a nowcast of the output of the ideal lowpass (this is not a trivial estimation problem, quite the contrary)
#       -Lag>0 means a backcast of the output of the ideal lowpass (typically used when revising historical estimates: in various applications this task is called 'smoothing' in contrast to 'filtering' i.e. nowcasting: Kalman-filter is the nowcast, Kalman-smoother is the backcast)
# 2. DFA-MSE optimization criterion 
#   Given L, the MSE-criterion tries to achieve an optimal mix of amplitude and time-shift fitting (of the target) 
#   The amplitude is related to noise leakage (more 'wrong' signals)
#   The shift is related to the delay of the real-time filter (relative to the target)
#   Ideally, the user would prefer to have no leakage and no delay
#   Unfortunately, both requirements cannot be met simultaneously
#     In general (though not always: see customization) a stronger noise suppression implies a larger delay (conversely: a smaller delay generally results in greater leakage) 
#   However, in contrast to classic time series approaches, DFA allows for a customization of the optimization criterion
#     See trilemma paper McElroy/Wildi
#     See later tutorial
# 3. Overfitting is potentially obtained when
#   1. the spectrum is noisy/spiky and/or
#   2. the target is noisy/spiky and
#   3. L is large
#   If the weighting-function (the spectrum) is noisy and L is large then we need additional 'constraints'
#     Later tutorial about regularization 
