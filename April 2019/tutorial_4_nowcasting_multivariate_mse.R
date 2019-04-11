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
# Example 1 VAR(1)-process
#   MSE
#   Bivariate



lenh<-500
len<-121
# Specify AR-process
a1<-0.9
a1<-0.1
# Generate series for each AR(1)-process
set.seed(10)
x_long<-arima.sim(list(ar=a1),n=lenh)
x_short<-x_long[lenh/2+(-len/2):((len/2)-1)]
periodicity<-6
cutoff<-pi/periodicity
# Ideal filter: as benchmark
# Order of approximation
M<-100
y<-ideal_filter_func(periodicity,M,x_long)$y

set.seed(12)
# Scaling of the idiosyncratic noise
scale_idiosyncratic<-0.
eps<-rnorm(lenh)
indicator<-x_long+scale_idiosyncratic*eps
# Data: first column=target, second column=x, 
#   third column=shifted (leading) indicator
data_matrix<-na.exclude(cbind(x_long,x_long,c(indicator[2:lenh],indicator[lenh])))
dimnames(data_matrix)[[2]]<-c("target","x","leading indicator")
# Extract 120 observations from the long sample
data_matrix_120<-data_matrix[lenh/2+(-len/2):((len/2)-1),]
head(round(data_matrix_120,4))
# Spectrum
# One can use either the function per
weight_func<-cbind(per(data_matrix_120[,1],T)$DFT,per(data_matrix_120[,2],T)$DFT,per(data_matrix_120[,3],T)$DFT)
# Or one can use spec_comp (the latter is slightly more general than per)
#   Here we can supply the full data matrix data_matrix_120 at once: spec_comp then returns the multivariate dft 
weight_func<-spec_comp(nrow(data_matrix_120), data_matrix_120, 0)$weight_func

K<-nrow(weight_func)-1

Gamma<-(0:(K))<=K*cutoff/pi
# Nowcast (Lag=0), Backcast (Lag>0) and Forecast (Lag<0)
Lag<-0
# Filter length: L/K should be 'small'
L<-2*periodicity


#weight_func<-spec_comp(120,data_matrix_120,0)$weight_func

# Estimate filter coefficients
mdfa_obj<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 
# Filter coefficients
b_mat<-mdfa_obj$b
dimnames(b_mat)[[2]]<-c("x","leading indicator")
dimnames(b_mat)[[1]]<-paste("Lag ",0:(L-1),sep="")#dim(b_mat)
head(b_mat)





yhat_bivariate_leading_indicator<-filt_func(data_matrix[,2:ncol(data_matrix)],b_mat)$yhat


weight_func_univariate<-weight_func[,1:2]
mdfa_obj<-MDFA_mse(L,weight_func_univariate,Lag,Gamma)$mdfa_obj 
# Filter coefficients
b<-mdfa_obj$b

yhat_univariate<-filt_func(as.matrix(data_matrix[,2]),b)$yhat


length(yhat_bivariate_leading_indicator)
length(y)


y_target_leading_indicator<-y
perf_mse<-matrix(c(mean(na.exclude((yhat_bivariate_leading_indicator-
                                      y_target_leading_indicator))^2),
                   mean(na.exclude((yhat_univariate-
                                      y_target_leading_indicator))^2)),nrow=1)
dimnames(perf_mse)[[2]]<-c("bivariate MDFA","DFA")
dimnames(perf_mse)[[1]]<-"Sample MSE"
round(perf_mse,3)
# Compare MSE of bivariate (time-domain) with criterion value (frequency-domain)
round(mdfa_obj$MS_error,3)


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


# Zoom into plot
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














