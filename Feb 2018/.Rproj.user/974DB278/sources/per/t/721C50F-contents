### R code from vignette source 'C:/wia_desktop/2018/Projekte/MDFA-Legacy/Sweave/Rnw/MDFA_Legacy.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: init
###################################################
#rm(list=ls())


###################################################
### code chunk number 2: init
###################################################
# Load packages: time series and xts
#library(tseries)
library(xts)
# Library for tables
library(Hmisc)
require(xtable)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
devtools::install_github("wiaidp/MDFA")
# MDFA package
library(MDFA)





path.pgm <- paste(getwd(),"/",sep="")


###################################################
### code chunk number 6: dft
###################################################
head(per,100)


###################################################
### code chunk number 7: dfa_ms
###################################################
# This function computes MSE DFA solutions 
# L is the length of the MA filter,
# periodogram is the frequency weighting function in the DFA
# Gamma is the transfer function of the symmetric filter (target) and
# Lag is the lag-parameter: Lag=0 implies real-time filtering, Lag=L/2
#     implies symmetric filter
# The function returns optimal coefficients as well as the transfer 
#     function of the optimized real-time filter
head(dfa_ms,100)


###################################################
### code chunk number 8: dfa_ms
###################################################
head(dfa_analytic)


###################################################
### code chunk number 9: dfa_ms
###################################################
set.seed(1)
len <- 100
target <- arima.sim(list(ar=0.9),n=len)
explanatory_2 <- target+rnorm(len)
explanatory <- cbind(target,explanatory_2)
x <- cbind(target,explanatory)
dimnames(x)[[2]] <- c("target","explanatory 1","explanatory 2")
head(x)


###################################################
### code chunk number 10: dfa_ms
###################################################
x<-cbind(x[,1],lag(x[,2:3],-1))
dimnames(x)[[2]]<-c("target","lagged explanatory 1","lagged explanatory 2")
head(x)


###################################################
### code chunk number 11: dfa_ms
###################################################
spec_comp


###################################################
### code chunk number 12: dfa_ms
###################################################
head(mdfa_analytic)


###################################################
### code chunk number 13: dfa_ms
###################################################
weight_func <- matrix(rep(1:6,2),ncol=2)
L <- 2


###################################################
### code chunk number 14: dfa_ms
###################################################
d<-0
lin_eta<-F
lambda<-0
Lag<-0
eta<-0
i1<-F
i2<-F
weight_constraint<-rep(1/(ncol(weight_func)-1),ncol(weight_func)-1)
lambda_cross<-lambda_smooth<-0
lambda_decay<-c(0,0)
lin_expweight<-F
shift_constraint<-rep(0,ncol(weight_func)-1)
grand_mean<-F
b0_H0<-NULL
c_eta<-F
weights_only<-F
weight_structure<-c(0,0)
white_noise<-F
synchronicity<-F
cutoff<-pi
lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)
troikaner<-F


###################################################
### code chunk number 15: dfa_ms
###################################################
source(file=paste(path.pgm,"control_default.r",sep=""))



###################################################
### code chunk number 17: dfa_ms
###################################################
head(MDFA_mse)
head(MDFA_mse_constraint)
head(MDFA_cust)
head(MDFA_cust_constraint)
head(MDFA_reg)
head(MDFA_reg_constraint)