# Todos
#   Link to trilemma paper




# Purpose of tutorial: 
# -Illustrate ATS-trilemma (see McElroy/Wildi)
# -Analyze filter characteristics (amplitude/shift) when emphasizing Timeliness (the T in the ATS-trilemma)
# -Analyze filter characteristics (amplitude/shift) when emphasizing Smoothness (the S in the ATS-trilemma)
# -Analyze filter characteristics (amplitude/shift) when emphasizing both T and S: power of trilemma (when compared to classic MSE-dilemma)
# -Compare univariate customized DFA to bivariate leading indicator (MSE-) MDFA of previous tutorial
# -Compare all designs to customized bivariate MDFA


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
source("Common functions/functions_trilemma.r")
#---------------------------------------------------------------------------------------------
# Example 1: emphasizing Timeliness
# Univariate


# Example 2: emphasizing Smoothness
# Univariate



# Example 3: emphasizing both Timeliness and Smoothness
# Univariate





#--------------------------------------------------------------------------------
# Example 4: compare bivariate leading indicator and univariate customized
#   See previous tutorial for a background to the play_bivariate_func function below
#     This function sets-up a MDFA-experiment based on a bivariate noisy leading-indicator design
# Weak autocorrelation (typical for log-returns of (positive) economic data)
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





