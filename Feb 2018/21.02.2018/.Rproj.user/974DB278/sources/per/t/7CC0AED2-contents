rm(list=ls())

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)




# Load the data (if it's not in the MDFA-package)

load(file=paste(getwd(),"/FX_mat.Rdata",sep=""))

tail(FX_mat)

start_date<-"2010-01-01/"


#--------------------------------------------------------------------------------------------
# Data transformations:
#   Invert pairs not going from USD
#   Apply log (transform ratios into differences)
#   Map into [0,1]
#   Apply differences

# Here we invert pairs which go into USD i.e. pairs which are not referenced from USD
FX_mat_USD<-FX_mat

for (i in 1:ncol(FX_mat))
{
  if (substr(colnames(FX_mat)[i],1,3)!="USD")
  {
    FX_mat_USD[,i]<-1/FX_mat[,i]
    colnames(FX_mat_USD)[i]<-paste(substr(colnames(FX_mat)[i],4,6),substr(colnames(FX_mat)[i],1,3),sep="")
  }
}
# All pairs now go from USD (into G10 currencies)
tail(FX_mat_USD)

# Take the log: that way the ratios (USD/EUR) become differences (log(USD)-log(EUR))
#   take log prior to [0,1] mapping (otherwise log=-\infty)
FX<-log(FX_mat_USD[start_date])
# Compute log-returns for trading below (note that signals of inverted pairs must be inverted in real-world)
FX_diff<-diff(FX)

# Maping into [0,1]: for purpose of direct comparison of pairs
FX_scale<-FX

for (i in 1:ncol(FX))
  FX_scale[,i]<-(FX[,i]-as.double(min(FX[,i])))/(as.double(max(FX[,i]))-as.double(min(FX[,i])))

tail(FX_scale)

# Kind of cointegration from 2000 on
colo<-rainbow(ncol(FX_scale))
plot(FX_scale[,1],main=paste("USD against G9, scaled"))
mtext(colnames(FX_scale)[1],col=colo[1],line=-1)
for (i in 2:ncol(FX_scale))
{  
  lines(FX_scale[,i],col=colo[i])
  mtext(colnames(FX_scale)[i],col=colo[i],line=-i)
}
lines(as.xts(apply(FX_scale,1,mean)),lwd=2,col="black")

FX_scale_diff<-diff(FX_scale)

# By aggragting USD over all pairs we obtain USD against a 'world basket': USD/World
#   Any other pair subtracted from US/world will provide the world-basket of the corresponding currency
#   So for example log(USD/world)-log(USDJPY)=log(JPY/world)
USDg9<-as.xts(apply(FX_scale,1,mean))

# Here we compute the ratio of var(diff)/var(level) (mean reversion vs. trending)
sqrt(var(diff(USDg9),na.rm=T)/var(USDg9))
for (i in 1:ncol(FX_scale))
  print(c(colnames(FX_scale)[i],sqrt(var(diff(FX_scale[,i]),na.rm=T)/var(FX_scale[,i]))))


#-----------------------------------------------------------------------------------------------
# Filtering and trading: sign and weight rules
# We use FX_diff as computed above: all pairs go from USD (in any g9-currency), log-returns
#   Load function
source(paste(getwd(),"/mdfa_mse_trade_func.r",sep=""))
source(paste(getwd(),"/plot_func.r",sep=""))
source(paste(getwd(),"/perf_coincident_signals.r",sep=""))
source(paste(getwd(),"/gain_picking.r",sep=""))
source(paste(getwd(),"/perf_as_func_of_day_in_trade.r",sep=""))


K<-240
# Spectrum: white noise assumption
weight_func<-matrix(rep(1,2*(K+1)),ncol=2)
#weight_func[1,]<-1000
periodicity<-10
Lag<-0
L<-10*periodicity
L<-100
L<-min(2*K,L)
plot_T<-F
asset<-"USDEUR"
lag_fx<-1
sign_day_of_week<-rep(1,5)
sign_day_of_week[5]<--1
x<-FX_diff[,asset]
#-------------------------------------------------
# MSE estimation and trading: sign rule and signal weighting

x<-FX_diff[,asset]
# This function is in mdfa_mse_trade_func.r
# MSE estimation and trading: sign rule and signal weighting
mdfa_mse_trade_obj<-mdfa_mse_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)

cum_perf_sign<-mdfa_mse_trade_obj$cum_perf_sign
# Weekday effect
plot_weekday_func(cum_perf_sign,asset)

# This function is in mdfa_mse_trade_func.r
# Loop over all G10 pairs
#   With weekday effect
trade_all_pairs_with_weekday.obj<-trade_all_pairs_mse_func(FX_diff,K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)
#   No weekday effect (sign_day_of_week=1 for all days)
sign_all_days_one<-rep(1,ncol(FX_diff))
trade_all_pairs_no_weekday.obj<-trade_all_pairs_mse_func(FX_diff,K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_all_days_one)

# perf_plot function is in plot_func.r
# Performances: long only
sharpe_long<-sqrt(250)*apply(FX_diff[-1,],2,mean)/sqrt(apply(FX_diff[-1,],2,var))
name<-"Buy&hold"
perf_obj<-perf_plot(as.xts(t(t(FX)-as.double(FX[1,]))),sharpe_long,name)
# Performances: weekday only (no filtering)
name<-"Weekday only"
trade_weekday.obj<-trade_weekday_sign_func(FX_diff,sign_day_of_week)
perf_obj<-perf_plot(trade_weekday.obj$cum_perf_weekday_sign,trade_weekday.obj$sharpe_weekday_sign,name)
# Performances: filtering without weekday effect
name<-"Filtering only: sign "
perf_sign<-trade_all_pairs_no_weekday.obj$perf_sign
sharpe_sign<-trade_all_pairs_no_weekday.obj$sharpe_sign
perf_obj<-perf_plot(perf_sign,sharpe_sign,name)
# Performances: filtering with weekday effect
name<-"Filtering with Weekday: sign"
perf_sign<-trade_all_pairs_with_weekday.obj$perf_sign
sharpe_sign<-trade_all_pairs_with_weekday.obj$sharpe_sign
perf_obj<-perf_plot(perf_sign,sharpe_sign,name)
# Performances: filtering without weekday effect
name<-"Filtering only: weight "
perf_weight<-trade_all_pairs_no_weekday.obj$perf_weight
sharpe_weight<-trade_all_pairs_no_weekday.obj$sharpe_weight
perf_obj<-perf_plot(perf_weight,sharpe_weight,name)
# Performances: filtering with weekday effect
name<-"Filtering with Weekday: weight "
perf_weight<-trade_all_pairs_with_weekday.obj$perf_weight
sharpe_weight<-trade_all_pairs_with_weekday.obj$sharpe_weight
perf_obj<-perf_plot(perf_weight,sharpe_weight,name)

#------------------------------------------------------------------
# Function perf_coincident_signals is in experimental.r
# It selects time points with coincident filter signals (for example at last 4 positive signals with respect to USD)
# and computes performances when restricted to these time points: the performances as well as the market exposure 
# (how often one is in the market) are computed

# Findings: 
# For good signals (periodicity=10) one finds that 
#   1. at-least performances (at least 6 or 8 coincident signs) are still very good even with reduced market exposure (30%)
#   2. at-most performances degrade rapidly i.e. conflicting signals are bad times
# Recall that all pairs are going out of USD here i.e. USD is the base currency for all 9 pairs
perf_coincident_signals(trade_all_pairs_no_weekday.obj)
  
perf_coincident_signals(trade_all_pairs_with_weekday.obj)

#--------------------------------------------------------------------
# Function gain_pick_func is in gain_picking.r
# data: 
#   either sign-trading or weighted trading
#   with or without weekday effect

# Apply weekday effects as specified in sign_day_of_week: 
#   This has a strong effect: doubling performances, either on long or short current sample, with or without picking or stop-rules
weekday<-T
# Trading rules: apply sign of filter or whole filter output (weight): 
#   weight=T is interesting particularly towards sample end (2017/2018).
#   The rule is differentiable and one can formulate a corresponding MDFA-target
#   Trading costs: new trade each day in order to follow signal-weight (but that could be discretized)
#   Possibly huge effects at turmoils due to implicit squaring (trim signal in such a case)
sign_rule<-F
# If T then we omit filtering i.e. we can evaluate effect of gain-picking alone (without filter): 
#   This is just for evaluation purposes: what is the added value of filter
#   Effect is strong (at least when periodicity=10) i.e. filter is very important
#   default is F
naked_gain_picking<-F
# Standard gain-pick threshold: 20 bp (note that mean absolute daily move is about 50bp)
# One could select threshold=0: perf are even better (for daily trading the movements are 'large' in comparison to spread; in contrast to intraday trading where threshold>20bP is required)
# Selecting threshold=-100 means that the condition is always immediately satisfied i.e. this is the performance 
#   of first day after a sign change occurred
threshold<-0.002*periodicity/10
# This parameter affects the partitioning of perf (performances) into episodes specified by sign changes of signal
#   Note that perf is based on lagged signal, already. So lag_signal does not affect perf.
#   Instead, lag_signal affects gain-picking: starting at sign change (or later as specified by lag_signal) and ending at next sign change (or later as specified by lag_signal)
#   Any positive value is feasable (causality)
# Note that one uses the signal without weekday effects for perf with or without weekday
#   i.e. one uses the pure filter output without additional sign-changes depending on weekday
lag_signal<-1
plot_T<-T

gain_pick_main_obj<-gain_pick_main_func(weekday,trade_all_pairs_with_weekday.obj,trade_all_pairs_no_weekday.obj,naked_gain_picking,FX_mat,start_date,threshold,lag_signal,plot_T,sign_rule)
  

#----------------------------------------------------------------
# Performances depending on day after signal change (first day, second day,...)
# The idea is: if performances are almost as good (or even better) when stoping trades earlier then we should do so
#   Note that if trade is stopped earlier (either by gain-pick or by next sign-change) then 
#   we just use the latest perf available i.e. we look at min(i,last day of k-th trade) where i runs through all possible days after a sign-change

# Here is the most interesting outcome when periodicity<-10 (the results are dependent on the signal, of course)
#   Without gain-picking the performance generally degrades with the age of the trade
#   With gain-picking this is much less the case because the gain-pick generally occurs early in the trade if the threshold is small
#   This result explains, also, why small thresholds should be used
#   More interestingly, though, is the fact that the performance is maxed at stop-time<-periodicity/2+1:
#     This corresponds more or less to the time-shift of the filter i.e. this is an intuitively very appealing result in the presence of mean-reversion

stop_time<-as.integer(periodicity/2+1)
#stop_time<-6

main_perf_as_function_of_day_after_trade_func(weekday,trade_all_pairs_with_weekday.obj,trade_all_pairs_no_weekday.obj,gain_pick_main_obj,FX,plot_T,stop_time)




  

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
# Old code






if (F)
{
# Brief check for the above piece of code  
  tail(gain_pick_obj$perf_pick_vec)
  tail(apply(perf_pick_mat,2,cumsum))

  tail(perf)
  tail(apply(perf_no_pick_mat,2,cumsum))
}

tail(perf)


#-------------------------------------------------------------------
# Weekday effect on aggregate (equally-weighted) portfolio

perf_agg<-perf_obj$perf_agg  

# Weekday effect
plot_weekday_func(perf_agg,"Equally-weighted")






#-----------------------------------------------------

# Alternating signs
acf(na.omit(diff(perf_agg)))

new_perf<-cumsum(na.omit((-1)^(1:length(diff(perf_agg)))*diff(perf_agg)))
sharpe_new<-sqrt(250)*mean(diff(new_perf),na.rm=T)/sqrt(var(diff(new_perf),na.rm=T))
plot(new_perf,main=round(sharpe_new,3))
#acf(na.omit(diff(new_perf)))


# Trade FX-pair defined in asset by aggregate filter output
x<-FX_diff[,asset]
agg_filt<-apply(yhat_mat,1,mean)
names(agg_filt)<-index(yhat_mat)
cum_perf_sign<-cumsum(na.omit(lag(sign(as.xts(agg_filt)),lag_fx)*as.xts(x)))
sharpe_sign<-sqrt(250)*mean(diff(cum_perf_sign),na.rm=T)/sqrt(var(diff(cum_perf_sign),na.rm=T))
plot(cum_perf_sign,main=round(sharpe_sign,3))


