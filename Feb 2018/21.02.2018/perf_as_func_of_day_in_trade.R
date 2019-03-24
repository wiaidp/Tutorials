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

main_perf_as_function_of_day_after_trade_func<-function(weekday,trade_all_pairs_with_weekday.obj,trade_all_pairs_no_weekday.obj,gain_pick_main_obj,FX,plot_T,stop_time)
{
  if (weekday)
  {
    perf_sign<-trade_all_pairs_with_weekday.obj$perf_sign
    perf_weight<-trade_all_pairs_with_weekday.obj$perf_weight
    # Note that yhat_mat changes sign depending on weekday effect which is not what we want when defining trading episodes
    #   depending on the length of the trades as specified by the filter output
    #   Therefore we rely on yhat_mat as provided by trade_all_pairs_no_weekday.obj$yhat_mat (without weekday effect)
    #yhat_mat<-trade_all_pairs_with_weekday.obj$yhat_mat
  } else
  {
    perf_sign<-trade_all_pairs_no_weekday.obj$perf_sign
    perf_weight<-trade_all_pairs_no_weekday.obj$perf_weight
    #  yhat_mat<-trade_all_pairs_no_weekday.obj$yhat_mat
  }
  # Note that yhat_mat changes sign depending on weekday effect which is not what we want when defining trading episodes
  #   depending on the length of the trades as specified by the filter output
  #   Therefore we rely on yhat_mat as provided by trade_all_pairs_no_weekday.obj$yhat_mat (without weekday effect)
  yhat_mat<-trade_all_pairs_no_weekday.obj$yhat_mat
  

  for (j in 1:ncol(perf_sign))#j<-2
  {
    # With gain-picking  
    perf_mat<-gain_pick_main_obj$perf_pick_list[[j]]
    date_mat<-gain_pick_main_obj$date_list[[j]]
    asset_name<-paste(colnames(FX)[j],": with gain-pick",sep="")
    
    par(mfrow=c(1,2))
    perf_as_function_of_day_after_trade_func.obj<-perf_as_function_of_day_after_trade_func(perf_mat,date_mat,plot_T,asset_name,stop_time)
    
    if (j==1)
    {
      perf_stop_mat_with_gain_pick<-perf_as_function_of_day_after_trade_func.obj$perf_stop_vec
      market_exposure_with_gain_picking<-perf_as_function_of_day_after_trade_func.obj$market_exposure
    } else
    {
      perf_stop_mat_with_gain_pick<-cbind(perf_stop_mat_with_gain_pick,perf_as_function_of_day_after_trade_func.obj$perf_stop_vec)
      market_exposure_with_gain_picking<-c(market_exposure_with_gain_picking,perf_as_function_of_day_after_trade_func.obj$market_exposure)
    }
    # Without gain-picking  
    perf_mat<-gain_pick_main_obj$perf_no_pick_list[[j]]
    asset_name<-paste(colnames(FX)[j],": without gain-pick",sep="")
    
    perf_as_function_of_day_after_trade_func.obj<-perf_as_function_of_day_after_trade_func(perf_mat,date_mat,plot_T,asset_name,stop_time)
    
    if (j==1)
    {
      perf_stop_mat_without_gain_pick<-perf_as_function_of_day_after_trade_func.obj$perf_stop_vec
      market_exposure_without_gain_picking<-perf_as_function_of_day_after_trade_func.obj$market_exposure
    } else
    {
      perf_stop_mat_without_gain_pick<-cbind(perf_stop_mat_without_gain_pick,perf_as_function_of_day_after_trade_func.obj$perf_stop_vec)
      market_exposure_without_gain_picking<-c(market_exposure_without_gain_picking,perf_as_function_of_day_after_trade_func.obj$market_exposure)
    }
  }
  names(market_exposure_without_gain_picking)<-names(market_exposure_with_gain_picking)<-colnames(perf_stop_mat_without_gain_pick)<-colnames(perf_stop_mat_with_gain_pick)<-colnames(FX)
  
  
  # Findings:
  #   Stop-rule without gain-picking provides slightly more performance at cost of higher market exposure and smaller sharpe
  # Trick: Impose full-weight (of trade) up to min(gain-pick,stop_time)
  #        Impose half-weight (of trade) from min(gain-pick,stop_time) to stop-time
  if (plot_T)
  {
    par(mfrow=c(1,2))
    perf_all_without<-as.xts(apply(perf_stop_mat_without_gain_pick,1,mean))
    sharpe_all_without<-sqrt(250)*mean(diff(perf_all_without),na.rm=T)/sqrt(var(diff(perf_all_without),na.rm=T))
    mean_market_exposure_without_gain_picking<-round(100*mean(market_exposure_without_gain_picking),1)
    print(plot(perf_all_without,main=paste("Aggregate without pick, weekday=",weekday,", lag_signal=",lag_signal,", stop-time=",stop_time,", sharpe=",round(sharpe_all_without,3),", exposure=",mean_market_exposure_without_gain_picking,"%",sep="")
    ))
    perf_all_with<-as.xts(apply(perf_stop_mat_with_gain_pick,1,mean))
    sharpe_all_with<-sqrt(250)*mean(diff(perf_all_with),na.rm=T)/sqrt(var(diff(perf_all_with),na.rm=T))
    mean_market_exposure_with_gain_picking<-round(100*mean(market_exposure_with_gain_picking),1)
    print(plot(perf_all_with,main=paste("Aggregate without pick, weekday=",weekday,", lag_signal=",lag_signal,", stop-time=",stop_time,", sharpe=",round(sharpe_all_with,3),", exposure=",mean_market_exposure_with_gain_picking,"%",sep="")
    ))
  }  
  
  
  # Here the comparison of all four strategies
  # 1. No gain-pick, no stop-time: top-left
  # 2. gain-pick, no stop-time: top-right
  # 3. No gain-pick, stop-time: bottom-left
  # 4. gain-pick, stop-time: bottom-right
  # Findings:
  # a. Stop-time effect roughly equal to gain-picking: compare bottom-left and top-right
  #     Slightly higher absolute return vs. slightly smaller sharpe vs. higher market exposure (but one could avoid some of the turmoils: Brexit, federal bank announcements)
  # b. Combination of gain-pick and stop-time (bottom right) 
  #     Best in terms of sharpe and market exposure (stop-time alone (bottom-left) hast slightly higher absolute returns)
  # c. The above findings with regards to relative performances of all four strategies remain true irrespective of weekday-effect
  #     Weekday=T doubles performances when starting in 2010 (absolute returns and sharpe)
  #     Weekday effect is no more relevant when starting in 2016
  #     Weekday effect is very relevant when starting in 2017 (doubles performances, once again)
  if (plot_T)
  {
    par(mfrow=c(2,2))
    perf_no_pick_agg<-as.xts(apply(perf_sign,1,mean))
    sharpe_no_pick_agg<-sqrt(250)*mean(diff(perf_no_pick_agg),na.rm=T)/sqrt(var(diff(perf_no_pick_agg),na.rm=T))
    mean_exposure<-100
    plot(perf_no_pick_agg,main=paste("Aggregate no pick, weekday=",weekday,", sign_rule=",sign_rule,", lag_signal=",lag_signal,", sharpe=",round(sharpe_no_pick_agg,3),", exposure=",mean_exposure,"%",sep=""))
    perf_pick_agg<-as.xts(apply(gain_pick_main_obj$perf_pick_mat,1,mean))
    sharpe_pick_agg<-sqrt(250)*mean(diff(perf_pick_agg),na.rm=T)/sqrt(var(diff(perf_pick_agg),na.rm=T))
    mean_exposure<-round(mean(gain_pick_main_obj$exposure_vec),3)*100
    plot(perf_pick_agg,main=paste("Aggregate with pick, weekday=",weekday,", sign_rule=",sign_rule,", lag_signal=",lag_signal,", sharpe=",round(sharpe_pick_agg,3),", exposure=",mean_exposure,"%",sep=""))
    perf_all_without<-as.xts(apply(perf_stop_mat_without_gain_pick,1,mean))
    sharpe_all_without<-sqrt(250)*mean(diff(perf_all_without),na.rm=T)/sqrt(var(diff(perf_all_without),na.rm=T))
    mean_market_exposure_without_gain_picking<-round(100*mean(market_exposure_without_gain_picking),1)
    plot(perf_all_without,main=paste("Aggregate without pick, weekday=",weekday,", sign_rule=",sign_rule,", lag_signal=",lag_signal,", stop-time=",stop_time,", sharpe=",round(sharpe_all_without,3),", exposure=",mean_market_exposure_without_gain_picking,"%",sep="")
    )
    perf_all_with<-as.xts(apply(perf_stop_mat_with_gain_pick,1,mean))
    sharpe_all_with<-sqrt(250)*mean(diff(perf_all_with),na.rm=T)/sqrt(var(diff(perf_all_with),na.rm=T))
    mean_market_exposure_with_gain_picking<-round(100*mean(market_exposure_with_gain_picking),1)
    plot(perf_all_with,main=paste("Aggregate without pick, weekday=",weekday,", sign_rule=",sign_rule,", lag_signal=",lag_signal,", stop-time=",stop_time,", sharpe=",round(sharpe_all_with,3),", exposure=",mean_market_exposure_with_gain_picking,"%",sep="")
    )
  }  
}










# Performances depending on day after signal change (first day, second day,...)
# The idea is: if performances are almost as good (or even better) when stoping trades earlier then we should do so
#   Note that if trade is stopped earlier (either by gain-pick or by next sign-change) then 
#   we just use the latest perf available i.e. we look at min(i,last day of k-th trade) where i runs through all possible days after a sign-change

# Here is the most interesting outcome when periodicity<-10 (the results are dependent on the signal, of course)
#   Without gain-picking the performance generally degrades with the age of the trade
#   With gain-picking this is much less the case because the gain-pick generally occurs early in the trade if the threshold is small
#   This result explains, also, why small thresholds should be used

perf_as_function_of_day_after_trade_func<-function(perf_mat,date_mat,plot_T,asset_name,stop_time)
{
  perf_mat_stop<-perf_mat_fill<-perf_mat  
# The variable k corresponds to sign changes (k-th sign change of signal: lag_signal above is accounted for, too)
  for (k in 1:nrow(perf_mat))#k<-1
  {
# Here we extend the rows by replacing NAs by the last observed value (the last value was obtained either as a gain pick or by the next sign-change)    
#   The i-th column of perf_mat then corresponds to performance at min(i,last day of k-th trade)
    if (length(which(is.na(perf_mat[k,])))>0)
    {
      first_na<-which(is.na(perf_mat[k,]))[1]
# Since perf_pick_mat has exactly the same structure (rows of identical lengths and filled parts of identical lengths) we
#   can use the same indices for the columns of perf_pick_mat (as for perf_no_pick_mat)      
      perf_mat_fill[k,first_na:ncol(perf_mat)]<-perf_mat[k,first_na-1]#perf_mat[k,]
      if (stop_time<first_na-1)
      {
# Perf up to stop_time        
        perf_stop<-perf_mat_stop[k,1:stop_time]
# Fill vector with stopped value up to end of trade i.e. first_na-1
#   This will give a vector of daily performances with stop-rule implemented
        perf_stop<-c(perf_stop,rep(perf_stop[length(perf_stop)],first_na-1-stop_time))
      } else
      {
        perf_stop<-perf_mat_stop[k,1:(first_na-1)]
      }
# Carry over the corresponding date      
      names(perf_stop)<-date_mat[k,1:(first_na-1)]
    } else
    {
# In this case the row is filled with observations (this is (one of) the longest connected trade) 
#   perf_mat does not have to be changed (it is filled up to the last column)      
      if (stop_time<ncol(perf_mat))
      {
# Perf up to stop_time        
        perf_stop<-perf_mat_stop[k,1:stop_time]
# Fill vector with stopped value up to end of trade i.e. first_na-1
        perf_stop<-c(perf_stop,rep(perf_stop[length(perf_stop)],ncol(perf_mat_stop)-stop_time))
      } else
      {
        perf_stop<-perf_mat[k,]
      }
      names(perf_stop)<-date_mat[k,]
    }
    if (k==1)
    {
      perf_stop_vec<-perf_stop
      names_perf_stop_vec<-names(perf_stop)
    } else
    {
# Don't forget to accumulate performances along the trading episodes
      perf_stop_vec<-c(perf_stop_vec,perf_stop_vec[length(perf_stop_vec)]+perf_stop)
      names_perf_stop_vec<-c(names_perf_stop_vec,names(perf_stop))
    }
  }
  names(perf_stop_vec)<-names_perf_stop_vec
  perf_stop_vec<-as.xts(perf_stop_vec)#tail(perf_stop_vec)
# Determine how often one is exposed to market  
  market_exposure<-1-length(which(diff(perf_stop_vec)==0))/length(perf_stop_vec)
  
# The i-th column of perf_mat corresponds to performance at min(i,last day of k-th trade)
# Cumsum is applied to obtain the cumulated performance
#   Note that the k-th row corresponds to the period between sign-change k and sign-change k+1
#   The length of the rows (dimension of matrix) corresponds to longest connected trade
#   Therefore, the cumulated performances (of the columns of perf_mat) are much shorter than the daily performances
#   The lengths (dimensions of perf_mat) depend on cutoff i.e. the signal specification  
  if (plot_T)
  {
    colo<-rainbow(ncol(perf_mat)*2)
    ts.plot(apply(perf_mat_fill,2,cumsum) ,col=colo,main=asset_name)
  }
  return(list(perf_stop_vec=perf_stop_vec,perf_mat_fill=perf_mat_fill,market_exposure=market_exposure))
    
}
