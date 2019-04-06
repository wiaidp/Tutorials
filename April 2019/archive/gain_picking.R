# Main gain-pick function
# Computes performances with/without gain picking per series and aggregate
gain_pick_main_func<-function(weekday,trade_all_pairs_with_weekday.obj,trade_all_pairs_no_weekday.obj,naked_gain_picking,FX_mat,start_date,threshold,lag_signal,plot_T,sign_rule)
{  
# Ensure causality
  lag_signal<-max(0,lag_signal)
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
  
  # Here one can plug in the original (log-transformed) FX_data
  #   This corresponds to the naked gain-picking based on ordering the data according to lengths of trades (as defined by yhat_mat)
  #   One can compare with perf_sign<-trade_all_pairs_no_weekday.obj$perf_sign (additional effect by applying trading rule without weekday effect)
  #   One can compare with perf_sign<-trade_all_pairs_with_weekday.obj$perf_sign (applying trading and weekday effect)
  
  if (naked_gain_picking)
    perf_sign<-perf_weight<-log(FX_mat[start_date])
  
  # Here we compute mean absolute daily movements: may be indicative for threshold
  apply(na.omit(abs(diff(perf_sign))),2,mean)
  apply(na.omit(abs(FX_diff)),2,mean)
  # We store the matrix performances in a vector of lists (because matrices have different dimensions depending on currency pair)
  perf_no_pick_list<-perf_pick_list<-date_list<-vector(mode="list")
  for (j in 1:ncol(perf_sign))#j<-1
  {
    yhat<-yhat_mat[,j]
    # Take differences
    if (sign_rule)
    {  
      perf<-perf_sign[,j]
    } else
    {
      perf<-perf_weight[,j]
    }
    gain_pick_obj<-gain_pick_func(threshold,yhat,perf,plot_T,lag_signal)
    if (j==1)
    {
      # The rows of the following matrix corresponds to sign-changes of signal (accounting for lag_signal above)
      #   The row are the cumulated performances from day 1 up to the last day of a connected trade (the number of columns corresponds to the longest trade)
      #     These are not the returns, but the cumulated perfs starting at zero at the begin of each connected trading period      
      perf_no_pick_list[[1]]<-gain_pick_obj$trade_mat_cumsum
      # As above but if a gain is picked, then the rows are filled with the gain up t the end of the trade    
      perf_pick_list[[1]]<-gain_pick_obj$trade_mat_pick
      perf_pick_mat<-gain_pick_obj$perf_pick_vec
      exposure_vec<-gain_pick_obj$market_exposure
      date_list[[1]]<-gain_pick_obj$date_mat
    } else
    {
      perf_no_pick_list[[j]]<-gain_pick_obj$trade_mat_cumsum
      perf_pick_list[[j]]<-gain_pick_obj$trade_mat_pick
      perf_pick_mat<-cbind(perf_pick_mat,gain_pick_obj$perf_pick_vec)
      exposure_vec<-c(exposure_vec,gain_pick_obj$market_exposure)
      date_list[[j]]<-gain_pick_obj$date_mat
    }
  }
  par(mfrow=c(1,2))
  perf_no_pick_agg<-as.xts(apply(perf_sign,1,mean))
  sharpe_no_pick_agg<-sqrt(250)*mean(diff(perf_no_pick_agg),na.rm=T)/sqrt(var(diff(perf_no_pick_agg),na.rm=T))
  mean_exposure<-100
  print(plot(perf_no_pick_agg,main=paste("Aggregate no pick, weekday=",weekday,", lag_signal=",lag_signal,", sharpe=",round(sharpe_no_pick_agg,3),", exposure=",mean_exposure,"%",sep="")))
  perf_pick_agg<-as.xts(apply(perf_pick_mat,1,mean))
  sharpe_pick_agg<-sqrt(250)*mean(diff(perf_pick_agg),na.rm=T)/sqrt(var(diff(perf_pick_agg),na.rm=T))
  mean_exposure<-round(mean(exposure_vec),3)*100
  print(plot(perf_pick_agg,main=paste("Aggregate with pick, weekday=",weekday,", lag_signal=",lag_signal,", sharpe=",round(sharpe_pick_agg,3),", exposure=",mean_exposure,"%",sep="")))
  
  return(list(perf_no_pick_list=perf_no_pick_list,perf_pick_list=perf_pick_list,perf_pick_mat=perf_pick_mat,exposure_vec=exposure_vec,date_list=date_list,perf_no_pick_agg=perf_no_pick_agg,sharpe_no_pick_agg=sharpe_no_pick_agg,perf_pick_agg=perf_pick_agg,sharpe_pick_agg=sharpe_pick_agg))
}






# These function applies gain-picking
# The natural events for applying gain picking are 'trades' or periods between sign changes of signal
# The threshold is positive and slightly below the mean absolute daily movements (say 20 bps but 0 should be fine too because daily moves are typically larger than spread; in contrast to intraday movements which require threshold>20bps, typically)
gain_pick_func<-function(threshold,yhat,perf,plot_T,lag_signal)
{
# Align filter outputs, performances and diff-performances, remove NAs  
# Some care is needed to understand what happens here
#   1. The performance in perf is based on lagged signal: the parameter lag_signal has no effect on this!!!!!
#   2. The effect of the parameter lag_signal is on the partitioning of the vector perf into time episodes defined by signal
#   3. If lag_signal=0 then perf is partitioned into episodes starting at sign change and ending at next sign change   
#       Note that the signal is already lagged by one time point when defining partitioning: which(sign(y_hat_al*lag(y_hat_al))==-1)
#       Therefore, lag_signal=0 is feasable (causality)
#   4. If lag_signal=-1 then perf would be partitioned into episodes starting one time point prior to sign change: this is not feasable
#   5. If lag_signal=1 then one waits one additional time point after sign change for starting the gain-picking sequence: this is feasable, of course
#   6. If performances are better if lag_signal=1 (than lag_signal=0 which is also feasable) then this means that 
#     first time point after sign change is not optimal for starting the gain-picking sequence
  
# Ensure causality
  if (lag_signal<0)
    print("Non causal setting: lag_signal is set to 0")
  lag_signal<-max(0,lag_signal)
  
  mat<-na.exclude(cbind(lag(yhat,lag_signal),(perf),diff(perf)))
  y_hat_al<-mat[,1]
  perf_al<-mat[,2]
  perf_al_diff<-mat[,3]
  length(y_hat_al)
  length(perf_al)
  length(perf_al_diff)
# These are the time points at which a sign of change is observed i.e. sign today is different from sign yesterday
  trade_vechh<-which(sign(y_hat_al*lag(y_hat_al))==-1)
# Add the last data-point if it's not already a trade: so the loop below will go up to the end of the data
  if (trade_vechh[length(trade_vechh)]<length(y_hat_al))
  {
    trade_vech<-c(trade_vechh,length(y_hat_al))
  } else
  {
    trade_vech<-trade_vechh
  }
# Add first non-NA data point
  trade_vec<-c(which(!is.na(y_hat_al))[1],trade_vech)               
# Maximum length of a trade: this determines the number of columns of the trade matrix
  number_col<-max(diff(trade_vec))
# We organize the data according to trades i.e. each time a new trade starts we assign the data to a new row of trade_mat until the trade ends (a new sign-change occurs)
#   The matrix date_mat stores the dates (this information is lost when slicing performances into trading episodes)  
  trade_mat_cumsum<-date_mat<-rep(NA,number_col)
# Loop over all trades
  for (j in 2:length(trade_vec))#j<-2
  {
# trade_row is the row of performances corresponding to the j-th trade  
#   It ends when one observes a change of sign i.e. perf_al[(trade_vec[j-1]+1):(trade_vec[j])]  
    trade_row_cumsum<-date_row<-rep(NA,number_col)
    trade_row_cumsum[1:(trade_vec[j]-trade_vec[j-1])]<-perf_al[(trade_vec[j-1]+1):(trade_vec[j])]-as.double(perf_al[trade_vec[j-1]])
    date_row[1:(trade_vec[j]-trade_vec[j-1])]<-substr(index(perf_al[(trade_vec[j-1]+1):(trade_vec[j])]),1,10)
    if (j==2)
    {
      trade_mat_cumsum<-trade_row_cumsum
      date_mat<-date_row
    } else
    {
      trade_mat_cumsum<-rbind(trade_mat_cumsum,trade_row_cumsum)
      date_mat<-rbind(date_mat,date_row)
    }
  }
# The following matrix trade_mat_pick fills the rows with the last value of perf_pick: this is useful when 
#   computing performances as a function of wait-time (wait-time as measured from the start of a trade)  
  trade_mat_pick<-trade_mat_cumsum
  
# Here we implement gain-picking withing trading episodes (between consecutive sign-changes of signal)
  for (i in 1:nrow(trade_mat_cumsum))#i<-2
  {
# Location of last observation in i-th trade: we llok for the last non-na (instead of looking for first NA which could a problem in the presence of NAs)
    non_na<-which(!is.na(trade_mat_cumsum[i,]))
    last_obs<-non_na[length(non_na)]
# Gain-picking or not: we use the cumulated performances (along rows) for determining gain-pick time
    if (max(trade_mat_cumsum[i,],na.rm=T)>threshold)#trade_mat[i,]
    {
# Location of first threshold trespassing
      max_loc<-which(trade_mat_cumsum[i,]>threshold)[1]
# Perf up to gain-pick:     
      perf_pick<-trade_mat_cumsum[i,1:max_loc]
# Fill with gain-picking value (not threshold because we have daily data i.e. we cannot trade at the time of threshold passing)
      if (max_loc<last_obs)
        perf_pick<-c(perf_pick,rep(trade_mat_cumsum[i,max_loc],last_obs-max_loc))
# If a gain was picked we fill the row in trade_mat_pick accordingly      
      trade_mat_pick[i,1:length(perf_pick)]<-perf_pick
    } else
    {
# No gain-picking: just fill in observed performance up to end of trade    
      perf_pick<-trade_mat_cumsum[i,1:last_obs]
    }
    if (i==1)
    {
# We add a zero for the first observation which is lost because of differences      
      perf_pick_vec<-c(0,perf_pick)
    } else
    {
      perf_pick_vec<-c(perf_pick_vec,perf_pick_vec[length(perf_pick_vec)]+perf_pick)
    }
  }
  names(perf_pick_vec)<-index(perf_al)
  perf_pick_vec<-as.xts(perf_pick_vec)
  length(perf_pick_vec)
  length(perf_al)
# Determine how often one is exposed to market  
  market_exposure<-1-length(which(diff(perf_pick_vec)==0))/length(perf_pick_vec)
  if (plot_T)
  {
    par(mfrow=c(2,1))
    sharpe_pick<-sqrt(250)*mean(diff(perf_pick_vec),na.rm=T)/sqrt(var(diff(perf_pick_vec),na.rm=T))
    plot(perf_pick_vec,main=paste(colnames(y_hat_al),": threshold: ",threshold,", Sharpe=",round(sharpe_pick,3),", exposure=",round(market_exposure,3)*100,"%",sep=""))
    sharpe<-sqrt(250)*mean(perf_al_diff,na.rm=T)/sqrt(var(perf_al_diff,na.rm=T))
    plot(perf_al,main=paste(colnames(y_hat_al)," no threshold, Sharpe=",round(sharpe,3),sep=""))
  }
  return(list(perf_pick_vec=perf_pick_vec,sharpe_pick=sharpe_pick,market_exposure=market_exposure,trade_mat_cumsum=trade_mat_cumsum,trade_mat_pick=trade_mat_pick,date_mat=date_mat)) 
}







#--------------------------------------------------------------------------------------------------------------
# Old code

# These function applies gain-picking
# The natural events for applying gain picking are 'trades' or periods between sign changes of signal
# The threshold is positive and slightly below the mean absolute daily movements (say 20 bps but 0 should be fine too because daily moves are typically larger than spread; in contrast to intraday movements which require threshold>20bps, typically)
perf_as_function_of_day_after_trade_func_old<-function(threshold,yhat,perf,plot_T,lag_signal)
{
  # Align filter outputs, performances and diff-performances, remove NAs  
  # Some care is needed to understand what happens here
  #   1. The performance in perf is based on lagged signal: the parameter lag_signal has no effect on this!!!!!
  #   2. The effect of the parameter lag_signal is on the partitioning of the vector perf into time episodes defined by signal
  #   3. If lag_signal=0 then perf is partitioned into episodes starting at sign change and ending at next sign change   
  #       Note that the signal is already lagged by one time point when defining partitioning: which(sign(y_hat_al*lag(y_hat_al))==-1)
  #       Therefore, lag_signal=0 is feasable (causality)
  #   4. If lag_signal=-1 then perf would be partitioned into episodes starting one time point prior to sign change: this is not feasable
  #   5. If lag_signal=1 then one waits one additional time point after sign change for starting the gain-picking sequence: this is feasable, of course
  #   6. If performances are better if lag_signal=1 (than lag_signal=0 which is also feasable) then this means that 
  #     first time point after sign change is not optimal for starting the gain-picking sequence
  
  # Ensure causality
  if (lag_signal<0)
    print("Non causal setting: lag_signal is set to 0")
  lag_signal<-max(0,lag_signal)
  
  mat<-na.exclude(cbind(lag(yhat,lag_signal),(perf),diff(perf)))
  y_hat_al<-mat[,1]
  perf_al<-mat[,2]
  perf_al_diff<-mat[,3]
  length(y_hat_al)
  length(perf_al)
  length(perf_al_diff)
  # These are the time points at which a sign of change is observed i.e. sign today is different from sign yesterday
  trade_vechh<-which(sign(y_hat_al*lag(y_hat_al))==-1)
  # Add the last data-point if it's not already a trade: so the loop below will go up to the end of the data
  if (trade_vechh[length(trade_vechh)]<length(y_hat_al))
  {
    trade_vech<-c(trade_vechh,length(y_hat_al))
  } else
  {
    trade_vech<-trade_vechh
  }
  # Add first non-NA data point
  trade_vec<-c(which(!is.na(y_hat_al))[1],trade_vech)               
  # Maximum length of a trade: this determines the number of columns of the trade matrix
  number_col<-max(diff(trade_vec))
  # We organize the data according to trades i.e. each time a new trade starts we assign the data to a new row of trade_mat
  trade_mat_cumsum<-rep(NA,number_col)
  # Loop over all trades
  for (j in 2:length(trade_vec))#j<-2
  {
    # trade_row is the row of performances corresponding to the j-th trade  
    #   It ends when one observes a change of sign i.e. perf_al[(trade_vec[j-1]+1):(trade_vec[j])]  
    trade_row_cumsum<-rep(NA,number_col)
    trade_row_cumsum[1:(trade_vec[j]-trade_vec[j-1])]<-perf_al[(trade_vec[j-1]+1):(trade_vec[j])]-as.double(perf_al[trade_vec[j-1]])
    if (j==2)
    {
      trade_mat_cumsum<-trade_row_cumsum
    } else
    {
      trade_mat_cumsum<-rbind(trade_mat_cumsum,trade_row_cumsum)
    }
  }
  # The following matrix trade_mat_pick fills the rows with the last value of perf_pick: this is useful when 
  #   computing performances as a function of wait-time (wait-time as measured from the start of a trade)  
  trade_mat_pick<-trade_mat_cumsum
  for (i in 1:nrow(trade_mat_cumsum))#i<-1
  {
    # Location of last observation in i-th trade: we llok for the last non-na (instead of looking for first NA which could a problem in the presence of NAs)
    non_na<-which(!is.na(trade_mat_cumsum[i,]))
    last_obs<-non_na[length(non_na)]
    # Gain-picking or not: we use the cumulated performances (along rows) for determining gain-pick time
    if (max(trade_mat_cumsum[i,],na.rm=T)>threshold)#trade_mat[i,]
    {
      # Location of first threshold trespassing
      max_loc<-which(trade_mat_cumsum[i,]>threshold)[1]
      # Perf up to gain-pick:     
      perf_pick<-trade_mat_cumsum[i,1:max_loc]
      # Fill with gain-picking value (not threshold because we have daily data i.e. we cannot trade at the time of threshold passing)
      if (max_loc<last_obs)
        perf_pick<-c(perf_pick,rep(trade_mat_cumsum[i,max_loc],last_obs-max_loc))
      # If a gain was picked we fill the row in trade_mat_pick accordingly      
      trade_mat_pick[i,1:length(perf_pick)]<-perf_pick
    } else
    {
      # No gain-picking: just fill in observed performance up to end of trade    
      perf_pick<-trade_mat_cumsum[i,1:last_obs]
    }
    if (i==1)
    {
      # We add a zero for the first observation which is lost because of differences      
      perf_pick_vec<-c(0,perf_pick)
    } else
    {
      perf_pick_vec<-c(perf_pick_vec,perf_pick_vec[length(perf_pick_vec)]+perf_pick)
    }
  }
  names(perf_pick_vec)<-index(perf_al)
  perf_pick_vec<-as.xts(perf_pick_vec)
  length(perf_pick_vec)
  length(perf_al)
  market_exposure<-1-length(which(diff(perf_pick_vec)==0))/length(perf_pick_vec)
  if (plot_T)
  {
    par(mfrow=c(2,1))
    sharpe_pick<-sqrt(250)*mean(diff(perf_pick_vec),na.rm=T)/sqrt(var(diff(perf_pick_vec),na.rm=T))
    plot(perf_pick_vec,main=paste(colnames(y_hat_al),": threshold: ",threshold,", Sharpe=",round(sharpe_pick,3),", exposure=",round(market_exposure,3)*100,"%",sep=""))
    sharpe<-sqrt(250)*mean(perf_al_diff,na.rm=T)/sqrt(var(perf_al_diff,na.rm=T))
    plot(perf_al,main=paste(colnames(y_hat_al)," no threshold, Sharpe=",round(sharpe,3),sep=""))
  }
  return(list(perf_pick_vec=perf_pick_vec,sharpe_pick=sharpe_pick,market_exposure=market_exposure,trade_mat_cumsum=trade_mat_cumsum,trade_mat_pick=trade_mat_pick)) 
}





