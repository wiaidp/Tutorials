# This function accounts for number of coincident filter signals (signs of filter outputs)
#   It selects all time points with either at most k coincident filter outputs or at least k coincident filter outputs 
#     Recall that all pairs start from USD (so that the concept of coincident signs makes sense)
#   Performances are then computed
#   If performances are still good despite the restricted market exposure then this is good news

perf_coincident_signals<-function(trade_all_pairs.obj)
{

  yhat_mat<-trade_all_pairs.obj$yhat_mat

  for (jj in 1+0:(ncol(yhat_mat)/2))  
  {
    # How many signals must coincide
    jk<-2*jj
    name<-paste("At most ",jk,sep="")
# Select time points at which number of coincident signals is at most jk
    int_jk<-which(abs(apply(sign(yhat_mat),1,sum))<=jk)
    
    compute_perf_coincident_signals(yhat_mat,int_jk,name,trade_all_pairs.obj)
  }
  for (jj in 0:(ncol(yhat_mat)/2))  
  {
    # How many signals must coincide
    jk<-2*jj
    name<-paste("At least ",jk,sep="")
    # Select time points at which number of coincident signals is at least jk
    int_jk<-which(abs(apply(sign(yhat_mat),1,sum))>=jk)
    
    compute_perf_coincident_signals(yhat_mat,int_jk,name,trade_all_pairs.obj)
  }
  
}


# This function belongs to perf_coincident_signals above
compute_perf_coincident_signals<-function(yhat_mat,int_jk,name,trade_all_pairs.obj)
{
# Ratio of retained observations to total number of observations
  market_exposure<-length(int_jk)/nrow(na.exclude(yhat_mat))
# Time index (dates) of corresponding time points: 
#   lag by lag_s (since we apply the rule to tomorrow's trade)
  lag_s<-1
  index_jk<-index(yhat_mat)[int_jk+lag_s]
# Select the corresponding performances 
  trade_select<-diff(trade_all_pairs.obj$perf_sign)[index_jk,]
# Compute cumulated performances
  perf_select<-cumsum(na.exclude(as.xts(apply(trade_select,1,mean))))
# Plot
  print(plot(perf_select,main=paste(name," coincident signs. Market exposure: ",100*round(market_exposure,3),"%",sep="")))
  if (F)
    acf(na.exclude(diff(perf_select)))
}



# Old code: This function belongs to perf_coincident_signals above
compute_perf_coincident_signals_old<-function(trade_all_pairs_no_weekday.obj,sign_day_of_week,int_jk,name)
{
  # Ratio of retained observations to total number of observations
  market_exposure<-length(int_jk)/nrow(na.exclude(trade_all_pairs_no_weekday.obj$yhat_mat))
  # Time index (dates) of corresponding time points: 
  #   lag by lag_s (since we apply the rule to tomorrow's trade)
  lag_s<-1
  index_jk<-index(trade_all_pairs_no_weekday.obj$yhat_mat)[int_jk+lag_s]
  # Select the corresponding performances 
  trade_select_no_weekday<-diff(trade_all_pairs_no_weekday.obj$perf_sign[index_jk,])
  # Shift weekday signs by one day
  sign_day_of_week_shift<-c(sign_day_of_week[5],sign_day_of_week[1:4])
  # The operation does not work with xts so that we have to rely on matrices...
  trade_select_with_weekday_mat<-trade_select_no_weekday_mat<-as.matrix(trade_select_no_weekday)
  # Loop over each weeday and change sign if necessary
  for (i in 1:length(sign_day_of_week))
    trade_select_with_weekday_mat[.indexwday(trade_select_no_weekday) %in% i,]<-sign_day_of_week_shift[i]*trade_select_no_weekday_mat[.indexwday(trade_select_no_weekday) %in% i]
  # Transform matrix in xts back... (no comment)
  trade_select_with_weekday<-as.xts(trade_select_with_weekday_mat)
  
  
  # Select either one
  if (F)
  {
    trade_select<-trade_select_no_weekday
  } else
  {
    trade_select<-trade_select_with_weekday
  }
  
  # Compute cumulated performances
  perf_select<-cumsum(na.exclude(as.xts(apply(trade_select,1,mean))))
  # Plot
  plot(perf_select,main=paste(name," coincident signs. Market exposure: ",round(market_exposure,3),sep=""))
  if (F)
    acf(na.exclude(diff(perf_select)))
}
