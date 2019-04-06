


ma_cross_mdfa_mse_trade_func<-function(K,periodicity_short,periodicity_long,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)
{
#-----------------
# Derived settings
  cutoff_short<-pi/periodicity_short
# Target 
  Gamma_short<-(0:(K))<=K*cutoff_short/pi+1.e-9
#----------------
# Estimation based on MDFA-MSE wrapper
  mdfa_obj_mse_short<-MDFA_mse(L,weight_func,Lag,Gamma_short)$mdfa_obj 
  
  b_short<-mdfa_obj_mse_short$b
# Normalize coefficients: amplitude in zero is 1: this is needed for MA-cross to cancel unit-root of original data  
  b_short<-b_short/sum(b_short)
# Plot of amplitude
  if (plot_T)
  {
    plot(abs(mdfa_obj_mse_short$trffkt),type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
         axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
    axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
  
  
  cutoff_long<-pi/periodicity_long
  # Target 
  Gamma_long<-(0:(K))<=K*cutoff_long/pi+1.e-9
  #----------------
  # Estimation based on MDFA-MSE wrapper
  mdfa_obj_mse_long<-MDFA_mse(L,weight_func,Lag,Gamma_long)$mdfa_obj 
  
  b_long<-mdfa_obj_mse_long$b
# Normalize coefficients: amplitude in zero is 1: this is needed for MA-cross to cancel unit-root of original data  
  b_long<-b_long/sum(b_long)
  # Plot of amplitude
  if (plot_T)
  {
    plot(abs(mdfa_obj_mse_long$trffkt),type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
         axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
    axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
  
# Filter data
  
  yhat<-filt_func(x,b_short-b_long)$yhat
#-------------
# Compute performances: sign-rule and signal weighting
  
# Trading performance: sign-rule (don't forget to lag the signal by one day)
  cum_perf_sign<-cumsum(na.omit(lag(sign((yhat)),lag_fx)*(diff(x))))
  sharpe_sign<-sqrt(250)*mean(diff(cum_perf_sign),na.rm=T)/sqrt(var(diff(cum_perf_sign),na.rm=T))
# Signal weigthing
#   We weight the trade by the filter-output (proxi for strength of signal (plausible indicator for confidence))
#   Normalize trading rule by scale of signal: divide by sqrt(var(yhat,na.rm=T)) (otherwise scale of performance is not meaningful)
  cum_perf_weight<-cumsum(na.omit(lag((yhat),lag_fx)*diff(x)))/as.double(sqrt(var(yhat,na.rm=T)))
  sharpe_weight<-sqrt(250)*mean(diff(cum_perf_weight),na.rm=T)/sqrt(var(diff(cum_perf_weight),na.rm=T))
  
  # Plot
  if (plot_T)
  {
    par(mfrow=c(1,1))
    print(plot(cbind(diff(x),yhat),main="Returns (black) and filtered series (red)"))
    par(mfrow=c(2,1))
    print(plot(cum_perf_sign,main=paste("Signum rule, ",colnames(x),", sharpe: ",round(sharpe_sign,3),sep="")))
    print(plot(cum_perf_weight,main=paste("Signal weighting, ",colnames(x),", sharpe: ",round(sharpe_weight,3),sep="")))
    par(mfrow=c(1,1))
  }
  return(list(cum_perf_weight=cum_perf_weight,cum_perf_sign=cum_perf_sign,yhat=yhat,mdfa_obj_mse_short=mdfa_obj_mse_short,mdfa_obj_mse_long=mdfa_obj_mse_long,sharpe_sign=sharpe_sign,sharpe_weight=sharpe_weight))
}






# Filter function: applies a filter b to a series x which can be xts or double
#   If x is xts then time ordering of b is reversed
filt_func<-function(x,b)
{
  L<-length(b)
  yhat<-x
  if (is.matrix(x))
  {  
    length_time_series<-nrow(x)
  } else
  {
    if (is.vector(x))
    {
      length_time_series<-length(x)
    } else
    {
      print("Error: x is neither a matrix nor a vector!!!!")
    }
  }
  for (i in L:length_time_series)
  {
    # If x is an xts object then we cannot reorder x in desceding time i.e. x[i:(i-L+1)] is the same as  x[(i-L+1):i]
    #   Therefore, in this case, we have to revert the ordering of the b coefficients.    
    if (is.xts(x))
    {
      yhat[i]<-as.double(b[L:1]%*%x[i:(i-L+1)])#tail(x) x[(i-L+1):i]
    } else
    {
      yhat[i]<-as.double(b%*%x[i:(i-L+1)])#tail(x) x[(i-L+1):i]
    }
  }
  #  names(yhat)<-index(x)#index(yhat)  index(x)
  #  yhat<-as.xts(yhat,tz="GMT")
  return(list(yhat=yhat))
}





mdfa_mse_trade_func<-function(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)
{
  #-----------------
  # Derived settings
  cutoff<-pi/periodicity
  # Target 
  Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
  #----------------
  # Estimation based on MDFA-MSE wrapper
  mdfa_obj_mse<-MDFA_mse(L,weight_func,Lag,Gamma)$mdfa_obj 
  
  b<-mdfa_obj_mse$b
  # Plot of amplitude
  if (plot_T)
  {
    plot_estimate_func(b,mdfa_obj_mse,weight_func)
  }

  #-----------
  # Filtering
  if (ncol(b)==1)
  {  
    yhat<-filt_func(x,b)$yhat
  } else
  {
    yhat_mat<-x
    for (j in 1:ncol(b))#j<-1
    {
      yhat_mat[,j]<-filt_func(x[,j],b[,j])$yhat
    }
    yhat<-as.xts(apply(yhat_mat,1,mean))
  }
  # Weekday effects
  if (length(sign_day_of_week)==5)
  {
    for (i in 1:5)
    {
      yhat[.indexwday(yhat) %in% i]<-sign_day_of_week[i]*yhat[.indexwday(yhat) %in% i]
    }
  }
  #-------------
  # Compute performances: sign-rule and signal weighting
  # First series in x is always target (the data in x is re-ordered that way)
  
  # Trading performance: sign-rule (don't forget to lag the signal by one day)
  cum_perf_sign<-cumsum(na.omit(lag(sign((yhat)),lag_fx)*x[,1]))
  sharpe_sign<-sqrt(250)*mean(diff(cum_perf_sign),na.rm=T)/sqrt(var(diff(cum_perf_sign),na.rm=T))
  # Signal weigthing
  #   We weight the trade by the filter-output (proxi for strength of signal (plausible indicator for confidence))
  #   Normalize trading rule by scale of signal: divide by sqrt(var(yhat,na.rm=T)) (otherwise scale of performance is not meaningful)
  cum_perf_weight<-cumsum(na.omit(lag((yhat),lag_fx)*x[,1]))/as.double(sqrt(var(yhat,na.rm=T)))
  sharpe_weight<-sqrt(250)*mean(diff(cum_perf_weight),na.rm=T)/sqrt(var(diff(cum_perf_weight),na.rm=T))
  
  # Plot
  if (plot_T)
  {
    par(mfrow=c(1,1))
    print(plot(cbind(x[,1],yhat),main="Returns (black) and filtered series (red)"))
    par(mfrow=c(2,1))
    print(plot(cum_perf_sign,main=paste("Signum rule, ",colnames(x),", sharpe: ",round(sharpe_sign,3),sep="")))
    print(plot(cum_perf_weight,main=paste("Signal weighting, ",colnames(x),", sharpe: ",round(sharpe_weight,3),sep="")))
    par(mfrow=c(1,1))
  }
  return(list(cum_perf_weight=cum_perf_weight,cum_perf_sign=cum_perf_sign,yhat=yhat,mdfa_obj_mse=mdfa_obj_mse,sharpe_sign=sharpe_sign,sharpe_weight=sharpe_weight))
}






# This function computes performances when trading according to weekday-effect only (no filter involved)
trade_weekday_sign_func<-function(x_mat,sign_day_of_week)
{
# Signal: long only  
  yhat<-rep(1,nrow(x_mat))
  names(yhat)<-index(x_mat)
  yhat<-as.xts(yhat)
  for (j in 1:ncol(x_mat))
  {
  # Weekday effects
    if (length(sign_day_of_week)==5)
    {
      for (i in 1:5)
      {
        yhat[.indexwday(yhat) %in% i]<-sign_day_of_week[i]*yhat[.indexwday(yhat) %in% i]
      }
    }
    x<-x_mat[,j]
    cum_perf_sign<-cumsum(na.omit(lag(sign((yhat)),lag_fx)*x))
    sharpe_sign<-sqrt(250)*mean(diff(cum_perf_sign),na.rm=T)/sqrt(var(diff(cum_perf_sign),na.rm=T))
    if (j==1)
    {
      cum_perf_weekday_sign<-cum_perf_sign
      sharpe_weekday_sign<-sharpe_sign
    } else
    {
      cum_perf_weekday_sign<-cbind(cum_perf_weekday_sign,cum_perf_sign)
      sharpe_weekday_sign<-c(sharpe_weekday_sign,sharpe_sign)
    }
  }
  colnames(cum_perf_weekday_sign)<-colnames(x_mat)
  return(list(cum_perf_weekday_sign=cum_perf_weekday_sign,sharpe_weekday_sign=sharpe_weekday_sign))
}




# This function applies MDFA-MSE trading function to all pairs
#   It accounts for weekday-effects (sign_day_of_week)
trade_all_pairs_mse_func<-function(FX_diff,K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)
{
  
  # Loop over all G10 pairs
  for (i in 1:ncol(FX_diff))#i<-1
  {
    x<-FX_diff[,i]
    
    mdfa_mse_trade_obj<-mdfa_mse_trade_func(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,sign_day_of_week)
    
    if (i==1)
    {
      perf_sign<-mdfa_mse_trade_obj$cum_perf_sign#head(perf_sign)
      perf_weight<-mdfa_mse_trade_obj$cum_perf_weight
      yhat_mat<-mdfa_mse_trade_obj$yhat#head(yhat_mat,102)
      sharpe_sign<-mdfa_mse_trade_obj$sharpe_sign
      sharpe_weight<-mdfa_mse_trade_obj$sharpe_weight
    } else
    {
      perf_sign<-cbind(perf_sign,mdfa_mse_trade_obj$cum_perf_sign)
      perf_weight<-cbind(perf_weight,mdfa_mse_trade_obj$cum_perf_weight)
      yhat_mat<-cbind(yhat_mat,mdfa_mse_trade_obj$yhat)
      sharpe_sign<-c(sharpe_sign,mdfa_mse_trade_obj$sharpe_sign)
      sharpe_weight<-c(sharpe_weight,mdfa_mse_trade_obj$sharpe_weight)
    }
  }
  colnames(yhat_mat)<-colnames(perf_sign)<-colnames(perf_weight)<-colnames(FX_diff)
  
  return(list(yhat_mat=yhat_mat,perf_sign=perf_sign,perf_weight=perf_weight,sharpe_sign=sharpe_sign,sharpe_weight=sharpe_weight))  
  
}