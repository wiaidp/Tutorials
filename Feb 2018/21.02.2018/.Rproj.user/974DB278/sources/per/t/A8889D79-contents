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
    plot(abs(mdfa_obj_mse$trffkt),type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
         axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
    axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
# We need ordinary time series (not xts) for filtering because of lag-structure when computing output
  x_fil<-as.double(x)
  names(x_fil)<-index(x)
#-----------
# Filtering
  yhat<-rep(NA,length(x_fil))
  for (i in L:length(x))
    yhat[i]<-sum(b*x_fil[i:(i-L+1)])
  names(yhat)<-names(x_fil)
  yhat<-as.xts(yhat)
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
  
# Trading performance: sign-rule (don't forget to lag the signal by one day)
  cum_perf_sign<-cumsum(na.omit(lag(sign((yhat)),lag_fx)*as.xts(x)))
  sharpe_sign<-sqrt(250)*mean(diff(cum_perf_sign),na.rm=T)/sqrt(var(diff(cum_perf_sign),na.rm=T))
# Signal weigthing
#   We weight the trade by the filter-output (proxi for strength of signal (plausible indicator for confidence))
#   Normalize trading rule by scale of signal: divide by sqrt(var(yhat,na.rm=T)) (otherwise scale of performance is not meaningful)
  cum_perf_weight<-cumsum(na.omit(lag((yhat),lag_fx)*as.xts(x)))/as.double(sqrt(var(yhat,na.rm=T)))
  sharpe_weight<-sqrt(250)*mean(diff(cum_perf_weight),na.rm=T)/sqrt(var(diff(cum_perf_weight),na.rm=T))
  
  # Plot
  if (plot_T)
  {
  par(mfrow=c(1,1))
  plot(as.xts(x),main="Returns (black) and filtered series (red)")
  lines(as.xts(yhat),col="red")
    par(mfrow=c(2,1))
    plot(cum_perf_sign,main=paste("Signum rule, ",colnames(x),", sharpe: ",round(sharpe_sign,3),sep=""))
    plot(cum_perf_weight,main=paste("Signal weighting, ",colnames(x),", sharpe: ",round(sharpe_weight,3),sep=""))
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