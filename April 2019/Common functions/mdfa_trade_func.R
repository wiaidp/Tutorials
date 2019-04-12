



# Filter function: applies a filter b to a series x which can be xts or double
#   If x is xts then time ordering of b is reversed
filt_func<-function(x,b)
{
  L<-nrow(b)
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
  if (is.xts(x))
  {
    yhat<-x[,1]
  } else
  {
    yhat<-rep(NA,length_time_series)
  }
  for (i in L:length_time_series)#i<-L
  {
    # If x is an xts object then we cannot reorder x in desceding time i.e. x[i:(i-L+1)] is the same as  x[(i-L+1):i]
    #   Therefore, in this case, we have to revert the ordering of the b coefficients.    
    if (is.xts(x))
    {
      if (ncol(b)>1)
      {
        yhat[i]<-as.double(sum(apply(b[L:1,]*x[i:(i-L+1),],1,sum)))
      } else
      {
        yhat[i]<-as.double(b[L:1,]%*%x[i:(i-L+1)])#tail(x) x[(i-L+1):i]
      }
    } else
    {
      if (ncol(b)>1)
      {
        yhat[i]<-as.double(sum(apply(b[1:L,]*x[i:(i-L+1),],1,sum)))
      } else
      {
        yhat[i]<-as.double(as.vector(b)%*%x[i:(i-L+1)])#tail(x) x[(i-L+1):i]
      }
    }
  }
  #  names(yhat)<-index(x)#index(yhat)  index(x)
  #  yhat<-as.xts(yhat,tz="GMT")
  return(list(yhat=yhat))
}




#x<-data_filter
mdfa_reg_trade_func<-function(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func,lambda_cross,lambda_decay,lambda_smooth,lambda,eta)
{
  #-----------------
  # Derived settings
  cutoff<-pi/periodicity
  # Target 
  Gamma<-(0:(K))<=K*cutoff/pi+1.e-9
  # MSE
  #----------------
  # Estimation based on MDFA-Regularization wrapper
  
  mdfa_obj_mse_reg<-MDFA_reg(L,weight_func,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth)$mdfa_obj
  
  b<-mdfa_obj_mse_reg$b
  # Plot of amplitude
  if (plot_T)
  {
    plot_estimate_func(mdfa_obj_mse_reg,weight_func,Gamma)
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
  #-------------
  # Compute performances: sign-rule and signal weighting
  # First series in x is always target (the data in x is re-ordered that way)
  
  # Trading performance: sign-rule (don't forget to lag the signal by one day)
  cum_perf_sign<-cumsum(na.omit(lag(sign((yhat)),lag_fx)*x[,1]))
  sharpe_sign<-sqrt(250)*mean(diff(cum_perf_sign),na.rm=T)/sqrt(var(diff(cum_perf_sign),na.rm=T))

  # Plot
  if (plot_T)
  {
    par(mfrow=c(1,1))
    print(plot(cbind(x[,1],yhat),main="Returns (black) and filtered series (red)"))
    print(plot(cum_perf_sign,main=paste(colnames(x)[1],", sharpe: ",round(sharpe_sign,3),sep="")))
    par(mfrow=c(1,1))
  }
  return(list(cum_perf_sign=cum_perf_sign,yhat=yhat,mdfa_obj_mse_reg=mdfa_obj_mse_reg,sharpe_sign=sharpe_sign))
}





#x<-data_filt
mdfa_mse_reg_trade_func<-function(K,periodicity,L,Lag,lag_fx,x,plot_T,weight_func)
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
    plot_estimate_func(mdfa_obj_mse,weight_func,Gamma)
  }
  yhat<-filt_func(x,b)$yhat
  #-----------
  # Filtering
#  if (ncol(b)==1)
  {  
    yhat<-filt_func(x,b)$yhat
  } 
  #-------------
  # Compute performances: sign-rule and signal weighting
  # First series in x is always target (the data in x is re-ordered that way)
  # Trading performance: sign-rule (don't forget to lag the signal by one day)
  cum_perf_sign<-cumsum(na.omit(as.xts(lag(sign((yhat)),lag_fx)*x[,1])))
  sharpe_sign<-sqrt(250)*mean(diff(cum_perf_sign),na.rm=T)/sqrt(var(diff(cum_perf_sign),na.rm=T))
  
  # Plot
  if (plot_T)
  {
    par(mfrow=c(1,1))
    print(plot(cbind(x[,1],yhat),main="Returns (black) and filtered series (red)"))
    print(plot(cum_perf_sign,main=paste("Signum rule, ",colnames(x),", sharpe: ",round(sharpe_sign,3),sep="")))
    par(mfrow=c(1,1))
  }
  return(list(cum_perf_sign=cum_perf_sign,yhat=yhat,mdfa_obj_mse=mdfa_obj_mse,sharpe_sign=sharpe_sign))
}






