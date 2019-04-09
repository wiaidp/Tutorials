ideal_filter_func<-function(periodicity,M,x)
{
  len<-length(x)  
  cutoff<-pi/periodicity
  # Compute coefficients gamma
  gamma<-c(cutoff/pi,(1/pi)*sin(cutoff*1:M)/(1:M))
  # Check 
  sum(gamma)+sum(gamma[2:M])
  y<-rep(NA,len)
  
  for (k in M:(len-M))#k<-0
  {
    y[k]<-gamma[1:M]%*%x[k+(0:(M-1))]+gamma[2:M]%*%x[k-(1:(M-1))]
  }
  mean_holding_time<-mean(diff(which(na.exclude(sign(y[1:(len-1)])!=sign(y[2:len])))))
  gamma<-c(gamma[M:1],gamma[2:M])
  return(list(y=y,mean_holding_time=mean_holding_time,gamma=gamma))
}
