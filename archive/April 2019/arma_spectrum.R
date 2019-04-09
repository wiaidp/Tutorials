

arma_spectrum_func<-function(a1,b1,K,plot_T)
{
  arma_spec<-rep(1,K+1)
  for (j in 0:K)
  {
# Frequency    
    omega_j<-j*pi/K
    # MA-part
    ma<-1
    if (length(b1)>0)
    {
      for (i in 1:length(b1))
      {
        ma<-ma+b1[i]*exp(1.i*i*omega_j)
      }
    }
# AR-part
    ar<-1
    if (length(a1)>0)
    {
      for (i in 1:length(b1))
      {
        ar<-ar-a1[i]*exp(1.i*i*omega_j)
      }
    }
    arma_spec[j+1]<-ma/ar
  }
  if (plot_T)
  {
    plot(abs(arma_spec),type="l",main=paste("ARMA-spectrum, a1=",a1,", b1=",b1,sep=""),
         axes=F,xlab="Frequency",ylab="Amplitude",col="black")
    # We take 2-nd colname from weight_func because the first column is the target        
    axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
    axis(2)
    box()
    
  }
  return(arma_spec=arma_spec)
}