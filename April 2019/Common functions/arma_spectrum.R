# a1<-ar_burg
#b1<-ma_burg

arma_spectrum_func<-function(a1,b1,K,plot_T)
{
  arma_spec<-rep(1,K+1)
  for (j in 0:K)#j<-0
  {
# Frequency    
    omega_j<-j*pi/K
    # MA-part
    ma<-1
    if (!is.null(b1)>0)
    {
      for (i in 1:length(b1))
      {
        ma<-ma+b1[i]*exp(1.i*i*omega_j)
      }
    }
# AR-part
    ar<-1
    if (!is.null(a1)>0)
    {
      for (i in 1:length(a1))
      {
        ar<-ar-a1[i]*exp(1.i*i*omega_j)
      }
    }
    arma_spec[j+1]<-ma/ar
  }
  if (plot_T)
  {
    plot(abs(arma_spec),type="l",main=paste("ARMA-spectrum, a1=",ifelse(is.null(a1),"NULL",round(a1,3)),", b1=",ifelse(is.null(b1),"NULL",round(b1,3)),sep=""),
         axes=F,xlab="Frequency",ylab="Amplitude",col="black")
    # We take 2-nd colname from weight_func because the first column is the target        
    axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
    axis(2)
    box()
    
  }
  return(list(arma_spec=arma_spec))
}



#a1<-rbind(c(0.5,0.3),c(0.6,0.2))
#a1<-rbind(c(0.5,0.),c(0.,0.2))
#b1<-rbind(c(0.5,0.3),c(0.6,0.2))
#b1<-rbind(c(0.,0.),c(0.,0.))
#sigma<-rbind(c(1,0.5),c(0.5,1))

varma_spectrum_func<-function(a1,b1,K,plot_T)
{
  varma_spec<-array(dim=c(K+1,dim(a1)))
  varma_spec[1,,]
  for (i in 1:(K+1))
    varma_spec[i,,]<-diag(1,dim(a1)[1])
  for (j in 0:K)#j<-1
  {
    # Frequency    
    omega_j<-j*pi/K
    # MA-part
    ma<-diag(1,dim(a1)[1])
    if (!is.null(b1)>0)
    {
      ma<-ma+b1*exp(1.i*omega_j)
    }
    # AR-part
    ar<-diag(1,dim(a1)[1])
    if (!is.null(a1)>0)
    {
      ar<-ar-a1*exp(1.i*omega_j)
    }
    varma_spec[j+1,,]<-solve(ar)%*%ma%*%sigma
  }
  if (plot_T)
  {
    for (i in 1:dim(a1)[1])#i<-10
    {
      plot(abs(varma_spec[,i,i]),type="l",main=paste("ARMA-spectrum, a1=",ifelse(is.null(a1),"NULL",round(a1,3)),", b1=",ifelse(is.null(b1),"NULL",round(b1,3)),sep=""),
           axes=F,xlab="Frequency",ylab="Amplitude",col="black")
      # We take 2-nd colname from weight_func because the first column is the target        
      axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                        "4pi/6","5pi/6","pi"))
      axis(2)
      box()
    }
  }
  return(list(varma_spec=varma_spec))
}

#varma_spec[10,,]



trffkt_func<-function(b,K,plot_T)
{
  trffkt<-rep(0,K+1)
  for (j in 0:K)#j<-0
  {
# Frequency    
    omega_j<-j*pi/K
    trffkt[j]<-b%*%exp(1.i*(0:(length(b)-1))*omega_j)
  }
  if (plot_T)
  {
    plot(abs(trffkt),type="l",main=paste("Amplitude, b=",ifelse(is.null(b),"NULL",round(b,3)),sep=""),
         axes=F,xlab="Frequency",ylab="Amplitude",col="black")
    # We take 2-nd colname from weight_func because the first column is the target        
    axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                      "4pi/6","5pi/6","pi"))
    axis(2)
    box()
    
  }
  return(list(trffkt=trffkt))
}
