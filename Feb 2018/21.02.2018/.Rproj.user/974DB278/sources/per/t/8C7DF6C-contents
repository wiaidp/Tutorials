KF_gen<-function(parma,y,opti,outofsample,maxlik,xi10,P10,Fm,H,time_varying,x_data)
{
  len<-length(y)
  if (length(parma)>2)
  {
    Q<-diag(parma[1:(length(parma)-1)]^2)
  } else
  {
    Q<-as.matrix(parma[1]^2)
  }
  R<-parma[length(parma)]^2
  Pttm1<-array(dim=c(dim(Q),len+1))
  Pttm1[,,1:(len+1)]<-P10
  Ptt<-array(dim=c(dim(Q),len))
  xttm1<-xi10
  logl<-0.
  xittm1<-matrix(nrow=len,ncol=(length(parma)-1))
  xittm1[1,]<-xi10
  xitt<-xittm1
  # If time_varying==T then we fill data into H (SSM with time-varying coefficients)
  if (time_varying)
  {
    if (is.null(x_data))
    {
      # For an autoregressive model we fill in past y's
      H<-c(y[1:dim(Q)[2]])
      # We need the first y[1:p] in H: therefore the first equation will be for t=p+1
      anf<-dim(Q)[2]+1
    } else
    {
      # For a regression model we fill in the explanatory data
      H<-x_data[1,]
      anf<-1
    }
  } else
  {
    anf<-1
  }
  # Kalman-Recursion: starts in i=dim(Q)[2] for a time series model
  # and in i=1 for a regression model
  
  for (i in anf:len)        #i<-1     H<-c(1,0)       xitt[,2]
  {
    # Kalman-Gain
    He<-(H%*%(Pttm1[,,i]%*%H))[1,1]+R
    epshatoutofsample<-y[i]-(H%*%xttm1)[1,1]
    xittm1[i,]<-xttm1
    xtt<-xttm1+Pttm1[,,i]%*%H*epshatoutofsample/He
    epshatinsample<-y[i]-(H%*%xtt)[1,1]
    xitt[i,]<-xtt
    xttm1<-Fm%*%xtt
    Ptt[,,i]<-Pttm1[,,i]-((Pttm1[,,i]%*%H)%*%(H%*%Pttm1[,,i]))/He
    Pttm1[,,i+1]<-Fm%*%Ptt[,,i]%*%t(Fm)+Q
    if (time_varying)
    {
      if (is.null(x_data))
      {
        # For an autoregressive model we fill past y's in H
        H<-c(y[i-dim(Q)[2]+1:(dim(Q)[2])])
      } else
      {
        # For a regression model we fill the explanatory data in H
        H<-x_data[min(len,i+1),]
      }
    }
    # Here we specify the optimization criterion: least-squares
    # or maximum likelihood, in-sample or out-of-sample
    if (outofsample)
    {
      if (maxlik)
      {
        logl<-logl+log(He)+epshatoutofsample^2/He
      } else
      {
        logl<-logl+epshatoutofsample^2
      }
    } else
    {
      if (maxlik)
      {
        logl<-logl+log(He)+epshatinsample^2/He
      } else
      {
        logl<-logl+epshatinsample^2
      }
    }
  }
  if (opti)
  {
    return(logl/len)
  } else
  {
    return(list(logl=logl/len,xitt=xitt,xittm1=xittm1,Ptt=Ptt,Pttm1=Pttm1))
  }
}






use_diff<-T
lag_k<-1
FX_pair<-"USDEUR"

if (use_diff)
{
  xhh<-as.matrix(lag(FX_scale,lag_k))
  x<-xh<-xhh[-(1:lag_k),]
  
  head(x)
  yh<-as.double(FX_scale[,FX_pair])
  y<-yh[-(1:lag_k)]
  
  x_diff<-apply(x,2,diff)
  y_diff<-diff(y)
  
  lm_obj<-lm(y_diff~x_diff)
  
  summary(lm_obj)
  R<-sqrt(var(lm_obj$res))
  x_mat<-cbind(rep(1,ncol(x_diff)),x_diff)
  y_dep<-y_diff          
} else
{
  xhh<-as.matrix(lag(FX_scale,lag_k))
  x<-xh<-xhh[-(1:lag_k),]
  
  head(x)
  yh<-as.double(FX_scale[,FX_pair])
  y<-yh[-(1:lag_k)]
  
  lm_obj<-lm(y~x)
  
  summary(lm_obj)
  R<-sqrt(var(lm_obj$res))
  x_mat<-cbind(rep(1,ncol(x),x))
  y_dep<-y          
          
}





Fm<-diag(dim(x_mat)[2])
H<-NULL

P10<-diag(rep(10000000*var(y_dep),dim(x_mat)[2]))   #ts.plot(y_dep)
xi10<-rep(0,dim(x_mat)[2])
xi10<-lm_obj$coef
maxlik<-T
outofsample<-T
parma<-c(rep(0,dim(x_mat)[2]),R)
time_varying<-T
opti<-T

if (opti)
{
  objopt<-nlminb(start=parma,objective=KF_gen,
               y=y_dep,opti=opti,outofsample=outofsample,maxlik=maxlik,xi10=xi10,P10=P10,Fm=Fm,H=H,time_varying=time_varying,x_data=x_mat)

  parma<-objopt$par
}

parma^2


opti<-F
obj<-KF_gen(parma,y_dep,opti,outofsample,maxlik,xi10,P10,Fm,H,time_varying,x_mat)

# Criterion value: compare with criterion value of fixed coefficient model
obj$logl

#----------------
# Exercise 4
# Residuals

res_adaptive<-y_dep-apply(obj$xittm1*x_mat,1,sum)

ts.plot(res_adaptive)

#-------------
# Exercise 5
var(res_adaptive)#0.0009516301


# Plot of coefficients 
len<-5000
enf<-nrow(obj$xitt)
anf<-max(1,enf-len)

colo<-rainbow(dim(x_mat)[2])
ts.plot(obj$xitt[anf:enf,1],ylim=c(min(obj$xitt[anf:enf,]),max(obj$xitt[anf:enf,])),col=colo[1],
        main=paste("maxlik=",maxlik,", outofsample=",outofsample,", R=",
                   round(parma[length(parma)]^2,3),", Q=",
                   round(parma[1]^2,3),
                   round(parma[2]^2,3),
                   round(parma[3]^2,3),
                   round(parma[4]^2,3),
                   round(parma[5]^2,3),sep=""),
        xlab="",ylab="")
mtext(paste("Variable ",1,sep=""), side = 3, line = -1,at=len/2,col=colo[1])
for (i in 2:dim(x_mat)[2])
{
  lines(obj$xitt[anf:enf,i],col=colo[i])
  mtext(colnames(x_mat)[i], side = 3, line = -i,at=len/2,col=colo[i])
}
abline(h=0,lwd=2)



# Statistical significance: add 90% intervals

sig<-1.6
colo<-rainbow(dim(x_mat)[2])

ts.plot(obj$xitt[anf:enf,1],ylim=c(min(obj$xitt[anf:enf,]),max(obj$xitt[anf:enf,])),col=colo[1],
        main=paste("maxlik=",maxlik,", outofsample=",outofsample,", R=",
                   round(parma[length(parma)]^2,3),", Q=",
                   round(parma[1]^2,3),", ",
                   round(parma[2]^2,3),", ",
                   round(parma[3]^2,3),", ",
                   round(parma[4]^2,3),", ",
                   round(parma[5]^2,3),sep=""),
        xlab="",ylab="")
lines(obj$xitt[anf:enf,1]+sig*sqrt(obj$Ptt[1,1,anf:enf]),col=colo[1],lty=2)
lines(obj$xitt[anf:enf,1]-sig*sqrt(obj$Ptt[1,1,anf:enf]),col=colo[1],lty=2)
mtext(paste("Variable ",1,sep=""), side = 3, line = -1,at=len/2,col=colo[1])
for (i in 2:dim(x_mat)[2])
{
  lines(obj$xitt[anf:enf,i],col=colo[i])
  lines(obj$xitt[anf:enf,i]+sig*sqrt(obj$Ptt[i,i,anf:enf]),col=colo[i],lty=2)
  lines(obj$xitt[anf:enf,i]-sig*sqrt(obj$Ptt[i,i,anf:enf]),col=colo[i],lty=2)
  mtext(colnames(x_mat)[i], side = 3, line = -i,at=len/2,col=colo[i])
}
abline(h=0)


#Plot absolute sum of coefficients: spread of data
#   Larger values are indicative for colinearity
ts.plot(apply(abs(obj$xitt[anf:enf,]),1,sum))



#--------------------------------------------------------------------------
# Trading

par(mfrow=c(2,1))
ts.plot(cumsum(sign(apply(obj$xittm*x_mat,1,sum))*y_dep))
acf(sign(apply(obj$xittm*x_mat,1,sum))*y_dep)
ts.plot(cumsum(y_dep))
acf(y_dep)

