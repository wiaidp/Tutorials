rm(list=ls())

###################################################
### code chunk number 18: exercise_dfa_ms_1
###################################################
# Generate series of length 2000
lenh<-2000
len<-120
# Specify the AR-coefficients
a_vec<-c(0.9,0.1,-0.9)
xh<-matrix(nrow=lenh,ncol=length(a_vec))
x<-matrix(nrow=len,ncol=length(a_vec))
yhat<-x
y<-x
# Generate series for each AR(1)-process
for (i in 1:length(a_vec))
{
  # We want the same random-seed for each process  
  set.seed(10)
  xh[,i]<-arima.sim(list(ar=a_vec[i]),n=lenh)
}


# Compute the coefficients of the symmetric target filter
cutoff<-pi/6
# Order of approximation
ord<-900
# Filter weights ideal trend (See DFA)
gamma<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))

ts.plot(gamma[1:36])

# Extract 120 observations in the midddle of the longer series
x<-xh[lenh/2+(-len/2):((len/2)-1),]
# Compute the outputs yt of the (truncated) symmetric target filter
for (i in 1:length(a_vec))
{
  for (j in 1:120)
  {
    y[j,i]<-gamma[1:900]%*%xh[lenh/2+(-len/2)-1+(j:(j-899)),i]+
      gamma[2:900]%*%xh[lenh/2+(-len/2)+(j:(j+898)),i]
  }
}



par(mfrow=c(3,1))
for (i in 1:3)   #i<-1
{
  ymin<-min(y[,i])
  ymax<-max(y[,i])
  ts.plot(x[,i],main=paste("Data and (blue) and target (red), a1 = ",a_vec[i],sep=""),col="blue",
          ylim=c(ymin,ymax),
          gpars=list(xlab="", ylab=""))
  lines(y[,i],col="red")
  mtext("data", side = 3, line = -1,at=len/2,col="blue")
  mtext("target", side = 3, line = -2,at=len/2,col="red")
}

