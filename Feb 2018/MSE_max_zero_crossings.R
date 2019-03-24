
# Input series
omega<-pi/61
len<-120
x<-cos(1:len*omega)

# Filter weights: equally-weighted filter of length 12

L<-12
bkh<-rep(1/L,L)
# Normalization (so that scale of output is equal to scale of input)
norm_c<-abs(sum(bkh*exp(1.i*(0:(L-1))*omega)))
bk<-bkh/norm_c

# Filter series
yhat<-rep(NA,len)
for (i in L:len)
  yhat[i]<-sum(bk*x[i:(i-L+1)])

# Plot series, filtered series and squared error
par(mfrow=c(2,1))
ts.plot(x,col="blue",main="Series(blue), filtered series (red) and squared error (black)")
lines(yhat,col="red")
abline(h=0)
ts.plot((x-yhat)^2)