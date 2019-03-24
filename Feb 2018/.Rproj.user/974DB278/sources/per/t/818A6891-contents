

###################################################
### code chunk number 20: exercise_dfa_ms_3
###################################################
plot_T<-F
periodogram<-matrix(ncol=3,nrow=len/2+1)
trffkt<-periodogram
perf_mat<-matrix(nrow=3,ncol=2)
dimnames(perf_mat)[[2]]<-c("Criterion Value",
                           "Mean-Square Sample Filter Error")
dimnames(perf_mat)[[1]]<-c("a1=0.9","a1=0.1","a1=-0.9")
# Filter length
L<-12
# Real-time design
Lag<-0
# Target ideal trend
Gamma<-c(1,(1:(len/2))<len/12)

plot(Gamma,type="l",main="Target",
     axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                     "4pi/6","5pi/6","pi"))
box()

b<-matrix(nrow=L,ncol=3)
# Compute real-time filters
for (i in 1:3)#i<-1
{
  # Compute the periodogram based on the data (length 120)  
  periodogram[,i]<-per(x[,i],plot_T)$per
  # Optimize filters
  filt<-dfa_ms(L,periodogram[,i],Lag,Gamma)
  trffkt[,i]<-filt$trffkt
  b[,i]<-filt$b
  # Compute real-time outputs (we can use the longer series in order 
  # to obtain estimates for time points t=1,...,11)
  for (j in 1:len)
    yhat[j,i]<-filt$b%*%xh[lenh/2+(-len/2)-1+j:(j-L+1),i]
}


###################################################
### code chunk number 21: exercise_dfa_ms_3
###################################################
for (i in 1:3)
{
  # Compute criterion values
  perf_mat[i,1]<-(2*pi/length(Gamma))*
    abs(Gamma-trffkt[,i])^2%*%periodogram[,i]
}
perf_mat[,1]


###################################################
### code chunk number 22: exercise_dfa_ms_4
###################################################
# Compute time-domain MSE
mse<-apply(na.exclude((yhat-y))^2,2,mean)
perf_mat[,2]<-mse
round(perf_mat[,2],3)


###################################################
### code chunk number 23: z_dfa_ar1_output.pdf
###################################################
library(Hmisc)
require(xtable)
#latex(cor_vec, dec = 1, , caption = "Example of using latex to create table",
#center = "centering", file = "", floating = FALSE)
xtable(perf_mat, dec = 1,digits=rep(3,dim(perf_mat)[2]+1),
       paste("Criterion values vs. sample (mean-square) filter errors",sep=""),
       label=paste("perf_mat",sep=""),
       center = "centering", file = "", floating = FALSE)


###################################################
### code chunk number 24: z_dfa_ar1_output.pdf
###################################################
par(mfrow=c(3,1))
for (i in 1:3)   #i<-1
{
  ymin<-min(min(y[,i]),min(na.exclude(yhat)[,i]))
  ymax<-max(max(y[,i]),max(na.exclude(yhat)[,i]))
  ts.plot(yhat[,i],main=paste("Time-domain MSE = ",
                              round(mse[i],3)," , Frequency-domain MSE = ",
                              round(perf_mat[i,1],3),", a1 = ",a_vec[i],sep=""),col="blue",
          ylim=c(ymin,ymax),
          gpars=list(xlab="", ylab=""))
  lines(y[,i],col="red")
  mtext("Real-time", side = 3, line = -1,at=len/2,col="blue")
  mtext("target", side = 3, line = -2,at=len/2,col="red")
}


###################################################
### code chunk number 25: z_dfa_ar1_output.pdf
###################################################
file = paste("z_dfa_ar1_sym_output", sep = "")
cat("\\begin{figure}[H]")
cat("\\begin{center}")
cat("\\includegraphics[height=6in, width=6in]{", file, "}\n",sep = "")
cat("\\caption{Real-time filter output (blue) vs. targets (red) for a1=0.9 (top), a1=0.1 (middle) and a1=-0.9 (bottom)", sep = "")
cat("\\label{z_dfa_ar1_sym_output}}", sep = "")
cat("\\end{center}")
cat("\\end{figure}")


###################################################
### code chunk number 26: z_dfa_ar1_output.pdf
###################################################

omega_k<-pi*0:(len/2)/(len/2)
par(mfrow=c(2,2))
amp<-abs(trffkt)
shift<-Arg(trffkt)/omega_k
plot(amp[,1],type="l",main="Amplitude functions",
     axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,1))
lines(amp[,2],col="orange")
lines(amp[,3],col="green")
lines(Gamma,col="violet")
mtext("Amplitude a1=0.9", side = 3, line = -1,at=len/4,col="black")
mtext("Amplitude a1=0.1", side = 3, line = -2,at=len/4,col="orange")
mtext("Amplitude a1=-0.9", side = 3, line = -3,at=len/4,col="green")
mtext("Target", side = 3, line = -4,at=len/4,col="violet")
axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                     "4pi/6","5pi/6","pi"))
axis(2)
box()
plot(shift[,1],type="l",main="Time-shifts",
     axes=F,xlab="Frequency",ylab="Shift",col="black",
     ylim=c(0,max(na.exclude(shift[,3]))))
lines(shift[,2],col="orange")
lines(shift[,3],col="green")
lines(rep(0,len/2+1),col="violet")
mtext("Shift a1=0.9", side = 3, line = -1,at=len/4,col="black")
mtext("Shift a1=0.1", side = 3, line = -2,at=len/4,col="orange")
mtext("Shift a1=-0.9", side = 3, line = -3,at=len/4,col="green")
mtext("Target", side = 3, line = -4,at=len/4,col="violet")
axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                     "4pi/6","5pi/6","pi"))
axis(2)
box()
plot(periodogram[,1],type="l",main="Periodograms",
     axes=F,xlab="Frequency",ylab="Periodogram",col="black",
     ylim=c(0,max(periodogram[,3])/6))
lines(periodogram[,2],col="orange")
lines(periodogram[,3],col="green")
mtext("Periodogram a1=0.9", side = 3, line = -1,at=len/4,col="black")
mtext("Periodogram a1=0.1", side = 3, line = -2,at=len/4,col="orange")
mtext("Periodogram a1=-0.9", side = 3, line = -3,at=len/4,col="green")
axis(1,at=c(0,1:6*len/12+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                     "4pi/6","5pi/6","pi"))
axis(2)
box()




