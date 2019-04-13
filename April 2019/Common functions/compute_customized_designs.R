

compute_customized_designs_func<-function(lambda_vec,eta_vec,L,weight_func,Lag,Gamma,cutoff)
{
  for (i in 1:length(lambda_vec))#i<-1
  {
    lambda<-lambda_vec[i]
    eta<-eta_vec[i]
    
    mdfa_obj<-MDFA_cust(L,weight_func,Lag,Gamma,cutoff,lambda,eta)$mdfa_obj
    # Keep track of transferfunctions, filter coefficients and of filtered outputs 
    if (i==1)
    {
      trffkt_mat<-mdfa_obj$trffkt
      b_mat<-mdfa_obj$b
      yhat_mat<-filt_func(x,mdfa_obj$b)$yhat
      
    } else
    {
      trffkt_mat<-cbind(trffkt_mat,mdfa_obj$trffkt)
      b_mat<-cbind(b_mat,mdfa_obj$b)
      yhat_mat<-cbind(yhat_mat,filt_func(x,mdfa_obj$b)$yhat)
    } 
  }
  
  # Compare amplitude and time shifts
  par(mfrow=c(2,2))
  mplot<-abs(trffkt_mat)
  dimnames(mplot)[[2]]<-paste("Amplitude (",lambda_vec,",",eta_vec,")",sep="")
  ax<-rep(NA,ncol(mplot))
  ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
  plot_title<-"Amplitude functions"
  insamp<-1.e+90
  title_more<-dimnames(mplot)[[2]]
  colo<-rainbow(ncol(mplot))
  mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
  mplot<-abs(trffkt_mat)
  # Scale amplitudes for easier visual inspection
  for (i in 1:ncol(mplot))
    mplot[,i]<-mplot[,i]/max(mplot[,i])
  dimnames(mplot)[[2]]<-paste("Scaled amplitude (",lambda_vec,",",eta_vec,")",sep="")
  ax<-rep(NA,ncol(mplot))
  ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
  plot_title<-"Scaled amplitude functions"
  insamp<-1.e+90
  title_more<-dimnames(mplot)[[2]]
  colo<-rainbow(ncol(mplot))
  mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
  mplot<-Arg(trffkt_mat)/(pi*(0:K)/K)
  dimnames(mplot)[[2]]<-paste("Time-shifts (",lambda_vec,",",eta_vec,")",sep="")
  ax<-rep(NA,ncol(mplot))
  ax[1+(0:6)*((nrow(mplot)-1)/6)]<-c(0,"pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
  plot_title<-"Time-shifts"
  insamp<-1.e+90
  title_more<-dimnames(mplot)[[2]]
  colo<-rainbow(ncol(mplot))
  mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
  
  
  # Compare filter outputs
  anf<-250
  enf<-300
  par(mfrow=c(1,2))
  mplot<-yhat_mat[anf:enf,]
  dimnames(mplot)[[2]]<-paste("Filter outputs (",lambda_vec,",",eta_vec,")",sep="")
  ax<-rep(NA,ncol(mplot))
  ax[1+(0:6)*((nrow(mplot)-1)/6)]<-as.integer(1+(0:6)*((nrow(mplot)-1)/6))
  plot_title<-"Filter outputs"
  insamp<-1.e+90
  title_more<-dimnames(mplot)[[2]]
  colo<-rainbow(ncol(mplot))
  mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
  
  
  # Compare scaled filter outputs
  mplot<-yhat_mat[anf:enf,]
  # Scale outputs for easier visual inspection
  for (i in 1:ncol(mplot))
    mplot[,i]<-mplot[,i]/max(abs(trffkt_mat)[,i])
  dimnames(mplot)[[2]]<-paste("Filter outputs (",lambda_vec,",",eta_vec,")",sep="")
  ax<-rep(NA,ncol(mplot))
  ax[1+(0:6)*((nrow(mplot)-1)/6)]<-as.integer(1+(0:6)*((nrow(mplot)-1)/6))
  plot_title<-"Scaled filter outputs"
  insamp<-1.e+90
  title_more<-dimnames(mplot)[[2]]
  colo<-rainbow(ncol(mplot))
  mplot_func(mplot, ax, plot_title, title_more, insamp, colo)
  
  
}