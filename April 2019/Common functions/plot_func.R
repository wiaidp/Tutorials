







mplot_func <- function(mplot, ax, plot_title, title_more, insamp, colo) {
  # if signal forecasts are computed then revision matrix is longer than data 
  # matrix
  ymin <- min(mplot, na.rm = TRUE)
  ymax <- max(mplot, na.rm = TRUE)
  # color vector for the lines and the text
  if(length(colo) == 0) {
    colo <- rainbow(ncol(mplot))
  }
  # plot the first data column
  plot(mplot[, 1], ylim = range(mplot, na.rm = TRUE), type = "l", axes = FALSE,
       col = colo[1], main = plot_title, xlab = "", ylab = "")
  abline(h = 0)
  # plot additional titles
  if(length(title_more) > 0) {
    mtext(title_more[1], col = colo[1])
  }
  # plot the remaining data columns
  if(ncol(mplot) > 1) {
    for(j in 2 : ncol(mplot)) {
      lines(mplot[, j], col = colo[j])
      # plot the remaining additional titles
      if(length(title_more) > 0) {
        mtext(title_more[j], col = colo[j], line = -(j - 1))
      }
      # plot a vertical line at the end of the in-sample window
      if((insamp[1] < 1.e+90) & (length(insamp) == ncol(mplot) - 1)) {
        abline(v = insamp[j - 1], col = colo[j])
        text(x = insamp[j - 1], y = -2 - (j - 1) * 0.3, col = colo[j],
             "End of In-Sample Window")
      }
    }
  }
  # add the x- and y-axis and draw a box around the current plot
  axis(1, at = 1 : nrow(mplot), labels = ax)
  axis(2)
  box()
}







plot_estimate_func<-function(mdfa_obj,weight_func,Gamma)
{
  par(mfrow=c(1,1))
  b<-mdfa_obj$b
  colo<-rainbow(ncol(b))
  plot(b[,1],type="l",main=paste("Filter coefficients",sep=""),
       axes=F,xlab="Lag",ylab="Coef",ylim=c(min(b),max(b)),col="black")
  mtext(colnames(weight_func)[2],line=-1,col="black")
  if (ncol(b)>1)
  {
    for (i in 2:ncol(b))
    {
      lines(b[,i],col=colo[i])
      # We take the i+1 colname from weight_func because the first column is the target        
      mtext(colnames(weight_func)[i+1],line=-i,col=colo[i])
      
    }
  }
  axis(1,at=1:L,labels=0:(L-1))
  axis(2)
  box()    
  
  
  plot(Arg(mdfa_obj$trffkt[,1])/((0:(nrow(weight_func)-1))*pi/(nrow(weight_func)-1)),type="l",main=paste("Time-shift concurrent, denseness=",K,sep=""),
       axes=F,xlab="Frequency",ylab="Time shift",col="black")
  # We take 2-nd colname from weight_func because the first column is the target        
  mtext(colnames(weight_func)[2],line=-1,col="black")
  if (ncol(abs(mdfa_obj$trffkt))>1)
  {
    for (i in 2:ncol(abs(mdfa_obj$trffkt)))
    {
      lines(abs(mdfa_obj$trffkt)[,i],col=colo[i])
      # We take the i+1 colname from weight_func because the first column is the target        
      mtext(colnames(weight_func)[i+1],line=-i,col=colo[i])
      
    }
  }
  axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
  axis(2)
  box()
  plot(abs(mdfa_obj$trffkt)[,1],type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
       axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,max(1,max(abs(mdfa_obj$trffkt)))))
  # We take 2-nd colname from weight_func because the first column is the target        
  mtext(colnames(weight_func)[2],line=-1,col="black")
  lines(Gamma,col="violet")
  mtext("Target Gamma",line=-2,col="violet")
  lines(0.5*abs(weight_func[,1])/max(abs(weight_func[,1])),col="red",lwd=3)
  mtext("Spectrum (red bold)",line=-3,col="red")
  if (ncol(abs(mdfa_obj$trffkt))>1)
  {
    for (i in 2:ncol(abs(mdfa_obj$trffkt)))
    {
      lines(abs(mdfa_obj$trffkt)[,i],col=colo[i])
      # We take the i+1 colname from weight_func because the first column is the target        
      mtext(colnames(weight_func)[i+1],line=-i-2,col=colo[i])
      
    }
  }
  axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
  axis(2)
  box()
  
}




perf_plot<-function(perf,sharpe,name)
{
  par(mfrow=c(1,1))
  colo<-rainbow(ncol(perf))
  print(plot(perf,main=paste("USD against G9, ",name),ylim=c(min(perf),max(perf))))
#  lines(perf[,1],col=colo[1])
  mtext(paste(colnames(perf)[1],", sharpe=",round(sharpe[1],3),sep=""),col=colo[1],line=-1)
  for (i in 2:ncol(perf))
  {  
#    lines(perf[,i],col=colo[i])
    mtext(paste(colnames(perf)[i],", sharpe=",round(sharpe[i],3),sep=""),col=colo[i],line=-i)
  }
  perf_agg<-as.xts(apply(perf,1,mean))
  sharpe_agg<-sqrt(250)*mean(diff(perf_agg),na.rm=T)/sqrt(var(diff(perf_agg),na.rm=T))
  lines(perf_agg,lwd=2,col="black")
  mtext(paste("Equally-weighted, sharpe=",round(sharpe_agg,3),sep=""),line=-i-1)
  return(list(perf_agg=perf_agg))
}


perf_plot_old<-function(perf,sharpe,name)
{
  par(mfrow=c(1,1))
  colo<-rainbow(ncol(perf))
  plot(perf[,1],main=paste("USD against G9, ",name),ylim=c(min(perf),max(perf)))
  lines(perf[,1],col=colo[1])
  mtext(paste(colnames(perf)[1],", sharpe=",round(sharpe[1],3),sep=""),col=colo[1],line=-1)
  for (i in 2:ncol(perf))
  {  
    lines(perf[,i],col=colo[i])
    mtext(paste(colnames(perf)[i],", sharpe=",round(sharpe[i],3),sep=""),col=colo[i],line=-i)
  }
  perf_agg<-as.xts(apply(perf,1,mean))
  sharpe_agg<-sqrt(250)*mean(diff(perf_agg),na.rm=T)/sqrt(var(diff(perf_agg),na.rm=T))
  lines(perf_agg,lwd=2,col="black")
  mtext(paste("Equally-weighted, sharpe=",round(sharpe_agg,3),sep=""),line=-i-1)
  return(list(perf_agg=perf_agg))
}





plot_compare_two_DFA_designs<-function(mdfa_obj,mdfa_obj_1,weight_func,weight_func_1,Gamma)
{
  
  par(mfrow=c(2,1))
  plot(abs(mdfa_obj$trffkt)[,1],type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
       axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,max(1,max(abs(mdfa_obj$trffkt)))))
  # We take 2-nd colname from weight_func because the first column is the target        
  mtext(colnames(weight_func)[2],line=-1,col="black")
  lines(Gamma,col="violet")
  mtext("Target Gamma",line=-2,col="violet")
  lines(0.5*abs(weight_func[,1])/max(abs(weight_func[,1])),col="red",lwd=3)
  mtext("Spectrum (red bold)",line=-3,col="red")
  if (ncol(abs(mdfa_obj$trffkt))>1)
  {
    for (i in 2:ncol(abs(mdfa_obj$trffkt)))
    {
      lines(abs(mdfa_obj$trffkt)[,i],col=colo[i])
      # We take the i+1 colname from weight_func because the first column is the target        
      mtext(colnames(weight_func)[i+1],line=-i-2,col=colo[i])
      
    }
  }
  axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
  axis(2)
  box()

  
  plot(abs(mdfa_obj_1$trffkt)[,1],type="l",main=paste("Amplitude concurrent, denseness=",K,sep=""),
       axes=F,xlab="Frequency",ylab="Amplitude",col="black",ylim=c(0,max(1,max(abs(mdfa_obj$trffkt)))))
  # We take 2-nd colname from weight_func because the first column is the target        
  mtext(colnames(weight_func_1)[2],line=-1,col="black")
  lines(Gamma,col="violet")
  mtext("Target Gamma",line=-2,col="violet")
  lines(0.5*abs(weight_func_1[,1])/max(abs(weight_func_1[,1])),col="red",lwd=3)
  mtext("Spectrum (red bold)",line=-3,col="red")
  if (ncol(abs(mdfa_obj_1$trffkt))>1)
  {
    for (i in 2:ncol(abs(mdfa_obj_1$trffkt)))
    {
      lines(abs(mdfa_obj_1$trffkt)[,i],col=colo[i])
      # We take the i+1 colname from weight_func because the first column is the target        
      mtext(colnames(weight_func_1)[i+1],line=-i-2,col=colo[i])
      
    }
  }
  axis(1,at=c(0,1:6*K/6+1),labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
  axis(2)
  box()
  
  par(mfrow=c(1,1))
  
  
    
}