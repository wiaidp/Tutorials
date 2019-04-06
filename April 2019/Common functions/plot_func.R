# Plots performances as a function of weekday
# series.xts<-perf_ss[start_date]
# name<-NULL
plot_weekday_func<-function(series.xts,name)
{
  mplot_diff<-diff(as.xts(series.xts))
  weekday<-c("Monday","Tuesday","Wednesday","Thursday","Friday")
  j<-0
  for (i in 1:5)#i<-4
  {
    new_mplot_diff<-mplot_diff[.indexwday(mplot_diff) %in% i]
    colo<-"blue"
    if (length(new_mplot_diff)>2)
    {
      j<-j+1
      mplot<-cumsum(na.exclude(new_mplot_diff))
      ax<-rownames(mplot)
      sharpe_vec<-sqrt(252/4)*apply(apply(mplot,2,diff),2,mean)/sqrt(apply(apply(mplot,2,diff),2,var))
      plot_title<-paste(name,": day ",i,", sharpe: ",round(sharpe_vec,2),sep="")
      plot(x = mplot[,1], xlab = "Time", 
           main = plot_title,  major.ticks= "quarters",major.format="%Y-%b-%d ",
           minor.ticks = F, col = colo[1])
      if (j>1)
      {
        boxplot_mat<-rbind(boxplot_mat,cbind(new_mplot_diff,i))
      } else
      {
        boxplot_mat<-cbind(new_mplot_diff,i)
      }
      acf(na.exclude(new_mplot_diff),demean=F)
      
    } else
    {
      print(paste("No data for day ",i,sep=""))
    }
  }
  colnames(boxplot_mat)<-c("perf","day")
  boxplot(perf ~ day, data = boxplot_mat, col = "lightgray")
  
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
  axis(1,at=1:L,labels=1:L)
  axis(2)
  box()    
  
  
  plot(Arg(mdfa_obj$trffkt[,1])/((0:(nrow(weight_func)-1))*pi/(nrow(weight_func)-1)),type="l",main=paste("Time-shift concurrent, denseness=",K,sep=""),
       axes=F,xlab="Frequency",ylab="Amplitude",col="black")
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
  lines(0.5*abs(weight_func[,1])/max(abs(weight_func[,1])),col="black",lwd=3)
  mtext("Spectrum (black bold)",line=-3,col="black")
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
