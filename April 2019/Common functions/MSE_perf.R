
#mdfa_obj<-mdfa_obj_mixed
#yhat<-yhat_mixed
MSE_perf_func<-function(insamp,y,yhat,len,mdfa_obj,data_mat)
{


  perf_mat<-mean((y-yhat)[1:insamp]^2,na.rm=T)
  if (insamp<len)
  {
    perf_mat<-c(perf_mat,mean((y-yhat)[(insamp+1):len]^2,na.rm=T))
  } else
  {
    perf_mat<-c(perf_mat,NA)
  }
  if (!target_as_explanatory)
  {
    print("The next MSE is cheating in this configuration because it uses low-freq data as estimate")
  }  
  mean((y-data_mat[,1])^2,na.rm=T)
  perf_mat<-c(perf_mat,mdfa_obj$MS_error)
  #mean(diff(diff(yhat_mixed[L:length(yhat_mixed)]))^2)
  perf_mat<-t(as.matrix(as.vector(perf_mat),nrow=1,ncol=length(perf_mat)))
  colnames(perf_mat)<-c("In-sample MSE", "Out-sample MSE", "DFA criterion")
  return(list(perf_mat=perf_mat))
}