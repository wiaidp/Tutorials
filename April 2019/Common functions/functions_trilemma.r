# Copyright: Marc Wildi
# 30.10.2014
# http://blog.zhaw.ch/sef/
# http://www.mdfapartners.com/






# This function computes ATS-components and all relevant performance measures
#L_sym<-940*2

Performance_func<-function(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,
cutoff,L,L_sym,a1,mba,scaled_ATS,estim_MBA,Lag,i1,i2)
{
  K<-length(weight_func)-1
  omega_k<-(0:K)*pi/K
  #
  amp<-matrix(ncol=length(lambda_vec)+1,nrow=K+1)
  shift<-amp
  selectivity_insamp<-rep(NA,dim(amp)[2])
  mean_shift_insamp<-rep(NA,dim(amp)[2])
  curvature_insamp<-rep(NA,dim(amp)[2])
  selectivity_outsamp<-rep(NA,dim(amp)[2])
  mean_shift_outsamp<-rep(NA,dim(amp)[2])
  curvature_outsamp<-rep(NA,dim(amp)[2])

  b<-matrix(ncol=dim(amp)[2],nrow=L)
  # Determine passband
  omega_Gamma<-round(cutoff*K/pi)+1
  passband<-1:omega_Gamma
  stopband<-(omega_Gamma+1):(K+1)
  Accuracy<-rep(NA,dim(amp)[2])
  Smoothness<-Accuracy
  Timeliness<-Accuracy
  Residual<-Accuracy
  anf_cor<--10
  enf_cor<-10
  lag_cor<-anf_cor:enf_cor
  peak_cor_outsamp<-peak_cor_insamp<-matrix(1:(length(lag_cor)*(dim(amp)[2])),nrow=length(lag_cor),ncol=dim(amp)[2])
  lag_max_insamp<-rep(NA,dim(amp)[2])
  lag_max_outsamp<-rep(NA,dim(amp)[2])
  xf<-matrix(nrow=len1,ncol=dim(amp)[2])
  sample_mse<-matrix(ncol=4,nrow=dim(amp)[2])
  lambda<-0
  eta<-0
  ord<-as.integer((L_sym-1)/2)
# Filter weights ideal trend (See DFA)
  gamma<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))

# Compute outputs on long sample
  xf_sym<-rep(NA,len1)
#  for (j in (1+ord):(len1-ord))#j<-2+ord   length(j:(j-ord/2))
  for (j in (1+ord):min((len1-ord+len),len1))#j<-2+ord   length(j:(j-ord/2))
    xf_sym[j]<-gamma[1:((ord/2)+1)]%*%x1[j:(j-ord/2)]+gamma[2:((ord/2)+1)]%*%x1[(j+1):(j+ord/2)]#ts.plot(xf_sym)
  weight_func_best<-rep(1,K+1)/(abs(1-a1*exp(1.i*omega_k))^2*2*pi)               #a1<--0
#
# Estimation loop for all eta/lambda values
  for (i in 0:(dim(amp)[2]-1))         #i<-0
  {
    if (i==0)
    {
      if (mba&!estim_MBA)
      {
        dfa<-dfa_analytic(L,0,weight_func,Lag,Gamma,0,cutoff,i1,i2)
        weight_func_p<-weight_func
      } else
      {
        dfa<-dfa_analytic(L,0,weight_func_best,Lag,Gamma,0,cutoff,i1,i2)
        weight_func_p<-weight_func_best
      }
    } else
    {
      dfa<-dfa_analytic(L,lambda_vec[i],weight_func,Lag,Gamma,eta_vec[i],cutoff,i1,i2)#ts.plot(weight_func)
      weight_func_p<-weight_func
    }


    b[,i+1]<-dfa$b                   #
    amp[,i+1]<-abs(dfa$trffkt)
# Selectivity: immunized against scaling effects
    selectivity_insamp[i+1]<-(amp[1:(K*cutoff/pi),i+1]^2%*%weight_func_p[1:(K*cutoff/pi)])/
    (amp[(K*cutoff/pi+1):(K+1),i+1]^2%*%(weight_func_p[(K*cutoff/pi+1):(K+1)]*
    (1+(0:(K-K*cutoff/pi)))^2))
    shift[1+1:K,i+1]<-Arg(dfa$trffkt)[1+1:K]/((pi*1:K)/K)
# Shift in frequency zero
    shift[1,i+1]<-b[,i+1]%*%(0:(L-1))/amp[1,i+1]
# Mean-shift statistic: immunized against scaling effects
    mean_shift_insamp[i+1]<-mean(shift[1:(K*cutoff/pi),i+1])
    if (scaled_ATS&i>0)
    {
# Rescaling customized filters such that amplitude is not shrunken: it is brought up to mean levels of MSE filter
      amp_error<-((Gamma-amp[,i+1]*mean(amp[1:(omega_Gamma),1])/mean(amp[1:(omega_Gamma),i+1])))^2*weight_func_p/K
      shift_error<-4*Gamma*amp[,i+1]*(mean(amp[1:(omega_Gamma),1])/mean(amp[1:(omega_Gamma),i+1]))*sin(shift[,i+1]*omega_k/2)^2*weight_func_p/K
# Rescaling customized filters such that amplitude is not shrunken: it is brought up to mean level of one
      amp_error<-((Gamma-amp[,i+1]*1/mean(amp[1:(omega_Gamma),i+1])))^2*weight_func_p/K
      shift_error<-4*Gamma*amp[,i+1]*(1/mean(amp[1:(omega_Gamma),i+1]))*sin(shift[,i+1]*omega_k/2)^2*weight_func_p/K
    } else
    {
      amp_error<-((Gamma-amp[,i+1]))^2*weight_func_p/K
      shift_error<-4*Gamma*amp[,i+1]*sin(shift[,i+1]*omega_k/2)^2*weight_func_p/K
    }
    Accuracy[i+1]<-sum(amp_error[passband])
    Smoothness[i+1]<-sum(amp_error[stopband])
    Timeliness[i+1]<-sum(shift_error[passband])
    Residual[i+1]<-sum(shift_error[stopband])
    yhat<-rep(NA,len1)
    for (j in L:len1)
      yhat[j]<-b[,i+1]%*%x1[j:(j-L+1)]
    xf[,i+1]<-yhat
# Relative Curvature (immunized against scaling effects)
    curvature_insamp[i+1]<-mean(diff(yhat[max(L_sym/2,L)+1:len],diff=2)^2,na.rm=T)/
    var(yhat[max(L_sym/2,L)+1:len],na.rm=T)
    curvature_outsamp[i+1]<-mean(diff(yhat[len+max(L_sym/2,L)+1:len],diff=2)^2,na.rm=T)/
    var(yhat[len+max(L_sym/2,L)+1:len],na.rm=T)
# Peak Correlation
    i_cor<-0
    for (j in lag_cor)
    {
      i_cor<-i_cor+1
      if (j<=0)
      {
        peak_cor_insamp[i_cor,i+1]<-cor(xf_sym[max(L_sym/2,L)+(1-j):len],yhat[max(L_sym/2,L)+1:(len+j)])
        peak_cor_outsamp[i_cor,i+1]<-cor(xf_sym[len+max(L_sym/2,L)+(1-j):len],yhat[len+max(L_sym/2,L)+1:(len+j)])
      } else
      {
        peak_cor_insamp[i_cor,i+1]<-cor(xf_sym[max(L_sym/2,L)+1:(len-j)],yhat[max(L_sym/2,L)+(1+j):len])
        peak_cor_outsamp[i_cor,i+1]<-cor(xf_sym[len+max(L_sym/2,L)+1:(len-j)],yhat[len+max(L_sym/2,L)+(1+j):len])
      }
    }
    lag_max_insamp[i+1]<-ifelse(sum(is.na(peak_cor_insamp[,i+1]))==0,which(peak_cor_insamp[,i+1]==max(peak_cor_insamp[,i+1]))-1,NA)+min(0,anf_cor)
    lag_max_outsamp[i+1]<-which(peak_cor_outsamp[,i+1]==max(peak_cor_outsamp[,i+1]))-1+min(0,anf_cor)
# In-sample MSE
    sample_mse[i+1,1]<-mean(((xf_sym-yhat)[max(L_sym/2,L)+1:len])^2,na.rm=T) #xf_sym[max(L_sym/2,L)+1:len]
# Out-of-sample MSE
    sample_mse[i+1,2]<-mean(((xf_sym-yhat)[(len+max(L_sym/2,L))+1:len])^2,na.rm=T)#yhat[(len+max(L_sym/2,L))+1:len]
# determine in-sample scaling: only for DFA-filters (not for best MSE)
    if (i>0)
    {
      alpha_in<-cov(xf_sym[max(L_sym/2,L)+1:len],yhat[max(L_sym/2,L)+1:len],use="pairwise.complete.obs")/var(yhat[max(L_sym/2,L)+1:len],na.rm=T)
# In-sample scaled MSE
      sample_mse[i+1,3]<-mean(((xf_sym-alpha_in*yhat)[max(L_sym/2,L)+1:len])^2,na.rm=T)
# Out-of-sample scaled MSE
      sample_mse[i+1,4]<-mean(((xf_sym-alpha_in*yhat)[len+max(L_sym/2,L)+1:len])^2,na.rm=T)
    } else
    {
      sample_mse[i+1,3]<-sample_mse[i+1,1]
      sample_mse[i+1,4]<-sample_mse[i+1,2]
    }

  }
# Rownames
  if (!mba)
  {
    dim_names<-c("True MSE",paste("Lambda=",lambda_vec,", eta=",eta_vec,sep=""))
    if (prod(lambda_vec*eta_vec)==0)
    {
      dim_names[1+which(lambda_vec==0&eta_vec==0)]<-"DFA-MSE"
    }
  } else
  {
    dim_names<-c("MBA-MSE",paste("Lambda=",lambda_vec,", eta=",eta_vec,sep=""))
  }

  dimnames(xf)[[2]]<-dim_names
  dimnames(amp)[[2]]<-dim_names
  dimnames(shift)[[2]]<-dim_names
  dimnames(b)[[2]]<-dim_names
##
## We collect all frequency-domain performance statistics in corresponding tables
#  if (!mba)
#  {
## When using the periodogram we compute also curvature and peak-correlation in-sample because
    amp_shift_mat_insamp<-cbind(selectivity_insamp,mean_shift_insamp,
    curvature_insamp,lag_max_insamp,sample_mse[,c(1,3)])
    dimnames(amp_shift_mat_insamp)[[2]]<-c("Selectivity","Mean-shift",
    "Curv-in","Peak-Cor-in","MSE-in","scaled MSE-in")
    dimnames(amp_shift_mat_insamp)[[1]]<-dim_names
#  } else
#  {
#    amp_shift_mat_insamp<-cbind(selectivity_insamp,mean_shift_insamp)
#    dimnames(amp_shift_mat_insamp)[[2]]<-c("Selectivity","Mean-shift")
#    dimnames(amp_shift_mat_insamp)[[1]]<-dim_names
#  }
  amp_shift_mat_outsamp<-cbind(curvature_outsamp,lag_max_outsamp,sample_mse[,c(2,4)])
  dimnames(amp_shift_mat_outsamp)[[1]]<-dim_names
  if (mba)
  {
    if (estim_MBA)
    {
      dimnames(amp_shift_mat_outsamp)[[2]]<-c("Curv.out","Peak-Cor.out","MSE-out","Scaled sample MSE out")
    } else
    {
      dimnames(amp_shift_mat_outsamp)[[2]]<-c("Curv.","Peak-Cor.","Sample MSE","scaled Sample MSE")
    }
  } else
  {
    dimnames(amp_shift_mat_outsamp)[[2]]<-c("Curv.out", "Peak-Cor.out","MSE-out","Scaled sample MSE out")
    dimnames(amp_shift_mat_outsamp)[[1]]<-dim_names
  }

  ats_mat<-cbind(Accuracy,Timeliness,Smoothness,Residual,
  Accuracy+Timeliness+Smoothness+Residual)
  dimnames(ats_mat)[[2]][dim(ats_mat)[2]]<-"Total MSE"
  dimnames(ats_mat)[[1]]<-dim_names

  return(list(xf=xf,amp_shift_mat_insamp=amp_shift_mat_insamp,
  amp_shift_mat_outsamp=amp_shift_mat_outsamp,
  ats_mat=ats_mat,amp=amp,shift=shift,b=b,xf_sym=xf_sym,weight_func=weight_func))
}




for_sim_out<-function(a_vec,len1,len,cutoff,L,mba,estim_MBA,L_sym,Lag,i1,i2,scaled_ATS,lambda_vec,eta_vec,anzsim,M,dif)
{

  amp_sim<-array(dim=c(len/2+1,length(lambda_vec)+1,length(a_vec),anzsim))
  shift_sim<-amp_sim
  amp_shift_mat_sim<-array(dim=c(length(lambda_vec)+1,10,length(a_vec),anzsim))
  amp_shift_mat_e<-array(dim=c(length(lambda_vec)+1,10,length(a_vec)))
  amp_e<-array(dim=c(len/2+1,length(lambda_vec)+1,length(a_vec)))
  shift_e<-amp_e
  weight_func_e<-array(dim=c(len/2+1,length(a_vec)))
  xff_sim<-array(dim=c(len1,length(lambda_vec)+1,length(a_vec),anzsim))
  xff_sim_e<-array(dim=c(len1,length(lambda_vec)+1,length(a_vec)))
  xff_sim_sym<-array(dim=c(len1,length(lambda_vec)+1,length(a_vec),anzsim))
  xff_sim_sym_e<-array(dim=c(len1,length(lambda_vec)+1,length(a_vec)))
  ats_sym_e<-array(dim=c(length(lambda_vec)+1,5,length(a_vec)))
  ats_sym<-array(dim=c(length(lambda_vec)+1,5,length(a_vec),anzsim))
  amp_shift_mat_eper<-amp_shift_mat_e
  amp_sim_per<-array(dim=c(len/2+1,length(lambda_vec)+1,length(a_vec),anzsim))
  shift_sim_per<-amp_sim_per
  weight_func_sim_per<-array(dim=c(len/2+1,length(a_vec),anzsim))
  amp_eper<-array(dim=c(len/2+1,length(lambda_vec)+1,length(a_vec)))
  shift_eper<-amp_eper
  b_sim<-array(dim=c(L,length(eta_vec)+1,length(a_vec),anzsim))
  b_array<-array(dim=c(L,length(eta_vec)+1,length(a_vec)))

  for (i in 1:anzsim) #i<-1
  {
    for (ki in 1:length(a_vec))    #ki<-1
    {
  # The same seed is used for all three processes
      set.seed(i)

    # The longer sample is used for implementing the symmetric filter and
    # for computing peak-correlations
      if (dif)
      {
        x1<-diff(arima.sim(list(ar=a_vec[ki]),n=len1+1))
      } else
      {
        x1<-arima.sim(list(ar=a_vec[ki]),n=len1)
      }
  # In-sample data (first observations are discarded in order to initialize filters)
      x<-x1[max(L_sym/2,L)+1:len]
    # Length of symmetric target filter
    # model-based spectrum
      if (!mba)
      {
        plot_T<-F
        weight_func<-per(x,plot_T)$per          #ts.plot(weight_func)
        omega_k<-(0:(len/2))*pi/(len/2)
        Gamma<-(0:(len/2))<=((cutoff/pi)*len/2)

      } else
      {
        if (estim_MBA)
        {
          a1<-arima(x,order=c(1,0,0),include.mean=F)$coef
  # Ensure stability of ar(1)-spectrum
          a1<-min(0.99,a1)
        } else
        {
          a1<-a_vec[ki]
        }
        omega_k<-(0:(M))*pi/(M)
        weight_func<-1/abs(1-a1*exp(1.i*omega_k))^2
        Gamma<-(0:(M))<=((cutoff/pi)*M)

      }

    # Compute in/out-of-sample performances
      perf<-Performance_func(lambda_vec,eta_vec,len1,len,x1,weight_func,Gamma,cutoff,L,L_sym,a_vec[ki],mba,
                             scaled_ATS,estim_MBA,Lag,i1,i2)

      xf<-perf$xf
      xf_sym<-perf$xf_sym
      xff_sim_e[,,ki]<-xf
      xff_sim_sym_e[,,ki]<-xf_sym
    # Performance measures
      amp_shift_mat_e[,,ki]<-cbind(perf$amp_shift_mat_insamp,perf$amp_shift_mat_outsamp)       #amp_shift_mat_sim[,,ki,i]
    # Amplitude and shifts
      amp_e[,,ki]<-perf$amp       #amp_sim[,,ki,i]
      shift_e[,,ki]<-perf$shift       #amp_sim[,,ki,i]
      weight_func_e[,ki]<-perf$weight_func
      ats_sym_e[,,ki]<-perf$ats_mat
      b_array[,,ki]<-perf$b
    }
    dimnames(amp_shift_mat_e)[[1]]<-dimnames(cbind(perf$amp_shift_mat_insamp,perf$amp_shift_mat_outsamp))[[1]]
    dimnames(amp_shift_mat_e)[[2]]<-dimnames(cbind(perf$amp_shift_mat_insamp,perf$amp_shift_mat_outsamp))[[2]]
    dimnames(ats_sym_e)[[1]]<-dimnames(perf$ats_mat)[[1]]
    dimnames(ats_sym_e)[[2]]<-dimnames(perf$ats_mat)[[2]]
    amp_shift_mat_sim[,,,i]<-amp_shift_mat_e
    amp_sim_per[,,,i]<-amp_e
    shift_sim_per[,,,i]<-shift_e
    weight_func_sim_per[,,i]<-weight_func_e
    xff_sim[,,,i]<-xff_sim_e
    xff_sim_sym[,,,i]<-xff_sim_sym_e
    ats_sym[,,,i]<-ats_sym_e
    b_sim[,,,i]<-b_array
    dimnames(b_sim)[[2]]<-dimnames(perf$ats_mat)[[1]]
    dimnames(b_sim)[[1]]<-paste("Lag ",0:(L-1),sep="")
  }
  dim_names<-dimnames(amp_shift_mat_e)
  dimnames(amp_shift_mat_sim)[[1]]<-dim_names[[1]]
  dimnames(amp_shift_mat_sim)[[2]]<-dim_names[[2]]
  dimnames(ats_sym)[[1]]<-dimnames(ats_sym_e)[[1]]
  dimnames(ats_sym)[[2]]<-dimnames(ats_sym_e)[[2]]

  return(list(amp_shift_mat_sim=amp_shift_mat_sim,amp_sim_per=amp_sim_per, shift_sim_per= shift_sim_per,
              xff_sim=xff_sim,xff_sim_sym=xff_sim_sym,ats_sym=ats_sym,dim_names=dim_names,
              weight_func_sim_per=weight_func_sim_per,Gamma=Gamma,b_sim=b_sim,Gamma=Gamma,x1=x1))
}





#dim(xff_sim)
#dim(sample_id_out[[i]]$xff_sim_e)
#amp_shift_mat_sim[2,8,2,]-perf_out_sample[,2,1]   median(amp_shift_mat_sim[2,8,3,92])
#perf_out_sample[92,2,1]



#filter_output<-filter_output_out_sample
#ts.plot(filter_output,lty=1:4)
# ts.plot(Arg(dfa$trffkt)/(pi*(0:(length(dfa$trffkt)-1))/((length(dfa$trffkt)-1))))







perf_func_c_pc_mse<-function(filter_output)
{
  filter_output_second_diff<-apply(apply(filter_output[,1:(ncol(filter_output)-1)],2,diff),2,diff)

  curvature<-apply(filter_output_second_diff^2,2,mean)/(apply(filter_output[,1:(ncol(filter_output)-1)],2,var))
  #  names(curvature)<-c(paste("DFA-cust(",lambda_vec,",",eta_vec,")",sep=""),"MDFA leading indicator")
  anf_cor<--10
  enf_cor<-10
  lag_cor<-anf_cor:enf_cor
  peak_cor<-matrix(1:(length(lag_cor)*(length(curvature))),nrow=length(lag_cor),ncol=length(curvature))
  i_cor<-0
  for (j in lag_cor)   #j<-lag_cor[21]
  {
    i_cor<-i_cor+1
    for (i_cust in 1:(length(curvature)))#i_cust<-1
    {
      if (j<=0)
      {
        peak_cor[i_cor,i_cust]<-cor(filter_output[(1-j):len,ncol(filter_output)],
                                    filter_output[1:(len+j),i_cust])
      } else
      {
        peak_cor[i_cor,i_cust]<-cor(filter_output[1:(len-j),ncol(filter_output)],
                                    filter_output[(1+j):len,i_cust])
      }
    }
  }
  peak_lag<-1:(length(curvature))
  for (i_cust in 1:(length(curvature)))
  {
    peak_lag[i_cust]<-which(peak_cor[,i_cust]==max(peak_cor[,i_cust]))-1+min(0,anf_cor)
  }
  mse<-apply((filter_output[,1:(ncol(filter_output)-1)]-filter_output[,ncol(filter_output)])^2,2,mean)
  perf_mat_c_li<-rbind(curvature,peak_lag,mse)
  dimnames(perf_mat_c_li)[[1]]<-c("Curvature","Peak correlation","MSE")
  dimnames(perf_mat_c_li)[[2]]<-names(curvature)
  return(list(perf_mat_c_li=perf_mat_c_li))
}






mdfa_mse_leading_indicator_vs_dfa_customized<-function(anzsim,a1,cutoff,L,lambda_vec,
                                                       eta_vec,len1,len,i1,i2,Lag,lambda_mdfa,eta_mdfa,troikaner)
{

  lenh<-len1
  b_mat<-vector(mode="list")
  mdfa_list<-vector(mode="list")
  weight_func<-matrix(rep(0:(len/2),3),ncol=3)
  xh<-1:len1
  y<-yhat<-x<-1:len
  perf_out_sample<-perf_in_sample<-array(dim=c(anzsim,3,length(lambda_vec)+length(lambda_mdfa)))
  ord<-1000
  # Filter weights ideal trend (See DFA)
  gamma<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))

# MSE-settings for DFA and MDFA
  d<-0
  lin_eta<-F
  weight_constraint<-rep(1/(ncol(weight_func)-1),ncol(weight_func)-1)
  lambda_cross<-lambda_smooth<-0
  lambda_decay<-c(0,0)
  lin_expweight<-F
  shift_constraint<-rep(0,ncol(weight_func)-1)
  grand_mean<-F
  b0_H0<-NULL
  c_eta<-F
  weights_only<-F
  weight_structure<-c(0,0)
  white_noise<-F
  synchronicity<-F
# 26.11
  lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)

  for (i_loop in 1:anzsim)#i_loop<-100
  {

# Generate series
    set.seed(i_loop)
    xh<-arima.sim(list(ar=a1),n=len1)

# Use another set.seed in order to avoid colinearity of the idiodyncratic component
    set.seed(1+i_loop)
# Scaling of the idiosyncratic noise
    scale_idiosyncratic<-0.1
    eps<-rnorm(length(xh))
    indicator<-xh+scale_idiosyncratic*eps
# Data: first column=target, second column=x, third column=shifted (leading) indicator
    data_matrix<-cbind(xh,xh,c(indicator[2:length(xh)],NA))
    dimnames(data_matrix)[[2]]<-c("target","x","leading indicator")
# Extract 120 observations from the long sample
    data_matrix_120<-data_matrix[len1/2+(-len/2):((len/2)-1),]

    insample<-nrow(data_matrix_120)
# Compute DFT
    weight_func<-spec_comp(insample, data_matrix_120, d)$weight_func
    K<-nrow(weight_func)-1
# Target
    Gamma<-(0:K)<=as.integer(cutoff*K/pi)+1
# Estimate MDFA MSE filter coefficients  #Lag<-2
    for (i_mdfa in 1:length(lambda_mdfa))#i_mdfa<-1  ts.plot(cust_leading_obj$mdfa_list[[i_mdfa]]$b)
    {

      mdfa_obj<-mdfa_analytic(L,lambda_mdfa[i_mdfa],weight_func,Lag,Gamma,eta_mdfa[i_mdfa],cutoff,i1,i2,weight_constraint,
                            lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,b0_H0,
                            c_eta,weight_structure,white_noise,synchronicity,lag_mat,troikaner)
      mdfa_list[[i_mdfa]]<-mdfa_obj


# Filter coefficients
      b_mat[[i_mdfa]]<-mdfa_obj$b#ts.plot(mdfa_obj$b)
      dimnames(b_mat[[i_mdfa]])[[2]]<-c("x","leading indicator")
      dimnames(b_mat[[i_mdfa]])[[1]]<-paste("Lag ",0:(L-1),sep="")#dim(b_mat)
    }
# DFA
    periodogram<-abs(weight_func[,1])^2
# Estimate customized coefficients
    for (i_cust in 1:length(lambda_vec))#i_cust<-1
    {
      dfa<-dfa_analytic(L,lambda_vec[i_cust],periodogram,Lag,Gamma,eta_vec[i_cust],cutoff,i1,i2)
      if (i_cust==1)
      {
        b_cust<-dfa$b
      } else
      {
        b_cust<-cbind(b_cust,dfa$b)
      }
    }

    y_insample<-rep(NA,len)
    y_outsample<-rep(NA,len)
# In order to replicate code in McElroy-Wildi we set all coefficients to zero for lags>470 (reconciliation)
    gamma[471:length(gamma)]<-0
# Compute the outputs yt of the (truncated) symmetric target filter
    for (j in 1:120)# j<-120
    {
      y_insample[j]<-gamma[1:700]%*%data_matrix[len1/2+(-len/2)-1+(j:(j-699)),2]+gamma[2:700]%*%data_matrix[len1/2+(-len/2)+(j:(j+698)),2]
      y_outsample[j]<-gamma[1:700]%*%data_matrix[len+len1/2+(-len/2)-1+(j:(j-699)),2]+gamma[2:700]%*%data_matrix[len+len1/2+(-len/2)+(j:(j+698)),2]
    }

    y_target_leading_indicator_in_sample<-y_insample
    y_target_leading_indicator_out_sample<-y_outsample

    yhat_multivariate_leading_indicator_in_sample<-yhat_multivariate_leading_indicator_out_sample<-matrix(nrow=len,ncol=length(lambda_mdfa))
    for (i_mdfa in 1:length(lambda_mdfa))
    {
      for (j in 1:len)
      {
        yhat_multivariate_leading_indicator_in_sample[j,i_mdfa]<-sum(apply(b_mat[[i_mdfa]]*data_matrix[lenh/2+(-len/2)-1+j:(j-L+1),2:3],1,sum))
        yhat_multivariate_leading_indicator_out_sample[j,i_mdfa]<-sum(apply(b_mat[[i_mdfa]]*data_matrix[len+lenh/2+(-len/2)-1+j:(j-L+1),2:3],1,sum))
      }
    }
    yhat_cust_out_sample<-yhat_cust_in_sample<-matrix(nrow=len,ncol=length(lambda_vec))
    for (i_cust in 1:length(lambda_vec))
    {
      for (j in 1:len)
      {
        yhat_cust_in_sample[j,i_cust]<-b_cust[,i_cust]%*%data_matrix[lenh/2+(-len/2)-1+j:(j-L+1),2]
        yhat_cust_out_sample[j,i_cust]<-b_cust[,i_cust]%*%data_matrix[len+lenh/2+(-len/2)-1+j:(j-L+1),2]
      }
    }
    filter_output_in_sample<-cbind(yhat_cust_in_sample,yhat_multivariate_leading_indicator_in_sample,
                                   y_target_leading_indicator_in_sample)

    dimnames(filter_output_in_sample)[[2]]<-c(paste("DFA(",lambda_vec,",",eta_vec,")",sep=""),
                                              paste("MDFA(",lambda_mdfa,",",eta_mdfa,")",sep=""),"Target")
    filter_output_out_sample<-cbind(yhat_cust_out_sample,yhat_multivariate_leading_indicator_out_sample,
                                    y_target_leading_indicator_out_sample)
    dimnames(filter_output_out_sample)[[2]]<-dimnames(filter_output_in_sample)[[2]]
    head(filter_output_in_sample)
# In sample performances
    perf_obj<-perf_func_c_pc_mse(filter_output_in_sample)
    perf_in_sample[i_loop,,]<-perf_obj$perf_mat_c_li
# out-of-sample performances
    perf_obj_out<-perf_func_c_pc_mse(filter_output_out_sample)
    perf_out_sample[i_loop,,]<-perf_obj_out$perf_mat_c_li
  }
  return(list(perf_out_sample=perf_out_sample,perf_in_sample=perf_in_sample,filter_output_in_sample=filter_output_in_sample,filter_output_out_sample=filter_output_out_sample,mdfa_list=mdfa_list,weight_func=weight_func))
}

#perf_out_sample[,2,1]  median(perf_out_sample[,2,1])


