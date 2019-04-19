


tic_simulation_experiment_func<-function(a1,b1,true_model_order,lambda_decay,lambda_smooth,anzsim)
{
  
  
  # We generate anzsim realizations of length len of the arma-process
  set.seed(1)
  len<-1000
  mse_true_arma<-mse_dfa<-NULL
  # Length of in-sample span
  #   For illustration we here use a short sample
  #   In fact regularization is deemed useful when working with small samples (otherwise unconstrained designs work fine)
  in_sample<-100
  # Frequency grid for DFA based on true model
  K_true<-600
  # Lowpass target  
  periodicity<-10
  # Deliberately exagerate overfitting: we use filters of length 4*periodicity 
  #   -This filter length is used for all designs
  #   -Expected effect: when imposing regularization in short samples the performances improve out-of-sample
  L<-4*periodicity
  #   MSE design (no customization)
  lambda<-eta<-0
  # Cross sectional regularization is not activated here since design is univariate  
  lambda_cross<-0 
  # Nowcast
  Lag<-0
  # Length of ideal filter
  M<-100
  mse_true<-mse_dft<-mse_reg_dft<-tic_by_hand<-tic<-degrees_of_freedom<-NULL
  pb <- txtProgressBar(min = 1, max = anzsim, style = 3)
  # Loop through all simulations and collect out-of-sample forecast performances
  for (i in 1:anzsim)#i<-1
  {
    # Distinguish white noise  
    if (abs(a1)+abs(b1)>0)
    {
      # Generate series  
      x<-as.vector(arima.sim(n=len,list(ar=a1,ma=b1)))
    } else
    {
      x<-rnorm(len)
    }
    # Use in-sample span for model-estimation and for dft  
    x_insample<-x[1:in_sample]
    # True model: estimate model-parameters by relying on classic arima-function
    arima_true_obj<-arima(x_insample,order=true_model_order,include.mean=F)
    # Spectrum based on true model  
    spec<-arma_spectrum_func(ifelse(!is.na(arima_true_obj$coef["ar1"]),arima_true_obj$coef["ar1"],0),ifelse(!is.na(arima_true_obj$coef["ma1"]),arima_true_obj$coef["ma1"],0),K_true,F)$arma_spec
    weight_func<-cbind(spec,spec)
    colnames(weight_func)<-c("spectrum target","spectrum explanatory")
    weight_func_true<-weight_func
    cutoff<-pi/periodicity
    # target true model: frequnecy grid is not the same as for dft below i.e. K is different
    Gamma_true<-(0:(K_true))<=K_true*cutoff/pi+1.e-9
    mdfa_true_obj<-MDFA_mse(L,weight_func_true,Lag,Gamma_true)$mdfa_obj 
    b_true<-mdfa_true_obj$b
    # Use in-sample span for dft  
    weight_func_dft<-cbind(per(x_insample,F)$DFT,per(x_insample,F)$DFT)
    colnames(weight_func_dft)<-c("spectrum target","spectrum explanatory")
    K_dft<-nrow(weight_func_dft)-1
    Gamma_dft<-(0:(K_dft))<=K_dft*cutoff/pi+1.e-9
    # Compute unconstrained MSE-filter (same as in example 3 of tutorial 3)
    mdfa_dft_obj<-MDFA_mse(L,weight_func_dft,Lag,Gamma_dft)$mdfa_obj 
    b_dft<-mdfa_dft_obj$b
    # New advanced feature: we set troikaner to T in order to calculate all relevant additional statistics  
    troikaner<-T
    mdfa_reg_obj<-MDFA_reg(L,weight_func_dft,Lag,Gamma_dft,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth,troikaner)$mdfa_obj
    b_dft_reg<-mdfa_reg_obj$b
    # Here we compute the new information criterion 'by hand'
    #   -One can show that the penalty term 2*mdfa_reg_obj$edof/(2*K_dft)) is larger than the 
    #    spurious decrease of the (log of the) MSE-criterion log(mdfa_reg_obj$MS_error) when overparametrization is attained, see up-coming book
    tic_by_hand<-c(tic_by_hand,log(mdfa_reg_obj$MS_error)+2*mdfa_reg_obj$edof/(2*K_dft))
    # Alternatively the criterion is computed explicitly in the MDFA function  
    tic<-c(tic,mdfa_reg_obj$tic)
    degrees_of_freedom<-c(degrees_of_freedom,mdfa_reg_obj$edof)
    # Filter data
    # 1. ideal filter  
    id_obj<-ideal_filter_func(periodicity,M,x)
    output_ideal<-id_obj$y
    # 2. DFA true
    filt_true_obj<-filt_func(x,b_true)
    output_dfa_true<-filt_true_obj$yhat
    # 3. DFA dft
    filt_dft_obj<-filt_func(x,b_dft)
    output_dfa_dft<-filt_dft_obj$yhat
    # 4. DFA dft regularized
    filt_dft_reg_obj<-filt_func(x,b_dft_reg)
    output_dfa_dft_reg<-filt_dft_reg_obj$yhat
    # Mean-square out-of-sample filter error  
    mse_true<-c(mse_true,mean((output_ideal-output_dfa_true)[(in_sample+1):len-M]^2,na.rm=T))
    mse_dft<-c(mse_dft,mean((output_ideal-output_dfa_dft)[(in_sample+1):len-M]^2,na.rm=T))
    mse_reg_dft<-c(mse_reg_dft,mean((output_ideal-output_dfa_dft_reg)[(in_sample+1):len-M]^2,na.rm=T))
    setTxtProgressBar(pb, i)
  }
  

  # Compute the ratio of root mean-square forecast errors:
  #   The ratio cannot be larger than 1 asymptotically because our particular design distinguishes arma as the universally best possible design
  result_mat<-matrix(c(sqrt(mean(mse_true)/mean(mse_dft)),sqrt(mean(mse_true)/mean(mse_reg_dft)),mean(tic),round(mean(degrees_of_freedom),0)),nrow=1,ncol=4)
  colnames(result_mat)<-c("Unconstrained","Regularized","tic","degrees of freedom")
  return(list(result_mat=result_mat))
  
}