
jags.piecewise.wei_chng <-"

data{
  for(i in 1:N){
    zero[i] <- 0
  }
  zero_prior <- 0
}

model {
  #Constant for zerors trick
  C <- 10000
  
  #Prior for change-point locations - See Chapell
  cp[1] = 0  # mcp helper value. Should be zero
  cp[N_CP+2] = 9999  # mcp helper value. very large number
  
  for(k in 1:N_CP){
    unif[k] ~ dunif(0.001, min(max(time),MAXX))
  }
  cp_x[2:(N_CP+1)] <- sort(unif)
  
  for(i in 2:(N_CP+1)){
    cp[i] <- cp_x[i]*equals(cp_fixed,0) + cp_fix[i]*equals(cp_fixed,1)
  }
  
  for(i in 2:(N_CP+1)){
    diff[i-1] <- cp[i]-cp[i-1]
  }
  
  log_prior <- loggam(2*N_CP +2) + sum(log(diff)) - ((2*N_CP +1)*log(MAXX))
  zero_prior ~ dpois(C - log_prior)
  
  #Prior for the model parameters
  
  for(i in 1:2){
    sd[i]  <- 1 # ~ dunif(0.2,5)
    prec[i] <- pow(sd[i], -2)
  }
  
  for(k in 1:(N_CP+1)){
    for(j in 1:ncovar_cp){ #Always the number of basic parameters plus 1 for a covariate i.e. one treatment
      beta_cp_x[k,j] ~ dnorm(0,prec[1])
      beta_cp_anc_x[k,j] ~ dnorm(0,prec[1])
    }
  }
  
  for(k in 1:(N_CP+1)){
  beta_cp[k,1] <-  beta_cp_x[k,1]
  beta_cp_anc[k,1] <- beta_cp_anc_x[k,1]
  
    for(j in 2:ncovar_cp){
      beta_cp[k,j] <-  beta_cp_x[k,j]*equals(scenario_common[k],0) 
      beta_cp_anc[k,j] <- beta_cp_anc_x[k,j]*equals(scenario_common[k],0)
    }
  }
  
  #Typically all other ancillary parameters e.g shape for weibull will only vary by 
  #time and not be subject to covariates
  
  #For Factors that do not change with respect to time or change-point
  for(i in 1:ncovar_fixed){ #Will always be 2 just as a dummy with matrix of 0,0
    beta_covar_fixed[i] ~ dnorm(0,prec[2])
  }
  #Even if there are no covariates we will just initalize it as a matrix of 0,0
  
  
  linpred_fixed <- X_mat_fixed %*%beta_covar_fixed
  
  # Model and likelihood
  for (i in 1:N) {
    
    for(k in 1:(N_CP+1)){
      #variable which gives the difference between the two intervals if time[i]>cp[k+1]
      #(i.e. cp[k+1] -cp[k]) or time between time[i] and cp[k]
      X[i,k] = max(min(time[i], cp[k+1]) - cp[k],0) 
      
      #Indicator variable which highlights which interval time is in 
      X_ind[i,k] = step(time[i]-cp[k])*step(cp[k+1]-time[i])
      
      for(p in 1:ncovar_cp){
        linpred_cp_param[i,k,p] <- X_mat_trt[i,p]*beta_cp[k,p]
      }
      linpred_cp[i,k] <- sum(linpred_cp_param[i,k,])
      param_1[i,k] <- exp(linpred_cp[i,k] + linpred_fixed[i]) #lambda
      
    #param_2[i,k] <- exp(beta_cp_anc[k,1]) #Different shapes but common for both arms
    #param_2[i,k] <- exp(beta_cp_anc[1,1]) #Can be made constant across all if k = 1
    #param_2[i,k] <- exp(beta_cp_anc[k,1] + equals(k,1)*beta_cp_anc[k,2]*X_mat_trt[i,2]) #We allow an independent model for the first time-point i.e. different shapes for both treatments.
    #param_2[i,k] <- if scenario is 4 then param is zero i.e. exponential.
    
    param_2[i,k] <- exp((beta_cp_anc[k,1]*equals(scenario,1)) + 
                        (beta_cp_anc[1,1]*equals(scenario,2)) + 
                        (beta_cp_anc[k,1]+equals(k,1)*beta_cp_anc[k,2]*X_mat_trt[i,2])*equals(scenario,3))


      log_haz_seg[i,k] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))*X_ind[i,k] 
      cum_haz_seg[i,k] <- param_1[i,k]*pow(X[i,k]+cp[k],param_2[i,k]) -  param_1[i,k]*pow(cp[k],param_2[i,k])
      
    }
    
    log_haz[i] <- sum(log_haz_seg[i,])
    cum_haz[i] <- sum(cum_haz_seg[i,])
    
    loglik[i] = status[i]*log_haz[i] - cum_haz[i]
    zero[i] ~ dpois(C - loglik[i])
    
  }
  
  #Predict Survival
  for(i in 1:length(t_pred)){
    
    for(k in 1:(N_CP+1)){
      
      
      X_pred[i,k] = max(min(t_pred[i], cp[k+1]) - cp[k],0) 
      X_ind_pred[i,k] = step(t_pred[i]-cp[k])*step(cp[k+1]-t_pred[i])
      
      for(q in 1:length(index_pred)){
        
        for(p in 1:ncovar_cp){
          linpred_cp_param_pred[i,k,q,p] <- X_mat_trt[index_pred[q],p]*beta_cp[k,p]
        }
        linpred_cp_pred[i,k,q] <- sum(linpred_cp_param_pred[i,k,q,])
        
        param_1_pred[i,k,q] <- exp(linpred_cp_pred[i,k,q]  + linpred_fixed[index_pred[q]]) #lambda
        param_2_pred[i,k,q] <- exp((beta_cp_anc[k,1]*equals(scenario,1)) + 
                                   (beta_cp_anc[1,1]*equals(scenario,2)) + 
                                   (beta_cp_anc[k,1]+equals(k,1)*beta_cp_anc[k,2]*X_mat_trt[index_pred[q],2])*equals(scenario,3))

        log_haz_seg_pred[i,k,q] <-  log(param_2_pred[i,k,q]*param_1_pred[i,k,q]*pow(t_pred[i],param_2_pred[i,k,q]-1))*X_ind_pred[i,k] 
        cum_haz_seg_pred[i,k,q] <- param_1_pred[i,k,q]*pow(X_pred[i,k]+cp[k],param_2_pred[i,k,q]) -  param_1_pred[i,k,q]*pow(cp[k],param_2_pred[i,k,q])
        
      }
      
    }
    
    for(q in 1:length(index_pred)){
      cum_haz_pred[i,q] <- sum(cum_haz_seg_pred[i,,q])
    }
  } 
  
  
  total_loglik <- sum(loglik)
}

"

jags.piecewise.wei_chng_wane <-"

  data{
    for(i in 1:N){
    zero[i] <- 0
    }
    zero_prior <- 0
  }

  model {
  #Constant for zerors trick
  C <- 10000

 #Prior for change-point locations - See Chapell
  cp[1] = 0  # mcp helper value. Should be zero
  cp[N_CP+2] = 9999  # mcp helper value. very large number

  for(k in 1:N_CP){
 unif[k] ~ dunif(0.001, max(time))
 
 #unif[k] ~ dunif(0.001,max(time)*2)
 
  }
  cp[2:(N_CP+1)] <- sort(unif)

  for(i in 2:(N_CP+1)){
  diff[i-1] <- cp[i]-cp[i-1]
  }

  log_prior <- loggam(2*N_CP +2) + sum(log(diff)) - ((2*N_CP +1)*log(MAXX))
  zero_prior ~ dpois(C - log_prior)


  #Prior for the model parameters
  
  for(i in 1:2){
    sd[i]  <- 1 # ~ dunif(0.2,5)
    prec[i] <- pow(sd[i], -2)
  }

  for(k in 1:(N_CP+1)){
    for(j in 1:2){ #Always the number of basic parameters plus 1 for a covariate 
      beta_cp[k,j] ~ dnorm(0,prec[1])
      beta_cp_anc[k,j] ~ dnorm(0,prec[1])
    }
    
  }
  
  #Typically all other ancillary parameters e.g shape for weibull will only vary by 
  #time and not be subject to covariates
  
  #For Factors that do not change with respect to time or change-point
  for(i in 1:ncovar_fixed){ #Will always be 2 just as a dummy with matrix of 0,0
  beta_covar_fixed[i] ~ dnorm(0,prec[2])
  }
  #Even if there are no covariates we will just initalize it as a matrix of 0,0
  

  linpred_fixed <- X_mat_fixed %*%beta_covar_fixed

  # Model and likelihood
  for (i in 1:N) {
  
  for(k in 1:N_CP){
    #variable which gives the difference between the two intervals if time[i]>cp[k+1]
    #(i.e. cp[k+1] -cp[k]) or time between time[i] and cp[k]
    X[i,k] = max(min(time[i], cp[k+1]) - cp[k],0) 
   
    #Indicator variable which highlights which interval time is in 
    X_ind[i,k] = step(time[i]-cp[k])*step(cp[k+1]-time[i])
    
    for(p in 1:ncovar_cp){
      linpred_cp_param[i,k,p] <- X_mat_trt[i,p]*beta_cp[k,p]
    }
    linpred_cp[i,k] <- sum(linpred_cp_param[i,k,])
    param_1[i,k] <- exp(linpred_cp[i,k] + linpred_fixed[i]) #lambda
    param_2[i,k] <- exp(beta_cp_anc[k,1]*equals(scenario,1)) #Different shapes but common for both arms

    log_haz_seg[i,k] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))*X_ind[i,k] 
    cum_haz_seg[i,k] <- param_1[i,k]*pow(X[i,k]+cp[k],param_2[i,k]) -  param_1[i,k]*pow(cp[k],param_2[i,k])
    
  }
  
  
  for(k in (N_CP+1):(N_CP+1)){ #Treatment waning interval - Last interval
   
    X[i,k] = max(min(time[i], cp[k+1]) - cp[k],0) 
    X_ind[i,k] = step(time[i]-cp[k])*step(cp[k+1]-time[i])
    
    #Baseline Hazard Parameter --- Assume only 2 treatments
    linpred_cp_param[i,k,1] <- X_mat_trt[i,1]*beta_cp[k,1] 
    
    #Treatment Waning - Parameter
    #Beta_cp;  trt waning parameter
    initial_HR[i] <-   exp(beta_cp[k-1,2])
    lambda_wane[i] <- exp(beta_cp[k,2])

    HR_wane[i] <- 1-(1-initial_HR[i])*exp(-lambda_wane[i]*X[i,k])
    linpred_cp_param[i,k,2] <- 0 # This is redundant - Integral will evaluate the treatment effect
    
    
    linpred_cp[i,k] <- sum(linpred_cp_param[i,k,])
    param_1[i,k] <- exp(linpred_cp[i,k] + linpred_fixed[i]) #lambda
    param_2[i,k] <- exp(beta_cp_anc[k,1]*equals(scenario,1)) #If scenario is 999 then exponential 
    
    #log_haz_stn[i] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))*X_ind[i,k] 
    #log_haz_wane[i] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))*X_ind[i,k] 
    
    #(log(haz)+log(haz_wane)) == log(haz*haz_waner)
    log_haz_seg[i,k] <-  (log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1)) + log(HR_wane[i])*equals(X_mat_trt[i,2],1))*X_ind[i,k] 
   
    #Need to integrate this
    #Treatment waning... a is shape; m is scale
    #(amt^{a-1})*(1-(1-r)exp(-wt)
    #drop am -- need to add it back in at end
    #(t^{a-1})*(1-(1-r)exp(-w(t-q))
    #(1/lambda2-exp(-t*lambda2)/lambda2)*r*lambda


#incomplete_gamma_upper<- function(s,x){
#  exp(lgamma(s))*(1-pgamma(x,s,1))
#}

#int_1 <- (t_0^a)/a-(incomplete_gamma_upper(a,hr_wane_lambda*t_0)*(r-1)*exp(q*hr_wane_lambda)*(t_0^a))/((hr_wane_lambda*t_0)^a)
#int_2 <- (t_1^a)/a-(incomplete_gamma_upper(a,hr_wane_lambda*t_1)*(r-1)*exp(q*hr_wane_lambda)*(t_1^a))/((hr_wane_lambda*t_1)^a)

    upper_inc_gamma1[i] <- exp(loggam(param_2[i,k]))*(1-pgamma(lambda_wane[i]*time[i],param_2[i,k],1))
    upper_inc_gamma2[i] <- exp(loggam(param_2[i,k]))*(1-pgamma(lambda_wane[i]*cp[k],param_2[i,k],1))
    upper_int[i] <- (time[i]^param_2[i,k])/param_2[i,k]-(upper_inc_gamma1[i]*(initial_HR[i]-1)*exp(cp[k]*lambda_wane[i])*(time[i]^param_2[i,k]))/((lambda_wane[i]*time[i])^param_2[i,k])
    lower_int[i] <- (cp[k]^param_2[i,k])/param_2[i,k]-(upper_inc_gamma2[i]*(initial_HR[i]-1)*exp(cp[k]*lambda_wane[i])*(cp[k]^param_2[i,k]))/((lambda_wane[i]*cp[k])^param_2[i,k])

    cum_haz_stn[i] <- param_1[i,k]*pow(X[i,k]+cp[k],param_2[i,k]) -  param_1[i,k]*pow(cp[k],param_2[i,k])
    cum_haz_wane[i] <- ((upper_int[i] - lower_int[i])*param_2[i,k]*param_1[i,k])*step(time[i]-cp[k]) #Has to be treatment and cp < time
    cum_haz_seg[i,k] <- cum_haz_wane[i]*equals(X_mat_trt[i,2],1) +  cum_haz_stn[i]*equals(X_mat_trt[i,2],0)


 }

  log_haz[i] <- sum(log_haz_seg[i,])
  cum_haz[i] <- sum(cum_haz_seg[i,])
  
  loglik[i] = status[i]*log_haz[i] - cum_haz[i]
  zero[i] ~ dpois(C - loglik[i])

  }
   
  
  total_loglik <- sum(loglik)
  }

"





