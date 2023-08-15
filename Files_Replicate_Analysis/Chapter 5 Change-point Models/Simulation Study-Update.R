library("sfsmisc")
library("xtable")
pathway <- "~/Change-Point Simulation Studies/Simulation Study 2023/"
source(paste0(pathway,"Jags_codes.R" ))

#shape will be 0.7 and 1.2
#n.samp will be 500, 200 #and 100
#Simulation Study


weib <- TRUE

t_max <- 15
t_cens <- 4


n.samp <- c(200,500,1000)
shape_1 <- c(1.2,0.7)
scale_1 <- 0.3
initial_HR <- c(0.25,0.5, 0.75)

param_vals <- expand.grid(n.samp = n.samp, shape_1=shape_1,scale_1=scale_1, initial_HR =initial_HR)


n_sims <- 100

param_mods <- c("exp","weibullPH",  "lnorm", "gamma", 
                "gompertz", "llogis","gengamma", "rps")



n.samp <- 500
n.thin <- 1 
n.burnin <- 100


#Treatment Waning
incomplete_gamma_upper<- function(s,x){
  exp(lgamma(s))*(1-pgamma(x,s,1))
}
#Write it out in integral calc
cum_haz_seg <- function(t_0,t_1,a, m, initial_HR, lambda_wane  ){
  
  constant <- (initial_HR-1)*exp(t_0*lambda_wane)
  int_1_ind <- (t_0^a)/a-(incomplete_gamma_upper(a,lambda_wane*t_0)*(t_0^a)*constant)/((lambda_wane*t_0)^a)
  
  int_2 <- rep(NA, length(t_1))
  int_1 <- rep(int_1_ind, length(t_1))
  
  int_2 <- (t_1^a)/a-(incomplete_gamma_upper(a,lambda_wane*t_1)*(t_1^a)*constant)/((lambda_wane*t_1)^a)
  return((int_2-int_1)*a*m)
}


#Specific vals

if(FALSE){#Not running this scenario

df_mean_vals <- list()
for(c in 1:nrow(param_vals)){
  
  lambda_wane = 0.5
  cp <- 2
  time_1 <- seq(0, cp, by = 0.001)
  time_all <- seq(0, t_max, by = 0.001)
  t_0 <- cp
  t_1 = 5
  n.samp <- param_vals$n.samp[c]
  shape_1 <- param_vals$shape_1[c]
  scale_1 <- param_vals$scale_1[c]
  scale_2 <- param_vals$scale_1[c]*param_vals$initial_HR[c]

    
  #Important plot
  #t_vec <- 0:15
  #plot(t_vec, y = scale_1*(1-(1-initial_HR)*exp(-lambda_wane*(t_vec-cp))), xlab = "Time", ylab = "Hazard")
  #abline(v = cp)#, h = initial_HR*lambda)
  
  
  cum_Haz_baseline <- flexsurv:::HweibullPH(time_all,
                                            shape = shape_1,
                                            scale = scale_1 )
  
  cum_Haz_trt <- flexsurv:::HweibullPH(time_1,shape = shape_1, scale = scale_2)
  
  
  time_2 <- seq(cp, t_max, by = 0.001)
  time_2 <- tail(time_2, n = -1)
  cum_Haz_trt_cp <- cum_haz_seg(t_0,time_2,a = shape_1, m = scale_1, param_vals$initial_HR[c], lambda_wane)
  
  cum_haz_all_trt_cp <- (tail(cum_Haz_trt, n = 1)+cum_Haz_trt_cp)
  
  St_baseline <- exp(-cum_Haz_baseline)
  St_trt <- exp(-c(cum_Haz_trt,cum_haz_all_trt_cp))
  
  plot(x = time_all, y = St_baseline, ylim = c(0,1), type = "l", 
       ylab = "St", xlab = "time")
  lines(x =time_all,y = St_trt ,col = "blue")
  
  all_res <- list()
  
  for(s in 1:n_sims){
 
    n.samp_size <- param_vals$n.samp[c]
    sim_unif1 <- runif(n.samp_size)
    sim_unif1 <- sim_unif1[order(sim_unif1)]
    
    sim_unif2 <- sim_unif1#runif(n.samp)
    sim_unif2 <- sim_unif2[order(sim_unif2)]
    
    time_trt_samp <- time_all[sapply(sim_unif1, function(x){which.min(abs(St_trt-x))})]
    time_comp_samp <- time_all[sapply(sim_unif2, function(x){which.min(abs(St_baseline-x))})]
    
    df_all <- data.frame(time = c(time_trt_samp,time_comp_samp),
                         arm = c(rep(1, length(time_trt_samp)),
                                 rep(0, length(time_comp_samp)))) 
    
    df_all$time <- df_all$time +0.0001
    
    df_all <- df_all  %>%
      mutate(time_event = time,
             status_true = 1,
             status = ifelse(time > t_cens, 0,1),
             time = ifelse(time > t_cens, t_cens, time))
    
    
    
    km.all <- survfit(Surv(time, status)~arm, data = df_all)
    
    #plot(km.all)
    #lines(x =time_all,y = St_trt ,col = "blue")
    #lines(x =time_all,y = St_baseline ,col = "red")
    
    #coxph(Surv(time, status)~arm, data = df_all)
    
    data_jags <- list()
    data_jags$N <- nrow(df_all)
    data_jags$time <- df_all$time
    data_jags$status <- df_all$status
    data_jags$MAXX <- max(df_all$time)
    data_jags$N_CP <- 1
    data_jags$ncovar <- 2
    
    data_jags$N_CP <-1
    data_jags$ncovar_fixed <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix
    data_jags$ncovar_cp <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix
    
    data_jags$t_pred <- seq(0,t_max, by = 0.2)
    
    #X_mat_fixed is just a dummy matrix
    #PH Model with Changepoint
    data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N*2)),  ncol = data_jags$ncovar_fixed)
    #independent parameters 
    data_jags$X_mat_trt <- matrix(c(rep(1, data_jags$N),
                                    df_all$arm),  ncol = 2)
    
    data_jags$scenario <- 1
    
    #Change-inits
    inits <- function(){
      list(unif = runif(data_jags$N_CP,max = data_jags$MAXX/data_jags$N_CP)
           #,
           #beta_cp_0 = c(0),
           #beta_cp = rnorm(data_jags$N_CP+1)
      )
    }
  mod.cp <- runjags::run.jags(
      model = get(paste0("jags.piecewise.wei_chng_wane")),
      #model = jags.piecewise.expo_wane,
      data = data_jags,
      n.chains = 2,
      monitor = c("cp", "beta_cp", "total_loglik", "cum_haz",
                  "lower_int", "upper_int", "beta_cp_anc", "loglik"), 
      sample=n.samp, 
      thin = n.thin, 
      burnin = n.burnin,
      inits = inits, 
      method ='rjparallel')
    
    #add.summary(mod.cp, "cp")
    
    param_names <- colnames(mod.cp$mcmc[[1]]) #rownames(expo.mod_1chng$BUGSoutput$summary)
    mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])
    
    # loglik_mat <- mcmc_output[,grep("loglik",param_names)]
    # loglik_mat <- loglik_mat[,-grep("total_loglik",colnames(loglik_mat))]
    # waic <- loo::waic(loglik_mat)[["estimates"]][3, 1]
    # pml <- -2 * sum(log(nrow(loglik_mat)/colSums(1/exp(loglik_mat))))
    
    St_mat <- exp(-mcmc_output[,grep("cum_haz[",param_names, fixed = T)])
    St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))
    
    t_pred <- seq(0, t_max, by = 0.2)
    HR_mat <- St_pred2 <- St_pred <- matrix(ncol = length(t_pred),
                                            nrow = nrow(mcmc_output))
    #Fix plot
    for(i in 1:nrow(mcmc_output)){
      
      m1_1 <- exp(mcmc_output[i,"beta_cp[1,1]"] + mcmc_output[i,"beta_cp[1,2]"] )
      m1_2 <- exp(mcmc_output[i,"beta_cp[1,1]"])
      
      
        a1 <- exp(mcmc_output[i,"beta_cp_anc[1,1]"])
      
      St_pred_cp <- flexsurv::pweibullPH(mcmc_output[i,"cp[2]"], scale = m1_1, shape = a1, lower.tail = F, log = F)
      St_pred2_cp <- flexsurv::pweibullPH(mcmc_output[i,"cp[2]"], scale = m1_2, shape = a1, lower.tail = F, log = F)
      
      for(j in 1:length(t_pred)){
        if(t_pred[j] <= mcmc_output[i,"cp[2]"]){
          St_pred[i,j] <- flexsurv::pweibullPH(t_pred[j], scale = m1_1, shape = a1, lower.tail = F, log = F)
          St_pred2[i,j] <- flexsurv::pweibullPH(t_pred[j], scale = m1_2, shape = a1, lower.tail = F, log = F)
          HR_mat[i,j] <-  exp(mcmc_output[i,"beta_cp[1,2]"])
        }else{
          m2 <- exp(mcmc_output[i,"beta_cp[2,1]"])
          a2 <- exp(mcmc_output[i,"beta_cp_anc[2,1]"])
          
          X_eval <- t_pred[j] -mcmc_output[i,"cp[2]"]
          
          initial_HR_sim <-  exp(mcmc_output[i,"beta_cp[1,2]"])
          lambda_wane <- exp(mcmc_output[i,"beta_cp[2,2]"])
          HR_final <-1-(1-initial_HR_sim)*exp(-lambda_wane*X_eval)
          HR_mat[i,j] <-  HR_final
          
          t_0 <- mcmc_output[i,"cp[2]"]
          t_1 <- t_pred[j]
          
          St_pred[i,j] <- exp(-(-log(St_pred_cp) +cum_haz_seg(t_0 = t_0, t_1 = t_1,a = a2, m = m2, initial_HR = initial_HR_sim, lambda_wane = lambda_wane)))
          St_pred2[i,j] <- exp(-(-log(St_pred2_cp) +flexsurv:::HweibullPH(t_1, scale = m2, shape = a2)-flexsurv:::HweibullPH(t_0, scale = m2, shape = a2)))
          
        }
      }
    }
    
    plot(km.all, xlim = c(0, max(data_jags$t_pred)))
    # points(x = df_all$time[df_all$arm == 0], 
    #        y = St_mat_quantile[2,df_all$arm == 0], col = "red")
    # points(x = df_all$time[df_all$arm == 1],
    #        y = St_mat_quantile[2,df_all$arm == 1], col = "blue")
    lines(y = colMeans(St_pred), t_pred, col = "blue")
    lines(y = colMeans(St_pred2), t_pred, col = "red")
    #lines(y = St_baseline, time_all, col = "red")
    #lines(y = St_trt, time_all, col = "blue")
    
    rmst_diff_cp <- sfsmisc::integrate.xy(x = t_pred, fx = colMeans(St_pred)) -
      sfsmisc::integrate.xy(x = t_pred, fx = colMeans(St_pred2))
    
    time_all2 <- seq(0, t_max-0.2, by = 0.2)
    km.all_event <- survfit(Surv(time_event, status_true)~arm, data = df_all)
    St_km <- summary(km.all_event, t = time_all2)
    plot(km.all_event)
    lines(y = colMeans(St_pred), t_pred, col = "blue")
    lines(y = colMeans(St_pred2), t_pred, col = "red")
    
    St_baseline_km <- St_km$surv[St_km[["strata"]] == "arm=0"]
    St_trt_km <- St_km$surv[St_km[["strata"]] == "arm=1"]
    
    rmst_diff_true <- sfsmisc::integrate.xy(x = time_all2, fx = St_trt_km)-
      sfsmisc::integrate.xy(x = time_all2 , fx = St_baseline_km) 
    
    # quantile_HR <- apply(HR_mat, 2, quantile, probs = c(0.025, 0.5,0.975))
    # plot(x = 1:length(t_pred), y =quantile_HR[2,] ,xlim = c(0, 50), ylim = c(0,1.5), type = "l",
    #      xlab = "Time (Years)", ylab = "HR")
    # lines(x = 1:length(t_pred), y =quantile_HR[1,],lty= 2)
    # lines(x = 1:length(t_pred), y =quantile_HR[3,], lty = 2)
    
    
    df_res <- data.frame(Model = c("True", "Changepoint"), 
                         Est= c(rmst_diff_true,rmst_diff_cp))
    
    for(p in 1:length(param_mods)){
      
      if(param_mods[p] == "rps"){
        
        mle.mod <- try(flexsurvspline(Surv(time, status)~as.factor(arm),
                                      data = df_all,  k = 1), silent = TRUE)
        
      }else{
        mle.mod <- flexsurvreg(Surv(time, status)~as.factor(arm),
                               data = df_all, dist = param_mods[p],)
      }
      
      if(class(mle.mod) != "try-error"){
        rmst_diff_mod <- summary(mle.mod, t = t_max, type = "rmst")
        df_res <- rbind(df_res,c(param_mods[p], rmst_diff_mod[[1]]$est- rmst_diff_mod[[2]]$est))
        plot(mle.mod, t = seq(0:t_max), xlim = c(0, t_max))
      }else{
        
        df_res <- rbind(df_res,c(param_mods[p], NA))
        
      }
      
    }
    
    df_res$Est <- round(as.numeric(df_res$Est), digits = 2)
    df_res$RMST <- (df_res$Est - rmst_diff_true)^2
    df_res$Abs_perc <-(df_res$Est - rmst_diff_true)/rmst_diff_true
    all_res[[s]] <- df_res
  }
  
  res_all <- matrix(ncol = n_sims, nrow = nrow(all_res[[1]]))
  
  for(i in 1:n_sims){
    res_all[,i] <- all_res[[i]]$Abs_perc
  }
  
  df_mean_vals[[c]] <- data.frame(Model = all_res[[1]]$Model, Avg.Error =  rowMeans(res_all,na.rm = FALSE ))
  
  
}


colnames_all <- c("Changepoint",param_mods)

data.res_1 <- NULL

for(c in 1:nrow(param_vals)){
  data.res_1 <- rbind(data.res_1,as.numeric(t(df_mean_vals[[c]])[2,]))
}
data.res_1 <- data.res_1[,-1]
colnames(data.res_1) <- colnames_all
mod.names.fit <- c(names(sort(abs(colMeans(data.res_1))))[1:3], "Changepoint")
mod.names.fit <- unique(mod.names.fit)
data.res_1 <- cbind(param_vals,data.res_1[,mod.names.fit])

write.csv(data.res_1, paste0(pathway,"data.res_1.csv"))

print(xtable(data.res_1), include.rownames=FALSE)

}
## Treatment delay


df_mean_vals2 <- list()



for(c in 1:nrow(param_vals)){
  
  print(paste0("Scenario ", c))
  
  cp <- 1
  t_0 <- cp
  t_1 = 5
  time_all <- seq(0, t_max, by = 0.001)
  
  n.samp_size <- param_vals$n.samp[c]
  shape_1 <- param_vals$shape_1[c]
  scale_1 <- param_vals$scale_1[c]
  scale_2 <- param_vals$scale_1[c]*param_vals$initial_HR[c]
  
  
  cum_Haz_baseline <- flexsurv:::HweibullPH(time_all,
                                            shape = shape_1,
                                            scale = scale_1 )
  
  time_1 <- seq(0, cp, by = 0.001)
  #Trt Delay
  cum_Haz_trt <- flexsurv:::HweibullPH(time_1,shape = shape_1, scale = scale_1)
  
  time_2 <- seq(cp, t_max, by = 0.001)
  time_2 <- tail(time_2, n = -1)
  
  cum_Haz_trt_cp <- flexsurv:::HweibullPH(time_2,shape = shape_1, scale = scale_2)-flexsurv:::HweibullPH(cp,shape = shape_1, scale = scale_2)
  
  cum_haz_all_trt_cp <- (tail(cum_Haz_trt, n = 1)+cum_Haz_trt_cp)
  
  St_baseline <- exp(-cum_Haz_baseline)
  
  St_trt <- exp(-c(cum_Haz_trt,cum_haz_all_trt_cp))
  plot(time_all, St_baseline, ylim = c(0,1), type = "l")
  lines(time_all,St_trt ,col = "blue")
  
  all_res <- list()
  
  for(s in 1:n_sims){
    
    print(paste0("Simulation ", s))
    
    sim_unif1 <- runif(n.samp)#seq(0,1, length.out = n.samp_size)#runif(n.samp)
    sim_unif1 <- sim_unif1[order(sim_unif1)]
    
    sim_unif2 <- sim_unif1#runif(n.samp)
    sim_unif2 <- sim_unif2[order(sim_unif2)]
    
    time_trt_samp <- time_all[sapply(sim_unif1, function(x){which.min(abs(St_trt-x))})]
    time_comp_samp <- time_all[sapply(sim_unif2, function(x){which.min(abs(St_baseline-x))})]
    
    df_all <- data.frame(time = c(time_trt_samp,time_comp_samp),
                         arm = c(rep(1, length(time_trt_samp)),
                                 rep(0, length(time_comp_samp)))) 
    df_all <- df_all  %>%
      mutate(time_event = time,
             status_true = 1,
             status = ifelse(time > t_cens, 0,1),
             time = ifelse(time > t_cens, t_cens, time))
    
    df_all$time <- df_all$time +0.0001
    km.all <- survfit(Surv(time, status)~arm, data = df_all)
    plot(km.all)
    
    
    data_jags <- list()
    data_jags$N <- nrow(df_all)
    data_jags$time <- df_all$time
    data_jags$status <- df_all$status
    data_jags$MAXX <- quantile(df_all$time[df_all$status == 1], 0.9)
    data_jags$N_CP <- 1
    data_jags$ncovar <- 2
    data_jags$index_pred <- c(min(which(df_all$arm == 1)),min(which(df_all$arm == 0))) #needs to be a vector 
    
    
    data_jags$ncovar_fixed <- 2 # will always be at least 2 even if they are none.. just add a dummy matrix
    data_jags$ncovar_cp <- 2 # will always be at least 2 even if they are none.. just add a dummy matrix
    
    data_jags$t_pred <- seq(0,15, by = 0.1)
    
    #X_mat_trt <- array(NA, c(data_jags$N, 2, data_jags$N_CP +1))
    #X_mat_trt[,,1]<- matrix(c(rep(1, data_jags$N), rep(0, data_jags$N)),  ncol = 2)
    
    #X_mat_trt[,,2]<- matrix(c(rep(1, data_jags$N), TA347_df_all$arm),  ncol = 2)
    X_mat_trt <-  matrix(c(rep(1, data_jags$N), df_all$arm),  ncol = 2)
    
    data_jags$X_mat_trt <-X_mat_trt
    data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N*2)),  ncol = 2)
    #data_jags$X_mat_fixed <-  matrix(c(rep(1, data_jags$N), df_all$arm),  ncol = 2)
    
    data_jags$scenario_common <- c(0,0)
    data_jags$scenario <- 2  #
    data_jags$cp_fix <- c(0, cp, 5) 
    data_jags$cp_fixed <- 0 
    ##Need to fix this
    #Change-inits
    inits <- function(){
      list(unif = runif(data_jags$N_CP,max = data_jags$MAXX/data_jags$N_CP)
           #,
           #beta_cp_0 = c(0),
           #beta_cp = rnorm(data_jags$N_CP+1)
      )
    }
    
    # Tweak the model

    mod.cp <- runjags::run.jags(
      model = get("jags.piecewise.wei_chng"),
      #model = get(paste0("jags.piecewise.wei_chng_common_haz")),
      data = data_jags,
      n.chains = 2,
      monitor = c("cp", "beta_cp","beta_cp_anc",
                  "prec","sd", "cum_haz_pred", "cum_haz", "total_loglik","beta_covar_fixed", "loglik"), 
      sample=n.samp, 
      thin = n.thin, 
      burnin = n.burnin,
      inits = inits, 
      method ='rjparallel')
    
    #add.summary(mod.cp, "cp")
    add.summary(mod.cp, "total_loglik")
    
    param_names <- colnames(mod.cp$mcmc[[1]])
    #plot(density(mcmc_output[,"total_loglik"]))
    mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])
    
    plot(density(mcmc_output[,grep("cp[2]",param_names, fixed = T)]))
    
    St_mat <- exp(-mcmc_output[,grep("cum_haz_pred[",param_names, fixed = T)])
    #St_mat <- exp(-mcmc_output[,grep("cum_haz[",param_names, fixed = T)])
    St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))
 
    St_mat_trt <- St_mat[,1:(ncol(St_mat)/2)]
    St_mat_comp <- St_mat[,(ncol(St_mat)/2 +1 ):ncol(St_mat)]
    
    plot(km.all, xlim = c(0, max(data_jags$t_pred)), ylab = "Survival", xlab = "Time (Months)")
    lines(x =data_jags$t_pred , y =St_mat_quantile[2,1:(ncol(St_mat_quantile)/2)], col = "blue")
    lines(x =data_jags$t_pred , y =St_mat_quantile[2,(ncol(St_mat_quantile)/2 +1 ):ncol(St_mat_quantile)], col = "red")
    
    
    #Proof that it is correct
    
    # maxLL <-mcmc_output[which.max(mcmc_output[,"total_loglik" ]), grep("cp", colnames(mcmc_output))]
    # 
    # split_data <- survSplit(Surv(time, status)~factor(arm),
    #           data = df_all, cut = maxLL["cp[2]"], episode ="interval")
    # 
    # LL_vec <- rep(NA, nrow(split_data))
    # 
    # sum(LL_vec)
    # 
    # scale_1_max <- exp(maxLL["beta_cp[1,1]"])
    # shape_1_max <- exp(maxLL["beta_cp_anc[1,1]"])
    # 
    # shape_2_max <- exp(maxLL["beta_cp_anc[2,1]"])
    # scale_2_1_max <- exp(maxLL["beta_cp[2,1]"])
    # scale_2_2_max <- exp(maxLL["beta_cp[2,1]"] +maxLL["beta_cp[2,2]"])
    # 
    # for(i in 1:nrow(split_data)){
    #   
    #   if(split_data$interval[i] ==1){
    #     LL_vec[i] <- log(flexsurv::hweibullPH(split_data$time[i],scale = scale_1_max,shape = shape_1_max))*split_data$status[i] - flexsurv::HweibullPH(split_data$time[i],scale = scale_1_max,shape = shape_1_max)
    #   }
    #   
    #   if(split_data$interval[i] ==2){
    #     
    #     if(split_data["factor(arm)"][i,1] == 0){
    #       LL_vec[i] <- log(flexsurv::hweibullPH(split_data$time[i],scale = scale_2_1_max,shape = shape_2_max))*split_data$status[i] - 
    #         flexsurv::HweibullPH(split_data$time[i],scale = scale_2_1_max,shape = shape_2_max) +
    #         flexsurv::HweibullPH(maxLL["cp[2]"],scale = scale_2_1_max,shape = shape_2_max)
    #     }else{
    #       LL_vec[i] <- log(flexsurv::hweibullPH(split_data$time[i],scale = scale_2_2_max,shape = shape_2_max))*split_data$status[i] - 
    #          flexsurv::HweibullPH(split_data$time[i],scale = scale_2_2_max,shape = shape_2_max) +
    #         flexsurv::HweibullPH(maxLL["cp[2]"],scale = scale_2_2_max,shape = shape_2_max)
    #       
    #     }
    # 
    #   }
    #   
    #   
    # }
    
    #exp(-1.2012078)
    
    # add.summary(mod.cp, "cp")
    # add.summary(mod.cp, "beta_covar_fixed")
    # add.summary(mod.cp, "total_loglik")
    
    # loglik_mat <- mcmc_output[,grep("loglik",param_names)]
    # loglik_mat <- loglik_mat[,-grep("total_loglik",colnames(loglik_mat))]
    # 
    # 
    # mean_cp <- mean(mcmc_output[,grep("cp[2]",param_names,fixed = T)])
    # waic <- loo::waic(loglik_mat)[["estimates"]][3, 1]
    # pml <- -2 * sum(log(nrow(loglik_mat)/colSums(1/exp(loglik_mat))))
    
    
    rmst_diff_cp <- sfsmisc::integrate.xy(x = data_jags$t_pred, fx = colMeans(St_mat_trt)) -
      sfsmisc::integrate.xy(x = data_jags$t_pred, fx = colMeans(St_mat_comp))
    
    time_all2 <- seq(0, max(data_jags$t_pred), by = 0.2)
    km.all_event <- survfit(Surv(time_event, status_true)~arm, data = df_all)
    St_km <- summary(km.all_event, t = time_all2)
    plot(km.all_event)
    lines(y = colMeans(St_mat_trt), data_jags$t_pred, col = "blue")
    lines(y = colMeans(St_mat_comp), data_jags$t_pred, col = "red")
    
    St_baseline_km <- St_km$surv[St_km[["strata"]] == "arm=0"]
    St_trt_km <- St_km$surv[St_km[["strata"]] == "arm=1"]
    time_baseline <-  St_km$time[St_km[["strata"]] == "arm=0"]
    time_trt <-  St_km$time[St_km[["strata"]] == "arm=1"]
    
    rmst_diff_true <- sfsmisc::integrate.xy(x = time_trt, fx = St_trt_km)-
      sfsmisc::integrate.xy(x = time_baseline , fx = St_baseline_km) 
    
    # rmst_diff_true <- sfsmisc::integrate.xy(x = time_all, fx = St_trt)-
    #   sfsmisc::integrate.xy(x = time_all , fx = St_baseline) 
     
    param_mods <- c("exp","weibullPH",  "lnorm", "gamma", 
                    "gompertz", "llogis","gengamma", "rps")
    
    df_res <- data.frame(Model = c("True", "Changepoint"), Est= c(rmst_diff_true,rmst_diff_cp ))
    
    for(p in 1:length(param_mods)){
      
      if(param_mods[p] == "rps"){
        
        mle.mod <- try(flexsurvspline(Surv(time, status)~as.factor(arm),
                                      data = df_all,  k = 1), silent = TRUE)
        
      }else{
        mle.mod <- flexsurvreg(Surv(time, status)~as.factor(arm),
                               data = df_all, dist = param_mods[p],)
      }
      
      if(class(mle.mod) != "try-error"){
        rmst_diff_mod <- summary(mle.mod, t = max(data_jags$t_pred), type = "rmst")
        df_res <- rbind(df_res,c(param_mods[p], rmst_diff_mod[[1]]$est- rmst_diff_mod[[2]]$est))
        plot(mle.mod, t = seq(0:t_max), xlim = c(0, t_max))
      }else{
        
        df_res <- rbind(df_res,c(param_mods[p], NA))
        
      }
      
    }
    
    df_res$Est <- round(as.numeric(df_res$Est), digits = 2)
    df_res$RMST <- (df_res$Est - rmst_diff_true)^2
    df_res$Abs_perc <- abs(df_res$Est - rmst_diff_true)/rmst_diff_true
    all_res[[s]] <- df_res
  }
  
  res_all <- matrix(ncol = n_sims, nrow = nrow(all_res[[1]]))
  
  for(i in 1:n_sims){
    res_all[,i] <- all_res[[i]]$Abs_perc
  }
  
  df_mean_vals2[[c]] <- data.frame(Model = all_res[[1]]$Model, Avg.Error =  
                                     rowMeans(res_all,na.rm = FALSE ))
  
}

colnames_all <- c("Changepoint",param_mods)

data.res_2 <- NULL

for(c in 1:nrow(param_vals)){
  data.res_2 <- rbind(data.res_2,as.numeric(t(df_mean_vals2[[c]])[2,]))
}
data.res_2 <- data.res_2[,-1]
colnames(data.res_2) <- colnames_all
mod.names.fit <- names(sort(abs(colMeans(data.res_2,na.rm = T))))
mod.names.fit <- unique(mod.names.fit)
data.res_2 <- cbind(param_vals,data.res_2)

write.csv(data.res_2, paste0(pathway,Sys.Date(),"data.res_2.csv"))

#print(xtable(read.csv(paste0(pathway,"2023-07-29data.res_2.csv"))[,2:8]), include.rownames=FALSE)
print(xtable(read.csv(paste0(pathway,Sys.Date(),"data.res_2.csv"))[,2:8]), include.rownames=FALSE)

### Change-point after cp
#Confirm that the calculations are correct....

df_mean_vals3 <- list()

for(c in 1:nrow(param_vals)){
  
  print(paste0("Scenario ", c))
 all_res <- list()
  cp <- 2
  t_0 <- cp
  time_1 <- seq(0, cp, by = 0.001)
  time_all <- seq(0, t_max, by = 0.001)
  t_1 = 5
  n.samp_size <- param_vals$n.samp[c]
  shape_1 <- param_vals$shape_1[c]
  scale_1 <- param_vals$scale_1[c]
  scale_2 <- param_vals$scale_1[c]*param_vals$initial_HR[c]
  

  cum_Haz_baseline <- flexsurv:::HweibullPH(time_all,
                                          shape = shape_1,
                                          scale = scale_1 )
  cum_Haz_trt <- flexsurv:::HweibullPH(time_1,
                                       shape = shape_1,
                                       scale = scale_2)

  time_2 <- seq(cp, t_max, by = 0.001)
  time_2 <- tail(time_2, n = -1)

  cum_Haz_trt_cp <- flexsurv:::HweibullPH(time_2,shape = shape_1, 
                                          scale = scale_1)-flexsurv:::HweibullPH(cp,shape = shape_1,
                                                                                 scale = scale_1)

  cum_haz_all_trt_cp <- (tail(cum_Haz_trt, n = 1)+cum_Haz_trt_cp)

  St_baseline <- exp(-cum_Haz_baseline)
  St_trt <- exp(-c(cum_Haz_trt,cum_haz_all_trt_cp))
 plot(time_all, St_baseline, ylim = c(0,1), type = "l")
 lines(time_all,St_trt ,col = "blue")

for(s in 1:n_sims){
  
  print(paste0("Simulation ", s))
  
sim_unif1 <- runif(n.samp_size)#seq(0,1,length.out = n.samp_size) #runif(n.samp_size)
sim_unif1 <- sim_unif1[order(sim_unif1)]

sim_unif2 <- sim_unif1#runif(n.samp)
sim_unif2 <- sim_unif2[order(sim_unif2)]


time_trt_samp <- time_all[sapply(sim_unif1, function(x){which.min(abs(St_trt-x))})]
time_comp_samp <- time_all[sapply(sim_unif2, function(x){which.min(abs(St_baseline-x))})]


df_all <- data.frame(time = c(time_trt_samp,time_comp_samp),
                     arm = c(rep(1, length(time_trt_samp)),
                             rep(0, length(time_comp_samp)))) 

df_all$time <- df_all$time +0.0000001

df_all <- df_all  %>%
  mutate(time_event = time,
         status_true = 1,
         status = ifelse(time > t_cens, 0,1),
         time = ifelse(time > t_cens, t_cens, time))

km.all <- survfit(Surv(time, status)~arm, data = df_all)
plot(km.all)

#Change-point
# df_all2 <- survSplit(Surv(time, status) ~arm, df_all,
#                   cut=1, episode ="timegroup")
# fit2 <- coxph(Surv(tstart, time, status) ~ arm* strata(timegroup), data= df_all2)
# 
# c(time1= exp(coef(fit2)[1]),
#   time2  = exp(sum(coef(fit2)[c(1,2)])))
# 
# km.all2 <- survfit(Surv(time_event, status_true)~arm, data = df_all)
# plot(km.all2)


#print(km.all, print.rmean=TRUE)



data_jags <- list()
data_jags$N <- nrow(df_all)
data_jags$time <- df_all$time
data_jags$status <- df_all$status
data_jags$MAXX <- max(df_all$time)
data_jags$N_CP <- 1
data_jags$ncovar <- 2
data_jags$index_pred <- c(min(which(df_all$arm == 1)),min(which(df_all$arm == 0))) #needs to be a vector 


data_jags$N_CP <-1
data_jags$ncovar_fixed <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix
data_jags$ncovar_cp <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix

data_jags$t_pred <- seq(0.01,t_max - 0.2, by = 0.2)

#X_mat_fixed is just a dummy matrix
#PH Model with Changepoint
data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N), rep(0, data_jags$N)),  ncol = data_jags$ncovar_fixed)
#independent parameters 
data_jags$X_mat_trt <- matrix(c(rep(1, data_jags$N), df_all$arm),  ncol = 2)

data_jags$scenario_common <- c(0,0)
data_jags$scenario <- 2 #
data_jags$cp_fix <- c(0, cp, 5) 
data_jags$cp_fixed <-0

#Change-inits
inits <- function(){
  list(unif = cp,
       beta_cp_x = matrix(c(log(scale_1), log(param_vals$initial_HR[c]),
                          log(scale_1),0), nrow = 2, byrow =T),
       beta_cp_anc_x = matrix(c(log(shape_1), 0,
                              log(shape_1),0), nrow = 2, byrow =T)
       #runif(data_jags$N_CP,max = data_jags$MAXX/data_jags$N_CP)
       #,
       #beta_cp_0 = c(0),
       #beta_cp = rnorm(data_jags$N_CP+1)
  )
}



#n.samp <- 5000
#n.burnin <- 1000
mod.cp <- runjags::run.jags(
  model = get(paste0("jags.piecewise.wei_chng")),
  #model = jags.piecewise.expo_chng,
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp", "beta_cp_anc",
              "prec","sd", "cum_haz_pred", "cum_haz", "total_loglik","loglik","beta_covar_fixed"), 
  sample=n.samp, 
  thin = n.thin, 
  burnin = n.burnin,
  inits = inits, 
  method ='rjparallel')

#add.summary(mod.cp, "cp")
#add.summary(mod.cp, "total_loglik")

param_names <- colnames(mod.cp$mcmc[[1]])
mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])

plot(density(mcmc_output[,grep("total_loglik",param_names)]))
plot(density(mcmc_output[,grep("cp[2]",param_names, fixed = T)]))



St_mat <- exp(-mcmc_output[,grep("cum_haz_pred",param_names)])
St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))

St_mat_trt <- St_mat[,1:(ncol(St_mat)/2)]
St_mat_comp <- St_mat[,(ncol(St_mat)/2 +1 ):ncol(St_mat)]


summary_output <- mcmc_output[,grep("beta_cp|cp|beta_covar_fixed",param_names)]
summary_output <- summary_output[, c(2,6,7) ]

summary_output[,2:3] <- exp(summary_output[,2:3])
colMeans(summary_output)

#summary_output[, c(2,3, 10)] <- exp(summary_output[, c(2,3, 10)])


plot(km.all, xlim = c(0, max(data_jags$t_pred)), xlab = "Time (Years)", ylab = "S(t)")
#plot(km.all_event,  xlab = "Time (Years)", ylab = "S(t)")
lines(x =data_jags$t_pred , y =colMeans(St_mat_comp), col = "blue")
lines(x =data_jags$t_pred , y =colMeans(St_mat_trt), col = "red")
lines(x = time_all , y = St_baseline ,col = "red", lty = 2)
lines(x = time_all , y = St_trt, col = "blue", lty = 2)

rmst_diff_cp <- sfsmisc::integrate.xy(x = data_jags$t_pred, fx = colMeans(St_mat_trt)) -
  sfsmisc::integrate.xy(x = data_jags$t_pred, fx = colMeans(St_mat_comp))

time_all2 <- seq(0, max(data_jags$t_pred), by = 0.2)
km.all_event <- survfit(Surv(time_event, status_true)~arm, data = df_all)
St_km <- summary(km.all_event, t = time_all2)
plot(km.all_event)
lines(y = colMeans(St_mat_trt), data_jags$t_pred, col = "blue")
lines(y = colMeans(St_mat_comp), data_jags$t_pred, col = "red")

St_baseline_km <- St_km$surv[St_km[["strata"]] == "arm=0"]
St_trt_km <- St_km$surv[St_km[["strata"]] == "arm=1"]
time_baseline <-  St_km$time[St_km[["strata"]] == "arm=0"]
time_trt <-  St_km$time[St_km[["strata"]] == "arm=1"]

rmst_diff_true <- sfsmisc::integrate.xy(x = time_trt, fx = St_trt_km)-
  sfsmisc::integrate.xy(x = time_baseline , fx = St_baseline_km) 

# rmst_diff_true <- sfsmisc::integrate.xy(x = time_all, fx = St_trt)-
#   sfsmisc::integrate.xy(x = time_all , fx = St_baseline) 

df_res <- data.frame(Model = c("True", "Changepoint"), Est= c(rmst_diff_true,rmst_diff_cp ))


for(p in 1:length(param_mods)){
  
  if(param_mods[p] == "rps"){
    
    mle.mod <- try(flexsurvspline(Surv(time, status)~as.factor(arm),
                              data = df_all,  k = 1), silent = TRUE)
    
  }else{
    mle.mod <- flexsurvreg(Surv(time, status)~as.factor(arm),
                           data = df_all, dist = param_mods[p],)
  }
  
  if(class(mle.mod) != "try-error"){
    rmst_diff_mod <- summary(mle.mod, t = max(data_jags$t_pred), type = "rmst")
    df_res <- rbind(df_res,c(param_mods[p], rmst_diff_mod[[1]]$est- rmst_diff_mod[[2]]$est))
    plot(mle.mod, t = seq(0:t_max), xlim = c(0, t_max))
  }else{
      
    df_res <- rbind(df_res,c(param_mods[p], NA))
    
  }

}



df_res$Est <- as.numeric(df_res$Est)
df_res$RMST <- (df_res$Est - rmst_diff_true)^2
df_res$Abs_perc <- abs(df_res$Est - rmst_diff_true)/rmst_diff_true
all_res[[s]] <- df_res

}


res_all <- matrix(ncol = n_sims, nrow = nrow(all_res[[1]]))

  for(i in 1:n_sims){
    res_all[,i] <- all_res[[i]]$Abs_perc
  }

df_mean_vals3[[c]] <- data.frame(Model = all_res[[1]]$Model,
                                 Avg.Error =  rowMeans(res_all,na.rm = FALSE ))

}


colnames_all <- c("Changepoint",param_mods)

data.res_3 <- NULL

for(c in 1:nrow(param_vals)){
  data.res_3 <- rbind(data.res_3,as.numeric(t(df_mean_vals3[[c]])[2,]))
}
data.res_3 <- data.res_3[,-1]
colnames(data.res_3) <- colnames_all
mod.names.fit <- names(sort(abs(colMeans(data.res_3,na.rm = T))))
mod.names.fit <- unique(mod.names.fit)
data.res_3 <- cbind(param_vals,data.res_3[,mod.names.fit])

write.csv(data.res_3, paste0(pathway,Sys.Date(),"data.res_3.csv"))

#need to understand why this is not maximizing the likelihood - it is just because of lack of identifiability

#maybe get rid of zero prior for the data.frame
print(xtable(read.csv(paste0(pathway,"2023-08-05data.res_3.csv"))[,c(2:5, 7:9)]), include.rownames=FALSE)

print(xtable(read.csv(paste0(pathway,Sys.Date(),"data.res_3.csv"))[,c(2:6, 8:9)]), include.rownames=FALSE)


