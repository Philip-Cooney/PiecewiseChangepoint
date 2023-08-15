pathway <- "~/Change-Point Simulation Studies/Simulation Study 2023/"
library("rjags")
#library("survival")
library("flexsurv")
library("survminer")
library("runjags", lib.loc = "~/R-packages")
library("Epi", lib.loc = "~/R-packages/")
library("bshazard", lib.loc = "~/R-packages/")
library("dplyr")
library("xlsx")

source(paste0(pathway, "Jags_codes.R"))


E1690.dat <- read.table(paste0(pathway,"e1690.missing.dat"), 
                        header=TRUE)

# Ibhrahim data
#Drop PFS events with time equal zero
E1690.dat <- E1690.dat[-which(E1690.dat$FAILTIME ==0),]

#Convert to the correct notation for survival objects
E1690.dat[which(E1690.dat$FAILCENS == 1),"FAILCENS"] <-0
E1690.dat[which(E1690.dat$FAILCENS == 2),"FAILCENS"] <-1
E1690.dat[which(E1690.dat$SURVCENS == 1),"SURVCENS"] <-0
E1690.dat[which(E1690.dat$SURVCENS == 2),"SURVCENS"] <-1
E1690.dat$SEX <- as.numeric(E1690.dat$SEX)
E1690.dat$AGE <- as.numeric(E1690.dat$AGE)
E1690.dat$AGE[is.na(E1690.dat$AGE)] <- mean(E1690.dat$AGE, na.rm = T)
E1690.dat$AGE_scale <- as.numeric(scale(E1690.dat$AGE))
E1690.dat$AGE_scale_abs <- abs(E1690.dat$AGE_scale)
# fit.OS <- survfit(Surv(SURVTIME, SURVCENS)~TRT, 
#                   data = E1690.dat[which(E1690.dat$STUDY == "1684"),])

n.samp <- 5000
n.burnin <- 1000
n.thin <- 1
fit.OS <- survfit(Surv(SURVTIME, SURVCENS)~TRT, 
                  data = E1690.dat)

cox.fit.OS <- coxph(Surv(SURVTIME, SURVCENS)~as.factor(TRT) + scale(AGE), 
                    data = E1690.dat)

cox.zph(cox.fit.OS)

fit.mle <- flexsurvreg(Surv(SURVTIME, SURVCENS)~as.factor(TRT), 
                       data = E1690.dat, dist = "weibullPH")

fs3 <- flexsurvspline(Surv(SURVTIME, SURVCENS) ~ factor(TRT) , data = E1690.dat,
                      k = 1, anc = list(gamma1 = ~ as.factor(TRT)))
plot(fs3, t = c(0:50), xlim = c(0,50))
fit_trt1<-bshazard(Surv(SURVTIME, SURVCENS==1) ~ 1,data=E1690.dat %>% filter(TRT == 1))
fit_trt2<-bshazard(Surv(SURVTIME, SURVCENS==1) ~ 1,data=E1690.dat %>% filter(TRT == 2))

png(paste0(pathway, "Hazards-E1690.png"), width = 10, height = 5, units = 'in', res = 300)
plot(fit_trt1$time, fit_trt1$hazard,xlab='Time (Years)', ylab='Rate/person-years',type='l', ylim = c(0, 0.3), col = "blue",
     main = "Hazards for E1684 and E1690 trials")
#lines(fit_trt1$time, fit_trt1$lower.ci, lty=2, lwd=1, col = "blue")
#lines(fit_trt1$time, fit_trt1$upper.ci, lty=2, lwd=1, col = "blue")

#points(fit$raw.data$time,fit$raw.data$raw.hazard,cex=.3, lwd=3,col=1)
lines(fit_trt2$time, fit_trt2$hazard, lty=1, lwd=1, col = "red")
#lines(fit_trt2$time, fit_trt2$lower.ci, lty=2, lwd=1, col = "red")
#lines(fit_trt2$time, fit_trt2$upper.ci, lty=2, lwd=1, col = "red")

dev.off()


data_jags <- list()
data_jags$N <- nrow(E1690.dat)
data_jags$time <- E1690.dat$SURVTIME
data_jags$status <- E1690.dat$SURVCENS
data_jags$MAXX <- max(E1690.dat$SURVTIME)
data_jags$N_CP <- 1
data_jags$ncovar <- 2
data_jags$index_pred <- c(which(E1690.dat$AGE_scale_abs== min(E1690.dat$AGE_scale_abs[E1690.dat$TRT == 1])),
                          which(E1690.dat$AGE_scale_abs== min(E1690.dat$AGE_scale_abs[E1690.dat$TRT == 2]))) #needs to be a vector 

data_jags$N_CP <-1
data_jags$ncovar_fixed <- 3 # will always be at least 2 even if thery are none.. just add a dummy matrix
data_jags$ncovar_cp <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix

data_jags$t_pred <- seq(0,20, by = 0.2) #Some patients will be over 100

#X_mat_fixed is just a dummy matrix
#PH Model with Changepoint
data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N), E1690.dat$TRT-1, E1690.dat$AGE_scale),  ncol = data_jags$ncovar_fixed)
#independent parameters 
data_jags$X_mat_trt <- matrix(c(rep(1, data_jags$N),rep(0, data_jags$N)),  ncol = 2)

#Change-inits
inits <- function(){
  list(unif = runif(data_jags$N_CP,max = data_jags$MAXX/data_jags$N_CP)
       #,
       #beta_cp_0 = c(0),
       #beta_cp = rnorm(data_jags$N_CP+1)
  )
}

#Scenario 2 -- Common Hazards for final interval
data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N), E1690.dat$TRT-1, E1690.dat$AGE_scale),  ncol = data_jags$ncovar_fixed)


data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N*2),E1690.dat$AGE_scale),  ncol = data_jags$ncovar_fixed)
X_mat_trt <- array(NA, c(data_jags$N, 2, 2))
X_mat_trt[,,1]<- matrix(c(rep(1, data_jags$N), E1690.dat$TRT-1),  ncol = 2)
X_mat_trt[,,2]<- matrix(c(rep(1, data_jags$N), rep(0, data_jags$N)),  ncol = 2)

data_jags$X_mat_trt <-X_mat_trt

##Need to fix this




#Add in GPM 
Conditional_Death_df <- read.xlsx(paste0(pathway, "Conditional_Death_UK.xlsx"),1)
colnames(Conditional_Death_df) <- c("age", "Male_cond_death", "Female_cond_death")
Conditional_Death_df <- rbind(Conditional_Death_df,c(101, .99999, .99999)) #All patients die at 101
mean_age_baseline_example <- age_baseline_example <-mean(E1690.dat$AGE)
prop_male <- 1-mean(as.numeric(E1690.dat$SEX)-1, na.rm = T)
max_age <- 90
time_factor <- 1

time_horizon <- 100
df_temp <- Conditional_Death_df
df_temp$Male_cond_haz <- -log(1-df_temp$Male_cond_death)/time_factor 
df_temp$Female_cond_haz <- -log(1-df_temp$Female_cond_death)/time_factor

n_row_df_temp <- nrow(df_temp)
time_partial_vec <- rep((0:(time_factor-1))/time_factor)

df_temp <- do.call("rbind", replicate(time_factor, df_temp, simplify = FALSE)) %>%
  arrange(age)
df_temp$age <- df_temp$age+  rep(time_partial_vec,times = n_row_df_temp)

#Cohort approach
df_temp_cohort <- df_temp %>% filter(age >= mean_age_baseline_example & age <= time_horizon)
#Cohort - Static Gender proportion
df_temp_cohort[, "mix_haz_static"] <- df_temp_cohort[,"Male_cond_haz"]*prop_male + df_temp_cohort[,"Female_cond_haz"]*(1-prop_male)




df <- E1690.dat

df <- df %>% mutate(sex = SEX-1,
                    age = AGE)

df$SEX[is.na(df$SEX)] <- 0
GMP_haz_jags <- GMP_cum_haz_jags <- rep(NA,nrow(df))
for(i in 1:nrow(df)){
  age_start <- which.min(abs(df$age[i] - df_temp$age))
  age_selc <-  which.min(abs((df$age[i]+df$SURVTIME[i]/time_factor) - df_temp$age))
  
  if(df$SEX[i] == 0){
    GMP_haz_jags[i] <- df_temp[age_selc, "Male_cond_haz"] 
    GMP_cum_haz_jags[i] <- sum(df_temp[age_start:age_selc,"Male_cond_haz"])
    
  }else{
    GMP_haz_jags[i] <- df_temp[age_selc, "Female_cond_haz"] #assume rounding cancels ##aproximation
    GMP_cum_haz_jags[i] <- sum(df_temp[age_start:age_selc,"Female_cond_haz"])
  } 
}


data_jags$GMP_haz <-GMP_haz_jags
data_jags$GMP_cum_haz <-GMP_cum_haz_jags




jags.piecewise.wei_chng_common_haz_gmp <-"

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
  }
  cp[2:(N_CP+1)] <- sort(unif)
  
  for(i in 2:(N_CP+1)){
    diff[i-1] <- cp[i]-cp[i-1]
  }
  
  log_prior <- loggam(2*N_CP +2) + sum(log(diff)) - ((2*N_CP +1)*log(MAXX))
  zero_prior ~ dpois(C - log_prior)
  
  #Prior for the model parameters
  
  for(i in 1:2){
    sd[i]  <- 0.5 # ~ dunif(0.2,5)
    prec[i] <- pow(sd[i], -2)
  }
  
  for(k in 1:N_CP){ #Independent up to the last interval
    for(j in 1:2){ #Always the number of basic parameters plus 1 for a covariate 
      beta_cp[k,j] ~ dnorm(0,prec[1])
      beta_cp_anc[k,j] ~ dnorm(0,prec[1])
    }
   }
   for(j in 1:2){ #Always the number of basic parameters plus 1 for a covariate 
      beta_cp_common[j] ~ dnorm(0,prec[1])
      beta_cp_anc_common[j] ~ dnorm(0,prec[1])
    }
	
   for(k in (N_CP+1):(N_CP+1)){ #Independent up to the last interval
    for(j in 1:2){ #Always the number of basic parameters plus 1 for a covariate 
      beta_cp[k,j] <- beta_cp_common[j]
      beta_cp_anc[k,j]<- beta_cp_anc_common[j]
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
		   ########### CHANGE - X_mat_trt is now an array
        linpred_cp_param[i,k,p] <- X_mat_trt[i,p,k]*beta_cp[k,p]
      }
      linpred_cp[i,k] <- sum(linpred_cp_param[i,k,])
      param_1[i,k] <- exp(linpred_cp[i,k] + linpred_fixed[i]) #lambda
      param_2[i,k] <- exp(beta_cp_anc[k,1]) #We allow an independent model for the first time-point 
      
      log_haz_seg[i,k] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))*X_ind[i,k] 
      cum_haz_seg[i,k] <- param_1[i,k]*pow(X[i,k]+cp[k],param_2[i,k]) -  param_1[i,k]*pow(cp[k],param_2[i,k])
      
    }
    
    log_haz[i] <- sum(log_haz_seg[i,])
    haz[i] <- exp(log_haz[i]) + GMP_haz[i]
    log_haz_final[i] <- log(haz[i])
    cum_haz[i] <- sum(cum_haz_seg[i,])
    cum_haz_final[i] <- cum_haz[i] + GMP_cum_haz[i]    
    
    
    loglik[i] = status[i]*log_haz_final[i] - cum_haz_final[i]
    zero[i] ~ dpois(C - loglik[i])
    
  }
  

    #Predict Survival
  for(i in 1:length(t_pred)){
    
    for(k in 1:(N_CP+1)){
      
      
      X_pred[i,k] = max(min(t_pred[i], cp[k+1]) - cp[k],0) 
      X_ind_pred[i,k] = step(t_pred[i]-cp[k])*step(cp[k+1]-t_pred[i])
      
      for(q in 1:length(index_pred)){
        
        for(p in 1:ncovar_cp){
			########### CHANGE - X_mat_trt is now an array
          linpred_cp_param_pred[i,k,q,p] <- X_mat_trt[index_pred[q],p,k]*beta_cp[k,p]
        }
        linpred_cp_pred[i,k,q] <- sum(linpred_cp_param_pred[i,k,q,])
        
        param_1_pred[i,k,q] <- exp(linpred_cp_pred[i,k,q]  + linpred_fixed[index_pred[q]]) #lambda
        param_2_pred[i,k,q] <- exp(beta_cp_anc[k,1]) #Can be made constant across all if k = 1
        #param_2_pred[i,k,q] <- exp(beta_cp_anc[k,1] + equals(k,1)*beta_cp_anc[k,2]*X_mat_trt[index_pred[q],2]) #We allow an independent model for the first time-point 
      
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



mod.cp <- runjags::run.jags(
  #model = jags.piecewise.expo_chng_common_haz,
  model = get(paste0("jags.piecewise.wei_chng_common_haz_gmp")),
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp", "prec","sd", "cum_haz_pred", "cum_haz", "total_loglik","beta_covar_fixed", "loglik"), 
  sample=n.samp, 
  thin = n.thin, 
  burnin = n.burnin,
  inits = inits, 
  method ='rjparallel')


#add.summary(mod.cp)

mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])
param_names <- colnames(mod.cp$mcmc[[1]]) 

GMP_cum_haz_mat <- matrix(nrow = nrow(df), ncol = length(data_jags$t_pred))
for(i in 1:nrow(df)){
  age_start <- which.min(abs(df$age[i] - df_temp$age))
  
  if(df$SEX[i] == 0){
    GMP_cum_haz_mat[i,] <- cumsum(df_temp[age_start+data_jags$t_pred,"Male_cond_haz"])
  }else{
    GMP_cum_haz_mat[i,] <- cumsum(df_temp[age_start+data_jags$t_pred,"Female_cond_haz"])
  } 
  
  GMP_cum_haz_mat[,1] <- 0 #assume zero cumulative hazard at age at baseline
}

cum_haz_pred <- mcmc_output[,grep("cum_haz_pred[",param_names, fixed = T)]
cum_haz_pred_comp <- cum_haz_pred[,1:(ncol(cum_haz_pred)/2)]
cum_haz_pred_trt <-cum_haz_pred[,(ncol(cum_haz_pred)/2 +1 ):ncol(cum_haz_pred)]

#cum_haz_pred <- cum_haz_pred[,1:length(data_jags$t_pred)]
cum_haz_array_comp <- array(dim = c(nrow(cum_haz_pred_comp), ncol(cum_haz_pred_comp),sum(df$TRT == 1)))
cum_haz_array_trt <- array(dim = c(nrow(cum_haz_pred_trt), ncol(cum_haz_pred_trt),sum(df$TRT == 2)))
index_trt <- index_comp <- 1

for(i in 1:data_jags$N){
  
  
  if(df$TRT[i] == 1){
    GMP_curr <- matrix(GMP_cum_haz_mat[i,], ncol = length(GMP_cum_haz_mat[i,]), nrow = nrow(cum_haz_array_comp), byrow = T)
    
    cum_haz_array_comp[,,index_comp] <- cum_haz_pred_comp +  GMP_curr
    index_comp <- index_comp +1
  }
  if(df$TRT[i] == 2){
    GMP_curr <- matrix(GMP_cum_haz_mat[i,], ncol = length(GMP_cum_haz_mat[i,]), nrow = nrow(cum_haz_array_trt), byrow = T)
    
    cum_haz_array_trt[,,index_trt] <- cum_haz_pred_trt +  GMP_curr
    index_trt <- index_trt +1
  }
  # if(any(is.na(cum_haz_array[,,i]))){
  #   print(paste0("NA in ", i))
  # }
  
}
Surv_array_trt <- exp(-cum_haz_array_trt)
Surv_array_comp <- exp(-cum_haz_array_comp)


St_final_internal_trt <- colMeans(apply(Surv_array_trt,c(1,2),mean))
St_final_internal_comp <- colMeans(apply(Surv_array_comp,c(1,2),mean))

png(paste0(pathway, "Common-Final-Weibull-model-GPM.png"), width = 10, height = 5, units = 'in', res = 300)
plot(fit.OS,
     #xlim = c(0, max(data_jags$t_pred)),
     xlab = "Time (Years)", ylab = "S(t)", 
     main = "Weibull 1 change-point model (Equal final hazards) - General Population Mortality", xlim = c(0,15))
lines(x =data_jags$t_pred , y =St_final_internal_trt, col = "blue")
lines(x =data_jags$t_pred , y =St_final_internal_comp, col = "red")
dev.off()

s <- summary(fit.OS, times = fit.OS$time, extend = T)

plot_data <- tibble(
  'time' = s$time,
  'n.risk' = s$n.risk,
  'n.event' = s$n.event,
  'n.censor' = s$n.censor,
  'Survival' = s$surv,
  'std.error' = s$std.err,
  'Treatment' = s$strata
)
plot_data <- plot_data %>% mutate(Treatment = ifelse(Treatment == "TRT=2", "INF", "OBS")) 

Surv_extrap <- data.frame(time = rep(data_jags$t_pred,2),
                          Survival =c(St_final_internal_trt,St_final_internal_comp), 
                          Treatment= c(rep("INF",length(data_jags$t_pred)),
                                       rep("OBS",length(data_jags$t_pred))))
mean.cp <- mean(mcmc_output[,"cp[2]"])
dev <-abs(Surv_extrap$time - mean.cp)

p <- ggplot(plot_data, aes(x = time, y = Survival, color = Treatment))+
  geom_step( size = 1)+
  ylim(c(0,1))+
  xlim(c(0,15))+
  xlab("Time (Years)")+
  ggtitle("Weibull 1 change-point model (Proportional Hazards) - General Population Mortality")+
  geom_line(data = Surv_extrap,
            aes(x = time, y = Survival,color = Treatment ))+
  geom_point(data = data.frame(time = mean.cp,
                               Survival = mean(pull(Surv_extrap,Survival)[which(dev == min(dev))])),
             aes(x = time, Survival), shape = 23, fill = "green",
             color = "darkred", size = 4, inherit.aes = F, stroke  = 1.5)+
  theme_bw()
ggsave(paste0(pathway, "Common-Final-Weibull-model-GPM.png"), width = 10, height = 5, units = 'in')










param_names <- colnames(mod.cp$mcmc[[1]]) #rownames(expo.mod_1chng$BUGSoutput$summary)
mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])
St_mat <- exp(-mcmc_output[,grep("cum_haz_pred",param_names)])

St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))

mean_cp <- mean(mcmc_output[,grep("cp[2]",param_names,fixed = T)])

loglik_mat <- mcmc_output[,grep("loglik",param_names)]
loglik_mat <- loglik_mat[,-grep("total_loglik",colnames(loglik_mat))]

waic <- loo::waic(loglik_mat)[["estimates"]][3, 1]
pml <- -2 * sum(log(nrow(loglik_mat)/colSums(1/exp(loglik_mat))))


png(paste0(pathway, "Common-Final-Weibull-model.png"), width = 10, height = 5, units = 'in', res = 300)
plot(fit.OS, xlim = c(0, max(data_jags$t_pred)), xlab = "Time (Years)", ylab = "S(t)", 
     main = "Weibull 1 change-point model (Equal final hazards)",
     sub=paste0("WAIC - ", round(waic, digits = 2),": Mean change-point - ", round(mean_cp, digits= 2) ))
lines(x =data_jags$t_pred , y =St_mat_quantile[2,1:(ncol(St_mat_quantile)/2)], col = "blue")
lines(x =data_jags$t_pred , y =St_mat_quantile[2,(ncol(St_mat_quantile)/2 +1 ):ncol(St_mat_quantile)], col = "red")
dev.off()