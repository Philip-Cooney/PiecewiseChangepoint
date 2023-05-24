jags.piecewise.expo_chng <-"
  data{
    for(i in 1:N){
    zero[i] <- 0
    }
    zero_prior <- 0
  }
  model {
  #Constant for zerors trick
  C <- 10000
  #Prior for change-point locations - Continous version of the Fearnhead Prior
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
    sd  <- 100 # ~ dunif(0.2,5)
    prec <- pow(sd, -2)
  
  for(k in 1:(N_CP+1)){
    beta_cp[k] ~ dnorm(0,prec)
    
  }
  
  # Model and likelihood
  for (i in 1:N) {
  
  for(k in 1:(N_CP+1)){
    #variable which gives the difference between the two intervals if time[i]>cp[k+1]
    #(i.e. cp[k+1] -cp[k]) or time between time[i] and cp[k]
    X[i,k] = max(min(time[i], cp[k+1]) - cp[k],0) 
   
    #Indicator variable which highlights which interval time is in 
    X_ind[i,k] = step(time[i]-cp[k])*step(cp[k+1]-time[i])
    
    param_1[i,k] <- exp(beta_cp[k]) #lambda
    
    log_haz_seg[i,k] <-  log(param_1[i,k])*X_ind[i,k]
    cum_haz_seg[i,k] <- param_1[i,k]*X[i,k]
    
  }
  
  log_haz[i] <- sum(log_haz_seg[i,]) 
  haz[i] <- exp(log_haz[i]) + GMP_haz[i]
  log_haz_final[i] <- log(haz[i])
  cum_haz[i] <- sum(cum_haz_seg[i,])
  cum_haz_gmp[i] <- cum_haz[i] + GMP_cum_haz[i]
  
  loglik[i] = status[i]*log_haz_final[i] - cum_haz_gmp[i]
  zero[i] ~ dpois(C - loglik[i])
  }
   
  #Predict Survival
  for(i in 1:length(t_pred)){
   for(k in 1:(N_CP+1)){
   
    X_pred[i,k] = max(min(t_pred[i], cp[k+1]) - cp[k],0) 
    X_ind_pred[i,k] = step(t_pred[i]-cp[k])*step(cp[k+1]-t_pred[i])
    
    param_1_pred[i,k] <- exp(beta_cp[k]) #lambda
    cum_haz_seg_pred[i,k] <- param_1_pred[i,k]*X_pred[i,k]
  }
  cum_haz_pred[i] <- sum(cum_haz_seg_pred[i,]) # Need to add the GMP cum_haz
    
  }
  total_loglik <- sum(loglik)
  }
"
#Load the following packages

library("runjags")
library("PiecewiseChangepoint")
library("survival")
library("xlsx")
library("rjags")

#Load the following packages

library("ggmcmc", lib.loc = "~/R-packages/")
library("runjags", lib.loc = "~/R-packages/")
library("PiecewiseChangepoint", lib.loc = "~/R-packages/")
library("survival")
library("xlsx")
library("runjags")
#Data

if(!exists("path")){
  stop("Set path variable")
}

set.seed(123)
n_obs =300
max_time = 24 #months
rate = c(0.75,0.25)/12 #we want to report on months
t_change =12 #change-point at 12 months


time_factor <- 12
time_horizon <- 100
Conditional_Death_df <- read.xlsx(paste0(path, "Conditional_Death_UK.xlsx"),1)
colnames(Conditional_Death_df) <- c("age", "Male_cond_death", "Female_cond_death")
Conditional_Death_df <- rbind(Conditional_Death_df,c(101, .99999, .99999)) #All patients die at 101
age_baseline_example <-80
prop_male <- 0.5
max_age <- 90
#Have to have a minimum of 60 months with GPM survival avialbe (if t_pred is 60)


rpwexp<-function (n, lam, s){
  U = runif(n, 0, 1)
  X = rep(NA, n)
  haz_seg <- diff(c(0, s)) * lam[-length(lam)]
  cum_haz <- cumsum(haz_seg)
  St_thres <- exp(-cum_haz)
  for (i in 1:n) {
    int <- which(U[i] < St_thres)
    if (length(int) == 0) {
      X[i] = qexp(U[i], rate = lam[1], lower.tail = F)
    }
    else {
      X[i] = s[max(int)] + qexp(U[i]/St_thres[max(int)], 
                                rate = lam[max(int) + 1], lower.tail = F)
    }
  }
  return(data.frame(St = U, time = X))
}
n_sims <- 10000

df_pwexp_sims <- rpwexp(n_sims, rate, t_change)
df_pwexp_sims[which(df_pwexp_sims[,"time"] > 600), "time"] <- 600
df_pwexp_sims[which(df_pwexp_sims[,"time"] >= 600), "St"] <- 0
#600 = 50 years

df_pwexp_sims$Cum_haz_pwexp <- -log(df_pwexp_sims$St)

age_vec <- rnorm(n = n_obs, mean = age_baseline_example, sd =10)
age_vec <- round(pmin(age_vec,max_age ), digits = 2)
sex_vec <- rbinom(n_obs, size = 1, prob = 0.5)


df_temp <- Conditional_Death_df
df_temp$Male_cond_haz <- -log(1-df_temp$Male_cond_death)/time_factor 
df_temp$Female_cond_haz <- -log(1-df_temp$Female_cond_death)/time_factor

n_row_df_temp <- nrow(df_temp)
time_partial_vec <- rep((0:(time_factor-1))/time_factor)

df_temp <- do.call("rbind", replicate(time_factor, df_temp, simplify = FALSE)) %>%
  arrange(age)
df_temp$age <- df_temp$age+  rep(time_partial_vec,times = n_row_df_temp)

mean_age_baseline_example <- mean(age_vec)
#Cohort approach
df_temp_cohort <- df_temp %>% filter(age >= mean_age_baseline_example & age <= time_horizon)
#Cohort - Static Gender proportion
df_temp_cohort[, "mix_haz_static"] <- df_temp_cohort[,"Male_cond_haz"]*prop_male + df_temp_cohort[,"Female_cond_haz"]*(1-prop_male)

#Cohort - Dynamic Gender proportion
df_temp_cohort$male_St <- exp(-cumsum(df_temp_cohort$Male_cond_haz))
df_temp_cohort$female_St <- exp(-cumsum(df_temp_cohort$Female_cond_haz))

df_temp_cohort$prop_male <- df_temp_cohort$male_St/(df_temp_cohort$male_St +df_temp_cohort$female_St)
df_temp_cohort[, "mix_haz_dynamic"] <- df_temp_cohort[,"Male_cond_haz"]*df_temp_cohort$prop_male  + df_temp_cohort[,"Female_cond_haz"]*(1-df_temp_cohort$prop_male)

GMP_cum_haz_sim <- matrix(nrow = n_sims, ncol = n_obs)

for(i in 1:n_sims){
  
  for(j in 1:n_obs){
    age_start <- which.min(abs(age_vec[j] - df_temp$age))
    age_selc <-  which.min(abs((age_vec[j]+ df_pwexp_sims$time[i]/time_factor) - df_temp$age))
    
    if(sex_vec[j] == 1){
      
      GMP_cum_haz_sim[i,j] <- sum(df_temp[age_start:age_selc,"Male_cond_haz"], na.rm = T)
    }else{
      GMP_cum_haz_sim[i,j] <- sum(df_temp[age_start:age_selc,"Female_cond_haz"], na.rm = T)
    } 
    
  }
  
}


#?sweep
#sweep(data, 1, to_add, "+")
final_cum_haz <- sweep(GMP_cum_haz_sim,1,df_pwexp_sims$Cum_haz_pwexp, FUN = "+")
final_St_mat <- exp(-final_cum_haz)

unif_vec <- runif(n_obs)
time_vec_final <- rep(NA, n_obs)
time_surv <- df_pwexp_sims$time

for(i in 1:n_obs){
  time_vec_final[i] <- time_surv[which.min(abs(final_St_mat[,i]-unif_vec[i]))] 
}

df <- data.frame(time_event = time_vec_final) %>% mutate(time = ifelse(time_event < max_time,
                                                                       time_event, max_time),
                                                         status =ifelse(time_event < max_time,
                                                                        1, 0),
                                                         age = age_vec,
                                                         sex = sex_vec)

plot(survfit(Surv(time,status)~1, data = df))

#This graph can appear a bit paradoxical as the blue curve might not go below a certain point.
#This is because we don't model after 100 years due to lack of data.
# For example if someone is 95 they have a 9% chance of living to 100.
# While a person who is only 73 which have a lower than 9% chance of living to 100.
plot(y = final_St_mat[,which.min(df$age)], x =df_pwexp_sims$time/time_factor, xlab = "Time in Years", ylab = "Survival" )
points(y = final_St_mat[,which.max(df$age)], x =df_pwexp_sims$time/time_factor, col = "blue")


split_df <- survSplit(Surv(time, status) ~., df,
                      cut=c(5,10,15, 20, 50), episode ="period")

split_df %>% group_by(period) %>% summarize(mean_age = mean(age),
                                            sd_age = sd(age),
                                            n = n(),
                                            lcl = mean_age -1.96*sd_age/sqrt(n),
                                            ucl = mean_age +1.96*sd_age/sqrt(n))

fit.km <- survfit(Surv(time,status)~1, data = df)

cor(df$time_event, df$age)
plot(df$time_event, df$age)
plot(fit.km)

GMP_haz_jags <- GMP_cum_haz_jags <- rep(NA,nrow(df))
for(i in 1:nrow(df)){
  age_start <- which.min(abs(df$age[i] - df_temp$age))
  age_selc <-  which.min(abs((df$age[i]+df$time[i]/time_factor) - df_temp$age))
  
  if(df$sex[i] == 1){
    GMP_haz_jags[i] <- df_temp[age_selc, "Male_cond_haz"] 
    GMP_cum_haz_jags[i] <- sum(df_temp[age_start:age_selc,"Male_cond_haz"])
    
  }else{
    GMP_haz_jags[i] <- df_temp[age_selc, "Female_cond_haz"] #assume rounding cancels ##aproximation
    GMP_cum_haz_jags[i] <- sum(df_temp[age_start:age_selc,"Female_cond_haz"])
  } 
}


data_jags <- list()
data_jags$N <- nrow(df)
data_jags$time <- df$time
data_jags$status <- df$status
data_jags$MAXX <- max(df$time)
data_jags$N_CP <- 1 #Number of changepoints

data_jags$GMP_haz <-GMP_haz_jags
data_jags$GMP_cum_haz <-GMP_cum_haz_jags
data_jags$t_pred <- c(0:60) #Time horizon to predict

#Change-inits
inits <- function(){
  list(unif = runif(data_jags$N_CP,max = data_jags$MAXX/data_jags$N_CP))
}


expo.mod_1chng <- runjags::run.jags(
  model = jags.piecewise.expo_chng,
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp", "cum_haz_pred", "cum_haz", "total_loglik"), 
  sample=5000, 
  thin = 1, 
  burnin = 1000,
  inits = inits, 
  method ='rjparallel')


add.summary(expo.mod_1chng, "cp")
mcmc_output <- rbind(expo.mod_1chng$mcmc[[1]],expo.mod_1chng$mcmc[[2]])
param_names <- colnames(expo.mod_1chng$mcmc[[1]]) 

GMP_cum_haz_mat <- matrix(nrow = nrow(df), ncol = length(data_jags$t_pred))
for(i in 1:nrow(df)){
  age_start <- which.min(abs(df$age[i] - df_temp$age))
  
  if(df$sex[i] == 1){
    GMP_cum_haz_mat[i,] <- cumsum(df_temp[age_start+data_jags$t_pred,"Male_cond_haz"])
  }else{
    GMP_cum_haz_mat[i,] <- cumsum(df_temp[age_start+data_jags$t_pred,"Female_cond_haz"])
  } 
  
  GMP_cum_haz_mat[,1] <- 0 #assume zero cumulative hazard at age at baseline
}

cum_haz_pred <- mcmc_output[,grep("cum_haz_pred[",param_names, fixed = T)]
cum_haz_pred <- cum_haz_pred[,1:length(data_jags$t_pred)]
cum_haz_array <- array(dim = c(nrow(cum_haz_pred), ncol(cum_haz_pred),data_jags$N))

for(i in 1:data_jags$N){
  
  GMP_curr <- matrix(GMP_cum_haz_mat[i,], ncol = length(GMP_cum_haz_mat[i,]), nrow = nrow(cum_haz_array), byrow = T)
  cum_haz_array[,,i] <- cum_haz_pred +  GMP_curr
  # if(any(is.na(cum_haz_array[,,i]))){
  #   print(paste0("NA in ", i))
  # }
  
}
Surv_array <- exp(-cum_haz_array)
St_final_internal <- colMeans(apply(Surv_array,c(1,2),mean))

Collapsing_Model <- collapsing.model(df,
                                     n.iter = 5750,
                                     burn_in = 750,
                                     n.chains = 2,
                                     timescale = "months")


St_all <- t(get_Surv(Collapsing_Model, time = data_jags$t_pred, chng.num = 1))

#Simulation based approach

#We can just take the mean of the survival 

St_mean <- colMeans(St_all)
cum_Haz_mean <- -log(St_mean)

#rather than the posterior distbution of the survival -which will give same result to three decimal places
# cum_haz_pred2 <- -log(St_all)
# cum_haz_array2 <- array(dim = c(nrow(cum_haz_pred2), ncol(cum_haz_pred2),data_jags$N))

#we don't add GMP until after the trial ends

t_max_round <- max(df$time)
index_no_GPM <- which(data_jags$t_pred <=t_max_round)


GMP_cum_haz_mat2 <- matrix(nrow = nrow(df), ncol = length(data_jags$t_pred))
for(i in 1:nrow(df)){
  age_start <- which.min(abs(df$age[i] - df_temp$age))
  df_temp2 <- df_temp[age_start+data_jags$t_pred,]
  df_temp2$Male_cond_haz[index_no_GPM] <- 0 
  df_temp2$Female_cond_haz[index_no_GPM] <- 0 
  if(df$sex[i] == 1){
    GMP_cum_haz_mat2[i,] <- cumsum(df_temp2[,"Male_cond_haz"])
  }else{
    GMP_cum_haz_mat2[i,] <- cumsum(df_temp2[,"Female_cond_haz"])
  } 
}



cum_haz_mat2 <- GMP_cum_haz_mat2
for(i in 1:data_jags$N){
  cum_haz_mat2[i,] <- cum_Haz_mean + GMP_cum_haz_mat2[i,]
}

Surv_mat2 <- exp(-cum_haz_mat2)
St_final_mean_sim <- colMeans(Surv_mat2)

# Static/Dynamic based approach
#Dont' include GPM before extrapolation
df_temp_cohort[index_no_GPM, c("mix_haz_static","mix_haz_dynamic")] <- 0


final_haz_static <- df_temp_cohort[1:(max(data_jags$t_pred)+1),"mix_haz_static"]
final_haz_dynamic <- df_temp_cohort[1:(max(data_jags$t_pred)+1),"mix_haz_dynamic"]

St_final_static <- exp(-(cum_Haz_mean +cumsum(final_haz_static)))
St_final_dynamic <- exp(-(cum_Haz_mean +cumsum(final_haz_dynamic)))


#Plot output -- 


png(paste0("Surv-plot-",round(mean_age_baseline_example ,digits = 0),".png"), width = 10, height = 4, units = 'in', res = 300)
plot(fit.km, xlim = c(0, max(data_jags$t_pred)), xlab = "Time (months)", ylab = "Survival")
points(x = data_jags$t_pred, y = St_final_static, col = "red")
points(x = data_jags$t_pred, y = St_final_dynamic, col = "blue")
points(x = data_jags$t_pred, y = St_final_mean_sim, col = "purple")
points(x = data_jags$t_pred, y = St_final_internal, col = "green")
legend("topright", legend=c("Cohort: Static Gender","Cohort: Dynamic Gender", "Simulation", "Internal Additive"),
       col=c("red", "blue", "purple", "green"), lty=1, cex=0.8)
dev.off()
