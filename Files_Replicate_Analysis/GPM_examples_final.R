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

library("ggmcmc", lib.loc = "~/R-packages/")
library("runjags", lib.loc = "~/R-packages/")
library("PiecewiseChangepoint", lib.loc = "~/R-packages/")
library("survival")

#Data

set.seed(123)
n_obs =300
n_events_req=300
max_time = 24 #months
rate = c(0.75,0.25)/12 #we want to report on months
t_change =12 #change-point at 12 months
df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                   num.breaks = length(t_change),rate = rate ,
                   t_change = t_change, max_time = max_time)

fit.km <- survfit(Surv(time,status)~1, data = df)


path<-"~/PhD/KM_Piecewise_Major_Review_final/"
Conditional_Death_df <- read.xlsx(paste0(path, "Conditional_Death_UK.xlsx"),1)

time_factor <- 12
time_horizon <- 100
colnames(Conditional_Death_df) <- c("age", "Male_cond_death", "Female_cond_death")

df_temp <- Conditional_Death_df
age_baseline_example <- 75
prop_male <- 0.5

df$age <- rnorm(n = nrow(df), mean = age_baseline_example, sd =10)
df$age <- round(pmin(df$age, 85), digits = 2)
df$sex <- rbinom(nrow(df), size = 1, prob = 0.5)

df_temp$Male_cond_haz <- -log(1-df_temp$Male_cond_death)/time_factor 
df_temp$Female_cond_haz <- -log(1-df_temp$Female_cond_death)/time_factor
df_temp <- df_temp %>% filter(age >= age_baseline_example & age <= time_horizon)

n_row_df_temp <- nrow(df_temp)
time_partial_vec <- rep((0:(time_factor-1))/time_factor)

df_temp <- do.call("rbind", replicate(time_factor, df_temp, simplify = FALSE)) %>%
  arrange(age)
df_temp$age <- df_temp$age+  rep(time_partial_vec,times = n_row_df_temp)


#Cohort - Static Gender proportion
df_temp[, "mix_haz_static"] <- df_temp[,"Male_cond_haz"]*prop_male + df_temp[,"Female_cond_haz"]*(1-prop_male)

#Cohort - Dynamic Gender proportion
df_temp$male_St <- exp(-cumsum(df_temp$Male_cond_haz))
df_temp$female_St <- exp(-cumsum(df_temp$Female_cond_haz))

df_temp$prop_male <- df_temp$male_St/(df_temp$male_St +df_temp$female_St)
df_temp[, "mix_haz_dynamic"] <- df_temp[,"Male_cond_haz"]*df_temp$prop_male  + df_temp[,"Female_cond_haz"]*(1-df_temp$prop_male)

#gmp_haz_vec_example = df_temp[, "mix_haz_static"]


GMP_haz <- GMP_cum_haz <- rep(NA,nrow(df))
for(i in 1:nrow(df)){
  age_start <- which.min(abs(df$age[i] - df_temp$age))
  age_selc <-  which.min(abs((df$age[i]+df$time[i]/time_factor) - df_temp$age))
  
  if(df$sex[i] == 1){
    GMP_haz[i] <- df_temp[age_selc, "Male_cond_haz"] 
    GMP_cum_haz[i] <- sum(df_temp[age_start:age_selc,"Male_cond_haz"])
    
  }else{
    GMP_haz[i] <- df_temp[age_selc, "Female_cond_haz"] #assume rounding cancels ##aproximation
    GMP_cum_haz[i] <- sum(df_temp[age_start:age_selc,"Female_cond_haz"])
  } 
}


data_jags <- list()
data_jags$N <- nrow(df)
data_jags$time <- df$time
data_jags$status <- df$status
data_jags$MAXX <- max(df$time)
data_jags$N_CP <- 1 #Number of changepoints

data_jags$GMP_haz <-GMP_haz
data_jags$GMP_cum_haz <-GMP_cum_haz
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
cum_haz_mean_pred2 <- -log(colMeans(St_all))

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
  cum_haz_mat2[i,] <- cum_haz_mean_pred2 + GMP_cum_haz_mat2[i,]
}

Surv_mat2 <- exp(-cum_haz_mat2)
St_final_mean_sim <- colMeans(Surv_mat2)

# 
# for(i in 1:data_jags$N){
#   GMP_curr <- matrix(GMP_cum_haz_mat[i,], ncol = length(GMP_cum_haz_mat[i,]), nrow = nrow(cum_haz_array2), byrow = T)
#   cum_haz_array2[,,i] <- cum_haz_pred2 +  GMP_curr
# }
# 
# Surv_array2 <- exp(-cum_haz_array2)
# Surv_mat2_ind <- apply(Surv_array2,c(1,2),mean)
# St_final2 <- colMeans(Surv_mat2_ind)

#abs(St_final2 -St_final_mean2) < 0.001

# Static/Dynamic based approach

St_mean <- colMeans(St_all)
cum_Haz_mean <- -log(St_mean)

#data_jags$t_pred assume this includes 0
final_haz_static <- df_temp[1:(max(data_jags$t_pred)+1),"mix_haz_static"]
final_haz_static[1] <- 0
final_haz_dynamic <- df_temp[1:(max(data_jags$t_pred)+1),"mix_haz_dynamic"]
final_haz_dynamic[1] <- 0

St_final_static <- exp(-(cum_Haz_mean +final_haz_static))
St_final_dynamic <- exp(-(cum_Haz_mean +final_haz_dynamic))


#Plot output -- 

png(paste0("Surv-plot-",round(age_baseline_example,digits = 0),".png"), width = 10, height = 4, units = 'in', res = 300)
plot(fit.km, xlim = c(0, max(data_jags$t_pred)), xlab = "Time (months)", ylab = "Survival")
points(x = data_jags$t_pred, y = St_final_static, col = "red")
points(x = data_jags$t_pred, y = St_final_dynamic, col = "blue")
points(x = data_jags$t_pred, y = St_final_mean_sim, col = "purple")
points(x = data_jags$t_pred, y = St_final_internal, col = "green")
legend("topright", legend=c("Cohort: Static Gender","Cohort: Dynamic Gender", "Simulation", "Internal Additive"),
       col=c("red", "blue", "purple", "green"), lty=1, cex=0.8)
dev.off()

