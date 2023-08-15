pathway <- "~/Change-Point Simulation Studies/Simulation Study 2023/"

library("rjags")
library("survival")
library("flexsurv")
#library("survminer")
library("runjags")
library("dplyr")
library("xlsx")
#library("bshazard")

source(paste0(pathway, "Jags_codes.R"))
mod <- "wei"

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
E1690.dat$AGE <- as.numeric(E1690.dat$AGE)
E1690.dat$AGE[is.na(E1690.dat$AGE)] <- mean(E1690.dat$AGE, na.rm = T)
E1690.dat$AGE_scale <- as.numeric(scale(E1690.dat$AGE))
E1690.dat$AGE_scale_abs <- abs(E1690.dat$AGE_scale)
# fit.OS <- survfit(Surv(SURVTIME, SURVCENS)~TRT, 
#                   data = E1690.dat[which(E1690.dat$STUDY == "1684"),])



n.samp <- 500
n.burnin <- 100
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

t_pred <-seq(0, 15, by = 0.1)
spline_St <- summary(fs3,t = t_pred )
spline_St_trt <- spline_St[[1]]
spline_St_comp <- spline_St[[2]]


png(paste0(pathway, "Royston-Parmar-Spline-Model.png"), width = 10, height = 5, units = 'in', res = 300)
plot(fit.OS, xlim = c(0,15),main = "Survival estimated by Royston-Parmar model")
lines(spline_St_trt$time, spline_St_trt$est, lty=1, lwd=1, col = "blue")
lines(spline_St_comp$time, spline_St_comp$est, lty=1, lwd=1, col = "red")
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



Surv_extrap <- data.frame(time = rep(t_pred,2),
                          Survival =c(summary(fs3, t= t_pred)[[1]]$est,summary(fs3, t= t_pred)[[2]]$est), 
                          Treatment= c(rep("INF",length(t_pred)),
                                       rep("OBS",length(t_pred))))


p <- ggplot(plot_data, aes(x = time, y = Survival, color = Treatment))+
  geom_step( size = 1)+
  ylim(c(0,1))+
  xlim(c(0,15))+
  xlab("Time (Months)")+
  ggtitle("Royston-Parmar Spline Survival Model")+
  geom_line(data = Surv_extrap,
            aes(x = time, y = Survival,color = Treatment ))+
  theme_bw()
ggsave(paste0(pathway, "Royston-Parmar-Spline-Model.png"), width = 10, height = 5, units = 'in')

  
if(FALSE){

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

}
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

data_jags$t_pred <- seq(0,15, by = 0.2)

#X_mat_fixed is just a dummy matrix
#PH Model with Changepoint
data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N), E1690.dat$TRT-1, E1690.dat$AGE_scale),  ncol = data_jags$ncovar_fixed)
#independent parameters 
data_jags$X_mat_trt <- matrix(c(rep(1, data_jags$N),rep(0, data_jags$N)),  ncol = 2)


data_jags$scenario_common <- c(0,0)
data_jags$scenario <- 1  #
data_jags$cp_fix <- c(0, 1, 5) 
data_jags$cp_fixed <- 0 


#Change-inits
inits <- function(){
  list(unif = runif(data_jags$N_CP,max = data_jags$MAXX/data_jags$N_CP)
       #,
       #beta_cp_0 = c(0),
       #beta_cp = rnorm(data_jags$N_CP+1)
  )
}


if(FALSE){
  
mod.cp <- runjags::run.jags(
  model = get(paste0("jags.piecewise.",mod,"_chng")),
  #model = jags.piecewise.expo_chng,
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp","beta_cp_ancs", "prec","sd", "cum_haz_pred", "cum_haz", "total_loglik","beta_covar_fixed"), 
  sample=n.samp, 
  thin = n.thin, 
  burnin = n.burnin,
  inits = inits, 
  method ='rjparallel')

add.summary(mod.cp, "cp")

param_names <- colnames(mod.cp$mcmc[[1]])
mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])
St_mat <- exp(-mcmc_output[,grep("cum_haz_pred",param_names)])
St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))

plot(fit.OS, xlim = c(0,15))
lines(x =data_jags$t_pred , y =St_mat_quantile[2,1:(ncol(St_mat_quantile)/2)], col = "blue")
lines(x =data_jags$t_pred , y =St_mat_quantile[2,(ncol(St_mat_quantile)/2 +1 ):ncol(St_mat_quantile)], col = "red")


}

data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N*2), E1690.dat$AGE_scale),  ncol = data_jags$ncovar_fixed)
#data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N*2)),  ncol = 2)
data_jags$X_mat_trt <- matrix(c(rep(1, data_jags$N), E1690.dat$TRT-1),  ncol = 2)


mod.cp <- runjags::run.jags(
  #model = jags.piecewise.expo_chng,
  model = get(paste0("jags.piecewise.",mod,"_chng")),
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp", "prec","sd", "cum_haz_pred", "cum_haz", "total_loglik","beta_covar_fixed", "loglik"), 
  sample=n.samp, 
  thin = n.thin, 
  burnin = n.burnin,
  inits = inits, 
  method ='rjparallel')

param_names <- colnames(mod.cp$mcmc[[1]])
model_output_compare <- mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])
St_mat <- exp(-mcmc_output[,grep("cum_haz_pred",param_names)])
St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))


loglik_mat <- mcmc_output[,grep("loglik",param_names)]
loglik_mat <- loglik_mat[,-grep("total_loglik",colnames(loglik_mat))]

waic <- loo::waic(loglik_mat)[["estimates"]][3, 1]
pml <- -2 * sum(log(nrow(loglik_mat)/colSums(1/exp(loglik_mat))))
#add.summary(mod.cp, "beta_cp")

summary_output <- mcmc_output[,grep("beta_cp|cp|beta_covar_fixed",param_names)]
summary_output <- cbind(summary_output[,2],exp(summary_output[, c(6,7, 10) ]))
colMeans(summary_output)

#summary_output[, c(2,3, 10)] <- exp(summary_output[, c(2,3, 10)])
summary_output_df <- data.frame(Model = c("Changepoint", "HR 1", "HR 2", "HR Age"))
summary_output_df <- cbind(summary_output_df,t(apply(summary_output, 2, quantile, probs = c(0.025,0.5,0.975))))

write.csv(summary_output_df,paste0(pathway, "summary_output_df_weib_cox.csv"))
           

png(paste0(pathway, "PH-Weibull-model.png"), width = 10, height = 5, units = 'in', res = 300)
plot(fit.OS, xlim = c(0, max(data_jags$t_pred)), xlab = "Time (Years)", ylab = "S(t)", 
     main = "Weibull 1 change-point model (Proportional Hazards)",
     sub=paste0("WAIC - ", round(waic, digits = 2)))
lines(x =data_jags$t_pred , y =St_mat_quantile[2,1:(ncol(St_mat_quantile)/2)], col = "blue")
lines(x =data_jags$t_pred , y =St_mat_quantile[2,(ncol(St_mat_quantile)/2 +1 ):ncol(St_mat_quantile)], col = "red")
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
           Survival =St_mat_quantile[2,], 
           Treatment= c(rep("OBS",length(data_jags$t_pred)),
                        rep("INF",length(data_jags$t_pred)))
           )


mean.cp <- mean(mcmc_output[,"cp[2]"])
dev <-abs(Surv_extrap$time - mean.cp)

p <- ggplot(plot_data, aes(x = time, y = Survival, color = Treatment))+
   geom_step( size = 1)+
   ylim(c(0,1))+
   xlab("Time (Years)")+
  ggtitle("Weibull 1 change-point model (PH change)")+
  geom_line(data = Surv_extrap,
            aes(x = time, y = Survival,color = Treatment ))+
  geom_point(data = data.frame(time = mean.cp,
                               Survival = mean(pull(Surv_extrap,Survival)[which(dev == min(dev))])),
             aes(x = time, Survival), shape = 23, fill = "green",
             color = "darkred", size = 4, inherit.aes = F, stroke  = 1.5)+
  theme_bw()
ggsave(paste0(pathway, "PH-Weibull-model.png"), width = 10, height = 5, units = 'in')



#Scenario 2 -- Common Hazards for final interval
data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N), E1690.dat$TRT-1, E1690.dat$AGE_scale),  ncol = data_jags$ncovar_fixed)


data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N*2),E1690.dat$AGE_scale),  ncol = data_jags$ncovar_fixed)
X_mat_trt <- matrix(c(rep(1, data_jags$N), E1690.dat$TRT-1),  ncol = 2)

data_jags$X_mat_trt <-X_mat_trt

data_jags$scenario_common <- c(0,1)
data_jags$scenario <- 1  #
data_jags$cp_fix <- c(0, 1, 5) 
data_jags$cp_fixed <- 0 

mod.cp <- runjags::run.jags(
  #model = jags.piecewise.expo_chng_common_haz,
  model = get(paste0("jags.piecewise.",mod,"_chng")),
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp", "prec","sd", "cum_haz_pred", "cum_haz", "total_loglik","beta_covar_fixed", "loglik"), 
  sample=n.samp, 
  thin = n.thin, 
  burnin = n.burnin,
  inits = inits, 
  method ='rjparallel')


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
                          Survival =St_mat_quantile[2,], 
                          Treatment= c(rep("OBS",length(data_jags$t_pred)),
                                       rep("INF",length(data_jags$t_pred))))


mean.cp <- mean(mcmc_output[,"cp[2]"])
dev <-abs(Surv_extrap$time - mean.cp)

p <- ggplot(plot_data, aes(x = time, y = Survival, color = Treatment))+
  geom_step( size = 1)+
  ylim(c(0,1))+
  xlab("Time (Years)")+
  ggtitle("Weibull 1 change-point model (Equal final hazards)")+
  geom_line(data = Surv_extrap,
            aes(x = time, y = Survival,color = Treatment ))+
  geom_point(data = data.frame(time = mean.cp,
                               Survival = mean(pull(Surv_extrap,Survival)[which(dev == min(dev))])),
             aes(x = time, Survival), shape = 23, fill = "green",
             color = "darkred", size = 4, inherit.aes = F, stroke  = 1.5)+
  theme_bw()
ggsave(paste0(pathway, "Common-Final-Weibull-model.png"), width = 10, height = 5, units = 'in')




#Scenario 3

data_jags$X_mat_fixed <- matrix(c(rep(0, data_jags$N*2),E1690.dat$AGE_scale),  ncol = data_jags$ncovar_fixed)
X_mat_trt<- matrix(c(rep(1, data_jags$N), E1690.dat$TRT-1),  ncol = 2)
data_jags$X_mat_trt <-X_mat_trt
#data_jags$t_pred <- 0:20

mod.cp <- runjags::run.jags(
  model = get(paste0("jags.piecewise.",mod,"_chng_wane")),
  #model = jags.piecewise.expo_wane,
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp", "total_loglik", "cum_haz", "lower_int", "upper_int", "beta_cp_anc", "loglik"), 
  sample=n.samp, 
  thin = n.thin, 
  burnin = n.burnin,
  inits = inits, 
  method ='rjparallel')

param_names <- colnames(mod.cp$mcmc[[1]]) #rownames(expo.mod_1chng$BUGSoutput$summary)
mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])


loglik_mat <- mcmc_output[,grep("loglik",param_names)]
loglik_mat <- loglik_mat[,-grep("total_loglik",colnames(loglik_mat))]

waic <- loo::waic(loglik_mat)[["estimates"]][3, 1]
pml <- -2 * sum(log(nrow(loglik_mat)/colSums(1/exp(loglik_mat))))

St_mat <- exp(-mcmc_output[,grep("cum_haz[",param_names, fixed = T)])
St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))

plot(fit.OS, xlim = c(0, max(data_jags$t_pred)))
points(x = E1690.dat$SURVTIME[E1690.dat$TRT == 1], y = St_mat_quantile[2,E1690.dat$TRT == 1], col = "red")
points(x = E1690.dat$SURVTIME[E1690.dat$TRT == 2], y = St_mat_quantile[2,E1690.dat$TRT == 2], col = "blue")

t_pred <- seq(0, 50, by = 0.1)
HR_mat <- St_pred2 <- St_pred <- matrix(ncol = length(t_pred), nrow = nrow(mcmc_output))

  incomplete_gamma_upper<- function(s,x){
    exp(lgamma(s))*(1-pgamma(x,s,1))
  }
  #Write it out in integral calc
  cum_haz_seg <- function(t_0,t_1,a, m, initial_HR, lambda_wane  ){
    constant <- (initial_HR-1)*exp(t_0*lambda_wane)
    int_1 <- (t_0^a)/a-(incomplete_gamma_upper(a,lambda_wane*t_0)*(t_0^a)*constant)/((lambda_wane*t_0)^a)
    int_2 <- (t_1^a)/a-(incomplete_gamma_upper(a,lambda_wane*t_1)*(t_1^a)*constant)/((lambda_wane*t_1)^a)
    return((int_2-int_1)*a*m)
  }
  
  


if( FALSE){

  #Double check and plot out the function for exponential model
q <- 2
t_0 <- q
t_1 = 5
m = 2
initial_HR = 0.7
lambda_wane = 0.5
#cum_haz_seg(t_0 = t_0, t_1 = t_1, m = m, initial_HR = initial_HR, lambda_wane = lambda_wane, q = q)
lambda = 0.8

#Important plot
t_vec <- 0:15
plot(t_vec, y = lambda*(1-(1-initial_HR)*exp(-lambda_wane*(t_vec-q))), xlab = "Time", ylab = "Hazard")
abline(v = q)#, h = initial_HR*lambda)
#So we need to integrate to the right of the change-point.


haz_waning <- function(m, initial_HR, t, q){
  m*(1-(1-initial_HR)*exp(-lambda_wane*(t-q)))
}

integrate(haz_waning, lower = q, upper = t_1, m = m ,initial_HR= initial_HR, q = q)

}


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
      
      initial_HR <-  exp(mcmc_output[i,"beta_cp[1,2]"])
      lambda_wane <- exp(mcmc_output[i,"beta_cp[2,2]"])
      HR_final <-1-(1-initial_HR)*exp(-lambda_wane*X_eval)
      HR_mat[i,j] <-  HR_final
      
      t_0 <- mcmc_output[i,"cp[2]"]
      t_1 <- t_pred[j]
      
      St_pred[i,j] <- exp(-(-log(St_pred_cp) +cum_haz_seg(t_0 = t_0, t_1 = t_1,a = a2, m = m2, initial_HR = initial_HR, lambda_wane = lambda_wane)))
      St_pred2[i,j] <- exp(-(-log(St_pred2_cp) +flexsurv:::HweibullPH(t_1, scale = m2, shape = a2)-flexsurv:::HweibullPH(t_0, scale = m2, shape = a2)))
      
    }
  }
}
mean_cp <- mean(mcmc_output[,grep("cp[2]",param_names,fixed = T)])

png(paste0(pathway, "Converge-Haz-Weibull-model.png"), width = 10, height = 5, units = 'in', res = 300)
plot(fit.OS, xlim = c(0, max(data_jags$t_pred)), xlab = "Time (Years)", ylab = "S(t)", 
     main = "Weibull 1 change-point model (Converge Haz hazards)",
     sub=paste0("WAIC - ", round(waic, digits = 2),": Mean change-point - ", round(mean_cp, digits= 2) ))
lines(y = colMeans(St_pred), t_pred, col = "blue")
lines(y = colMeans(St_pred2), t_pred, col = "red")
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

Surv_extrap <- data.frame(time = rep(t_pred,2),
                          Survival =c(colMeans(St_pred),colMeans(St_pred2)), 
                          Treatment= c(rep("INF",length(t_pred)),
                                       rep("OBS",length(t_pred))))



mean.cp <- mean(mcmc_output[,"cp[2]"])
dev <-abs(Surv_extrap$time - mean.cp)

p <- ggplot(plot_data, aes(x = time, y = Survival, color = Treatment))+
  geom_step( size = 1)+
  ylim(c(0,1))+
  xlim(c(0,15))+
  xlab("Time (Years)")+
  ggtitle("Weibull 1 change-point model (Converging Hazards)")+
  geom_line(data = Surv_extrap,
            aes(x = time, y = Survival,color = Treatment ))+
  geom_point(data = data.frame(time = mean.cp,
                               Survival = mean(pull(Surv_extrap,Survival)[which(dev == min(dev))])),
             aes(x = time, Survival), shape = 23, fill = "green",
             color = "darkred", size = 4, inherit.aes = F, stroke  = 1.5)+
  theme_bw()
ggsave(paste0(pathway, "Converge-Haz-Weibull-model.png"), width = 10, height = 5, units = 'in')





quantile_HR <- apply(HR_mat, 2, quantile, probs = c(0.025, 0.5,0.975))

png(paste0(pathway, "Converge-HR-Weibull-model.png"), width = 10, height = 5, units = 'in', res = 300)
plot(x = t_pred, y =quantile_HR[2,] ,xlim = c(0, 10), ylim = c(0,1.5), type = "l", 
     xlab = "Time (Years)", ylab = "HR")
lines(x = t_pred, y =quantile_HR[1,],lty= 2)
lines(x = t_pred, y =quantile_HR[3,], lty = 2)
dev.off()

quantile_HR2 <- data.frame(time =t_pred, t(quantile_HR))

ggplot(quantile_HR2,
       aes(x = time, y = X50.))+
  geom_line()+
  geom_line(aes(x = time, y = X2.5.), linetype = "dashed" )+
  geom_line(aes(x = time, y = X97.5.), linetype = "dashed" )+
  ylim(c(0,1.5))+
  xlim(c(0 ,5))+
  xlab("Time (Years)")+
  ylab("Hazard Ratio")+
  ggtitle("Hazard Ratio for the Converging Hazards Model")+
  theme_bw()
ggsave(paste0(pathway, "Converge-HR-Weibull-model.png"), width = 10, height = 5, units = 'in')





#https://www.ispor.org/docs/default-source/euro2022/6367isporkamgarposter21oct2022-pdf.pdf?sfvrsn=f86dde27_0



# Treatment Delay
#TA347 

path_mod_eval <- pathway

# TA347_PFS_N_df <- read.table(paste0(path_mod_eval,"TA347_N_PFS_Initial","/IPDdata.txt"), header = T)
# TA347_PFS_N_df$arm <- 1#"Nintedanib"
# TA347_PFS_Doce_df <- read.table(paste0(path_mod_eval,"TA347_Doce_PFS_Initial","/IPDdata.txt"), header = T)
# TA347_PFS_Doce_df$arm <- 0#"Docetaxel"
TA347_OS_N_df <- read.table(paste0(path_mod_eval,"TA347_N_OS_Initial","/IPDdata.txt"), header = T)
TA347_OS_N_df$arm <- 1#"Nintedanib"
TA347_OS_Doce_df <- read.table(paste0(path_mod_eval,"TA347_Doce_OS_Initial","/IPDdata.txt"), header = T)
TA347_OS_Doce_df$arm <- 0#"Docetaxel"


TA347_df_all <- rbind(TA347_OS_N_df,TA347_OS_Doce_df)
#TA347_df_all <- rbind(TA347_PFS_N_df,TA347_PFS_Doce_df)

TA347_df_all$time <- TA347_df_all$time
TA347.km <- survfit(Surv(time,status)~arm, data = TA347_df_all)


if(FALSE){
  survp <- ggsurvplot(TA347.km, legend.labs = c("Docetaxel", "Nintedanib + Docetaxel"), xlab = "Time (Months)", conf.int = F)
  ggsave(file = paste0(pathway, "TA347-Kaplan-Meier.png"), width = 10)
  ggsurvplot(TA347.km, fun = "cumhaz")
  
}


data_jags <- list()
data_jags$N <- nrow(TA347_df_all)
data_jags$time <- TA347_df_all$time
data_jags$status <- TA347_df_all$status
data_jags$MAXX <- max(TA347_df_all$time)
data_jags$N_CP <- 1
data_jags$ncovar <- 2
data_jags$index_pred <- c(min(which(TA347_df_all$arm == 1)),min(which(TA347_df_all$arm == 0))) #needs to be a vector 


data_jags$ncovar_fixed <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix
data_jags$ncovar_cp <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix

data_jags$t_pred <- seq(0,50, by = 0.01)


data_jags$X_mat_trt <-X_mat_trt

data_jags$scenario_common <- c(0,1)
data_jags$scenario <- 1  #
data_jags$cp_fix <- c(0, 1, 5) 
data_jags$cp_fixed <- 0 

data_jags$X_mat_trt<- matrix(c(rep(1, data_jags$N), TA347_df_all$arm),  ncol = 2)
data_jags$X_mat_fixed<- matrix(c(rep(2, data_jags$N*2)),  ncol = 2)

data_jags$scenario_common <- c(1,0)
data_jags$scenario <- 1  #
data_jags$cp_fix <- c(0, 1, 5) 
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
  monitor = c("cp", "beta_cp","beta_cp_ancs", "prec","sd", "cum_haz_pred", "cum_haz", "total_loglik","beta_covar_fixed", "loglik"), 
  sample=n.samp, 
  thin = n.thin, 
  burnin = n.burnin,
  inits = inits, 
  method ='rjparallel')

param_names <- colnames(mod.cp$mcmc[[1]])
mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])
St_mat <- exp(-mcmc_output[,grep("cum_haz_pred[",param_names, fixed = T)])
#St_mat <- exp(-mcmc_output[,grep("cum_haz[",param_names, fixed = T)])
St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))

add.summary(mod.cp, "cp")
add.summary(mod.cp, "beta_covar_fixed")
add.summary(mod.cp, "total_loglik")

loglik_mat <- mcmc_output[,grep("loglik",param_names)]
loglik_mat <- loglik_mat[,-grep("total_loglik",colnames(loglik_mat))]


mean_cp <- mean(mcmc_output[,grep("cp[2]",param_names,fixed = T)])
waic <- loo::waic(loglik_mat)[["estimates"]][3, 1]
pml <- -2 * sum(log(nrow(loglik_mat)/colSums(1/exp(loglik_mat))))


png(paste0(pathway, "Tretament-Delay-TA347.png"), width = 10, height = 5, units = 'in', res = 300)
plot(TA347.km, xlim = c(0, max(data_jags$t_pred)), ylab = "Survival", xlab = "Time (Months)",
     sub=paste0("WAIC - ", round(waic, digits = 2),": Mean change-point - ", round(mean_cp, digits= 2)))
 lines(x =data_jags$t_pred , y =St_mat_quantile[2,1:(ncol(St_mat_quantile)/2)], col = "blue")
 lines(x =data_jags$t_pred , y =St_mat_quantile[2,(ncol(St_mat_quantile)/2 +1 ):ncol(St_mat_quantile)], col = "red")
 legend("topright", legend=c("Nintedanib + Docetaxel", "Docetaxel"),
        col=c("blue", "red"), lty=1, cex=0.8)
#points(x = TA347_df_all$time[TA347_df_all$arm == 1], y = St_mat_quantile[2,TA347_df_all$arm == 1], col = "red")
#points(x = TA347_df_all$time[TA347_df_all$arm == 0], y = St_mat_quantile[2,TA347_df_all$arm == 0], col = "blue")
dev.off()
 


s <- summary(TA347.km, times = TA347.km$time, extend = T)

plot_data <- tibble(
  'time' = s$time,
  'n.risk' = s$n.risk,
  'n.event' = s$n.event,
  'n.censor' = s$n.censor,
  'Survival' = s$surv,
  'std.error' = s$std.err,
  'Treatment' = s$strata
)
plot_data <- plot_data %>% mutate(Treatment = ifelse(Treatment == "arm=1", "Nintedanib + Docetaxel", "Docetaxel")) 

Surv_extrap <- data.frame(time = rep(data_jags$t_pred,2),
                          Survival =c(St_mat_quantile[2,1:(ncol(St_mat_quantile)/2)],
                                      St_mat_quantile[2,(ncol(St_mat_quantile)/2 +1 ):ncol(St_mat_quantile)]), 
                          Treatment= c(rep("Nintedanib + Docetaxel",length(data_jags$t_pred)),
                                       rep("Docetaxel",length(data_jags$t_pred))))
mean.cp <- mean(mcmc_output[,"cp[2]"])
dev <-abs(Surv_extrap$time - mean.cp)

p <- ggplot(plot_data, aes(x = time, y = Survival, color = Treatment))+
  geom_step( size = 1)+
  ylim(c(0,1))+
  xlim(c(0,50))+
  xlab("Time (Months)")+
  ggtitle("Weibull 1 change-point model (Common Hazards Post Change-point)")+
  geom_line(data = Surv_extrap,
            aes(x = time, y = Survival,color = Treatment ))+
  geom_point(data = data.frame(time = mean.cp,
                               Survival = mean(pull(Surv_extrap,Survival)[which(dev == min(dev))])),
             aes(x = time, Survival), shape = 23, fill = "green",
             color = "darkred", size = 4, inherit.aes = F, stroke  = 1.5)+
  theme_bw()
ggsave(paste0(pathway, "Tretament-Delay-TA347.png"), width = 10, height = 5, units = 'in')



if(FALSE){
  
library("survHE",lib.loc = "~/R-packages/" )
fit.hmc <- survHE::fit.models(Surv(time,status)~as.factor(arm), data = TA347_df_all,method = "hmc",  dist = "weibullPH")
plot(fit.hmc, add.km = T, t = seq(0,50, by = 0.01))
beta_eval <- rstan::extract(fit.hmc$models$`Weibull (PH)`)$beta 
data.stan <- fit.hmc$misc$data.stan[[1]]
linpred <- beta_eval %*% t(data.stan$X)

shape <- rstan::extract(fit.hmc$models$`Weibull (PH)`)$alpha 
logf <- matrix(unlist(lapply(1:nrow(linpred), function(i) {
  data.stan$d * log(hweibullPH(data.stan$t, shape =  shape[i],scale=  exp(linpred[i, 
  ]))) + log(1 - pweibullPH(data.stan$t,shape =  shape[i], scale = exp(linpred[i,])))})), nrow = nrow(linpred), byrow = T)


p.hmc = make.surv(fit.hmc,nsim=1000, t = seq(0, 50, seq = 0.1))

png(paste0(pathway, "Tretament-Delay-TA347-Weibull-No-CP.png"), width = 10, height = 5, units = 'in', res = 300)
plot(TA347.km, xlim = c(0, max(data_jags$t_pred)), ylab = "Survival", xlab = "Time (Months)",
     sub=paste0("WAIC - ", round(loo::waic(logf)[["estimates"]][3, 1], digits = 2)))
lines(x =p.hmc$S[[1]]$t , y =p.hmc$S[[1]]$S, col = "blue")
lines(x =p.hmc$S[[2]]$t , y =p.hmc$S[[2]]$S, col = "red")
legend("topright", legend=c("Nintedanib + Docetaxel", "Docetaxel"),
       col=c("blue", "red"), lty=1, cex=0.8)
#points(x = TA347_df_all$time[TA347_df_all$arm == 1], y = St_mat_quantile[2,TA347_df_all$arm == 1], col = "red")
#points(x = TA347_df_all$time[TA347_df_all$arm == 0], y = St_mat_quantile[2,TA347_df_all$arm == 0], col = "blue")
dev.off()


}
fs3 <- flexsurvreg(Surv(time, status) ~ factor(arm), data = TA347_df_all,
                   dist = "weibullPH")


s <- summary(TA347.km, times = TA347.km$time, extend = T)

plot_data <- tibble(
  'time' = s$time,
  'n.risk' = s$n.risk,
  'n.event' = s$n.event,
  'n.censor' = s$n.censor,
  'Survival' = s$surv,
  'std.error' = s$std.err,
  'Treatment' = s$strata
)
plot_data <- plot_data %>% mutate(Treatment = ifelse(Treatment == "arm=1", "Nintedanib + Docetaxel", "Docetaxel")) 

Surv_extrap <- data.frame(time = rep(t_pred,2),
                          Survival =c(summary(fs3, t= t_pred)[[1]]$est,summary(fs3, t= t_pred)[[2]]$est), 
                          Treatment= c(rep("Nintedanib + Docetaxel",length(t_pred)),
                                       rep("Docetaxel",length(t_pred))))


p <- ggplot(plot_data, aes(x = time, y = Survival, color = Treatment))+
  geom_step( size = 1)+
  ylim(c(0,1))+
  xlim(c(0,50))+
  xlab("Time (Months)")+
  ggtitle("Standard Weibull Survival Model")+
  geom_line(data = Surv_extrap,
            aes(x = time, y = Survival,color = Treatment ))+
  theme_bw()
ggsave(paste0(pathway, "Tretament-Delay-TA347-Weibull-No-CP.png"), width = 10, height = 5, units = 'in')




fit.coxph <- coxph(Surv(time,status)~as.factor(arm), data = TA347_df_all)

cox.zph(fit.coxph)

if(FALSE){
  plot(fit.hmc, t = 0:50 ,xlim = c(0,50), ci = F)
  plot(fit.mle, t = 0:50 ,xlim = c(0,50), type = "haz")
  
}



fs3 <- flexsurvreg(Surv(time, status) ~ arm + sigma(arm), data = TA347_df_all,
                   dist = "gengamma")
flexsurvreg(Surv(time, status) ~ factor(arm), data = TA347_df_all,
            dist = "weibullPH")

fs3 <- flexsurvspline(Surv(time, status) ~ factor(arm) , data = TA347_df_all,
                   k = 1)

TA347_df_all_test <- TA347_df_all
TA347_df_all_test$arm <- as.factor(TA347_df_all_test$arm)
fs3 <- flexsurvspline(Surv(time, status) ~ arm , data = TA347_df_all_test,
                      k = 1, anc = list(gamma1 = ~ arm))

plot(fs3, t = 0:50 ,xlim = c(0,50), ci = F)
plot(fs3, t = 0:50 ,xlim = c(0,50), type = "haz")


#Bagust and Beale hazard convergence.


path_BRIM <- paste0(pathway,"BRIM-3/BRIM-3/")

TA269_OS_D_df <- read.xlsx(paste0(path_BRIM,"Dacarbazine","/Pseudo_IPD.xlsx"), header = T, sheetIndex = 1)[,-1]
TA269_OS_D_df$arm <- 0

TA269_OS_V_df <- read.xlsx(paste0(path_BRIM,"Vemurafenib","/Pseudo_IPD.xlsx"), header = T, sheetIndex = 1)[,-1]
TA269_OS_V_df$arm <- 1

# TA269_PFS_D_df <- read.table(paste0(path_mod_eval,"TA269_D_PFS_Initial","/IPDdata.txt"), header = T)
# TA269_PFS_D_df$arm <- 0
# TA269_PFS_V_df <- read.table(paste0(path_mod_eval,"TA269_V_PFS_Initial","/IPDdata.txt"), header = T)
# TA269_PFS_V_df$arm <- 1
# TA269_OS_D_df <- read.table(paste0(path_mod_eval,"TA269_D_OS_Initial","/IPDdata.txt"), header = T)
# TA269_OS_D_df$arm <- 0
# TA269_OS_V_df <- read.table(paste0(path_mod_eval,"TA269_V_OS_Initial","/IPDdata.txt"), header = T)
# TA269_OS_V_df$arm <- 1

#TA269_df_all <- rbind(TA269_PFS_D_df,TA269_PFS_V_df)
TA269_df_all <- rbind(TA269_OS_V_df,TA269_OS_D_df)
TA269.km <- survfit(Surv(time,status)~as.factor(arm), data = TA269_df_all)
plot(TA269.km)
#survminer::ggsurvplot(TA269.km,conf.int = TRUE)
#TA269.km$lower

jags.piecewise.wei_chng_common2 <-"

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
    sd[i]  <- 100 # ~ dunif(0.2,5)
    prec[i] <- pow(sd[i], -2)
  }
  
  for(k in 1:(N_CP+1)){ #Independent up to the last interval
    for(j in 1:2){ #Always the number of basic parameters plus 1 for a covariate 
      beta_cp[k,j] ~ dnorm(0,prec[1])
      beta_cp_anc[k,j] ~ dnorm(0,prec[1]) # Weibull or Exponential model (set to 0)
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
  for (i in 1:N1) {
    
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
      param_1[i,k] <- exp(linpred_cp[i,k] + linpred_fixed[i]) 
      param_2[i,k] <- exp(beta_cp_anc[k,1]*equals(scenario_expo,0)) #We allow an independent model for the first time-point 
      
      log_haz_seg[i,k] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))*X_ind[i,k] 
      cum_haz_seg[i,k] <- param_1[i,k]*pow(X[i,k]+cp[k],param_2[i,k]) -  param_1[i,k]*pow(cp[k],param_2[i,k])
      
    }
    
    log_haz[i] <- sum(log_haz_seg[i,])
    cum_haz[i] <- sum(cum_haz_seg[i,])
    
    loglik[i] = status[i]*log_haz[i] - cum_haz[i]
    zero[i] ~ dpois(C - loglik[i])
    
  }
  
#Other arm  
for (i in (N1+1):(N1+N2)) {
  
  for(k in (N_CP+1):(N_CP+1)){
   
      for(p in 1:ncovar_cp){
       linpred_cp_param[i,k,p] <- X_mat_trt[i,p]*beta_cp[k,p]
      }
      
  param_1[i,k] <- exp(beta_cp[k,1]) 
  param_2[i,k] <- exp(beta_cp_anc[k,1]*equals(scenario_expo,0)) #We allow an independent model for the first time-point 
      
  log_haz[i] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))
  cum_haz[i] <- param_1[i,k]*pow(time[i],param_2[i,k]) 
    
    }   
    
    loglik[i] = status[i]*log_haz[i] - cum_haz[i]
    zero[i] ~ dpois(C - loglik[i])
    
  }
  

  total_loglik <- sum(loglik)
}

"



data_jags <- list()
data_jags$N2 <- sum(TA269_df_all$arm == 0)
data_jags$N1 <- sum(TA269_df_all$arm == 1)
data_jags$N <- data_jags$N2+data_jags$N1
data_jags$time <- TA269_df_all$time
data_jags$status <- TA269_df_all$status
data_jags$MAXX <- max(TA269_df_all$time)
data_jags$N_CP <- 1
data_jags$ncovar <- 2
data_jags$index_pred <- c(min(which(TA269_df_all$arm == 1)),min(which(TA269_df_all$arm == 0))) #needs to be a vector 


data_jags$ncovar_fixed <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix
data_jags$ncovar_cp <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix
data_jags$t_pred <- seq(0,50, by = 0.01)
X_mat_trt <- matrix(c(rep(1, data_jags$N), rep(0, data_jags$N)),  ncol = 2)


data_jags$X_mat_trt <-X_mat_trt
data_jags$X_mat_fixed <- matrix(c(rep(0, (data_jags$N)*2)),  ncol = 2)

data_jags$scenario_expo <- 0


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




mod.cp1 <- runjags::run.jags(
  model = get("jags.piecewise.wei_chng_common2"),
  #model = get(paste0("jags.piecewise.wei_chng_common_haz")),
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp","beta_cp_ancs" ,"prec","sd",  "cum_haz", "loglik","total_loglik","beta_covar_fixed", "beta_cp_anc"), 
  sample=n.samp, 
  thin = n.thin, 
  burnin = n.burnin,
  inits = inits, 
  method ='rjparallel')

data_jags$scenario_expo <-1
  
mod.cp2 <- runjags::run.jags(
  model = get("jags.piecewise.wei_chng_common2"),
  #model = get(paste0("jags.piecewise.wei_chng_common_haz")),
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp", "prec","sd",  "cum_haz", "loglik","total_loglik","beta_covar_fixed", "beta_cp_anc"), 
  sample=n.samp, 
  thin = n.thin, 
  burnin = n.burnin,
  inits = inits, 
  method ='rjparallel')

mod_vec <- c("Weibull", "Exponential")


for(mod in 1:2){
  assign("mod.cp", get(paste0("mod.cp",mod)))
  param_names <- colnames(mod.cp$mcmc[[1]])
  mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])
  #St_mat <- exp(-mcmc_output[,grep("cum_haz_pred[",param_names, fixed = T)])
  St_mat <- exp(-mcmc_output[,grep("cum_haz[",param_names, fixed = T)])
  St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))
  
  loglik_mat <- mcmc_output[,grep("loglik",param_names)]
  loglik_mat <- loglik_mat[,-grep("total_loglik",colnames(loglik_mat))]
  
  waic <- loo::waic(loglik_mat)[["estimates"]][3, 1]
  pml <- -2 * sum(log(nrow(loglik_mat)/colSums(1/exp(loglik_mat))))
  
  
  mean_cp <- mean(mcmc_output[,grep("cp[2]",param_names,fixed = T)])
  
  
  # add.summary(mod.cp, "cp")
  # add.summary(mod.cp, "total_loglik")
  
  
  
  data_jags$t_pred <- seq(0, 50, by = 0.2)
  St_pred2 <- St_pred <- matrix(ncol = length(data_jags$t_pred), nrow = nrow(mcmc_output))
  #Fix plot
  t_pred <- data_jags$t_pred
  for(i in 1:nrow(mcmc_output)){
    
    m1 <- exp(mcmc_output[i,"beta_cp[1,1]"])
    m2 <- exp(mcmc_output[i,"beta_cp[2,1]"])
    if(mod_vec[mod] == "Exponential"){
      a1 <- 1
      a2 <- 1
    }else{
      a1 <- exp(mcmc_output[i,"beta_cp_anc[1,1]"])
      a2 <- exp(mcmc_output[i,"beta_cp_anc[2,1]"])
      
    }
    
    St_pred_cp <- flexsurv::pweibullPH(mcmc_output[i,"cp[2]"], scale = m1, shape = a1, lower.tail = F, log = F)
    
    
    for(j in 1:length(t_pred)){
      if(t_pred[j] <= mcmc_output[i,"cp[2]"]){
        St_pred[i,j] <- flexsurv::pweibullPH(t_pred[j], scale = m1, shape = a1, lower.tail = F, log = F)
        St_pred2[i,j] <- flexsurv::pweibullPH(t_pred[j], scale = m2, shape = a2, lower.tail = F, log = F)
      }else{
        
        X_eval <- t_pred[j] -mcmc_output[i,"cp[2]"]
        
        t_0 <- mcmc_output[i,"cp[2]"]
        t_1 <- t_pred[j]
        
        cum_haz_seq <-  flexsurv:::HweibullPH(t_1, scale = m2, shape = a2) -flexsurv:::HweibullPH(t_0, scale = m2, shape = a2)
        St_pred[i,j] <- exp(-(-log(St_pred_cp) + cum_haz_seq))
        St_pred2[i,j] <- flexsurv::pweibullPH(t_pred[j], scale = m2, shape = a2, lower.tail = F, log = F)
        
        
      }
    }
    
    
  }
  
  index_trt_1 <- 1:TA269.km$strata[1]
  index_trt_2 <-(TA269.km$strata[1]+1):(TA269.km$strata[1]+TA269.km$strata[2])
  
  s <- summary(TA269.km, times = TA269.km$time, extend = T)
  
  plot_data <- tibble(
    'time' = s$time,
    'n.risk' = s$n.risk,
    'n.event' = s$n.event,
    'n.censor' = s$n.censor,
    'Survival' = s$surv,
    'std.error' = s$std.err,
    'Treatment' = s$strata,
    'lower' = s$lower,
    'upper' = s$upper,
  )
  plot_data <- plot_data %>% mutate(Treatment = ifelse(Treatment == "as.factor(arm)=0", "Dacarbazine", "Vemurafenib")) 
  
  Surv_extrap <- data.frame(time = rep(t_pred,2),
                            Survival =c(colMeans(St_pred),
                                        colMeans(St_pred2)), 
                            Treatment= c(rep("Vemurafenib",length(t_pred)),
                                         rep("Dacarbazine",length(t_pred))))
  mean.cp <- mean(mcmc_output[,"cp[2]"])
  dev <-abs(Surv_extrap$time - mean.cp)
  

  
  
  
  
  png(paste0(pathway, mod_vec[mod],"-Haz-Converge-TA269.png"), width = 10, height = 5, units = 'in', res = 300)
  plot(TA269.km, xlim = c(0, max(data_jags$t_pred)), xlab = "Time (Months)", ylab = "S(t)", 
       main = paste0(mod_vec[mod]," 1 change-point model (Converge Haz hazards)"),
       sub=paste0("WAIC - ", round(waic, digits = 2),": Mean change-point - ", round(mean_cp, digits= 2)))
  lines(x = TA269.km$time[index_trt_1], y = TA269.km$lower[index_trt_1], lty= 2)
  lines(x = TA269.km$time[index_trt_1], y = TA269.km$upper[index_trt_1], lty= 2)
  lines(x = TA269.km$time[index_trt_2], y = TA269.km$lower[index_trt_2], lty= 2)
  lines(x = TA269.km$time[index_trt_2], y = TA269.km$upper[index_trt_2], lty= 2)
  #lines(x =data_jags$t_pred , y =St_mat_quantile[2,1:(ncol(St_mat_quantile)/2)], col = "blue")
  #lines(x =data_jags$t_pred , y =St_mat_quantile[2,(ncol(St_mat_quantile)/2 +1 ):ncol(St_mat_quantile)], col = "red")
  #points(x = TA269_df_all$time[TA269_df_all$arm == 1], y = St_mat_quantile[2,TA269_df_all$arm == 1], col = "red")
  #points(x = TA269_df_all$time[TA269_df_all$arm == 0], y = St_mat_quantile[2,TA269_df_all$arm == 0], col = "blue")
  lines(x = data_jags$t_pred, y = colMeans(St_pred2), col = "blue")
  lines(x = data_jags$t_pred, y = colMeans(St_pred), col = "red")
  dev.off()
  #Independent change-points
  
  
  p <- ggplot(plot_data, aes(x = time, y = Survival, color = Treatment))+
    geom_step( size = 1)+
    ylim(c(0,1))+
    xlim(c(0,50))+
    xlab("Time (Months)")+
    ggtitle(paste0(mod_vec[mod]," 1 change-point model (Converge Haz hazards)"))+
    geom_line(data = Surv_extrap,
              aes(x = time, y = Survival,color = Treatment))+
    geom_line(data = plot_data,
              aes(x = time, y = lower,color = Treatment), linetype = "dashed")+
    geom_line(data = plot_data,
              aes(x = time, y = upper,color = Treatment), linetype = "dashed")+
    geom_point(data = data.frame(time = mean.cp,
                                 Survival = mean(pull(Surv_extrap,Survival)[which(dev == min(dev))[1]])),
               aes(x = time, Survival), shape = 23, fill = "green",
               color = "darkred", size = 4, inherit.aes = F, stroke  = 1.5)+
    theme_bw()
  
  
  ggsave(paste0(pathway, mod_vec[mod],"-Haz-Converge-TA269.png"), width = 10, height = 5, units = 'in')
}



#Cox Change-point Model


library("nimble", lib.loc =  "~/R-packages")
#https://stats.stackexchange.com/questions/15279/change-point-in-cox-survival-model

head(E1690.dat)

E1690.dat_mod <- E1690.dat[, c("SURVTIME","SURVCENS", "TRT", "AGE_scale")]
LL_cox_ph <-  function(params_beta,params_cp,  data){
  library("survival")
  params_cp <- head(params_cp, n= -1)
  data <- as.data.frame(data)
  colnames(data) <- c("time","status","trt", "age") 
  data$trt <- as.factor(data$trt)
  
  data_split <- survSplit(data,cut=params_cp,
                          end="time",start="start",
                          event="status", episode="period")
  
  cox_fit <- coxph(Surv(start, time, status) ~ 
                     I((trt==2)&(period==1)) + I((trt==2)&(period==2)) + age, 
                   data=data_split,
                   init =params_beta, #
                   control = coxph.control(iter.max = 0))
  #model.matrix(cox_fit)
  
  return(cox_fit$loglik[2])
  
}

R_LL_cox_ph <- nimbleRcall(function(params_beta = double(1),
                                    params_cp = double(1),
                                    data = double(2)){},
                           Rfun = 'LL_cox_ph',
                           returnType = double(0))


mod <- nimbleCode({
  
  for(i in 1:n_param_covar){
    params_beta[i] ~ dunif(-10,10)
    # params_beta[i] <- 1
    # params_cp[i] <- 10
  }
    
  
  for(i in 1:n_param){
   params_cp[i] ~ dunif(0,max_cp)
     # params_beta[i] <- 1
    # params_cp[i] <- 10
  }
  
  C <- 1000
  
  LL <- R_LL_cox_ph(params_beta[1:n_param_covar],
                    params_cp[1:n_param],
                    data[1:df_nrow,1:df_ncol])
  
  zero ~ dpois(-LL + C)
  
  
})

data_mod <- as.matrix(E1690.dat_mod)
#data_mod <- E1690.dat_mod

#data_mod$x = as.numeric(as.factor(data_mod$x))
params_beta = c(1,1) #coefficient for hazard(s) and changepoint
params_cp = c(10, -999) #change-point padding the last val

#data_mod <- matrix(c(data_mod$SURVTIME, data_mod$SURVCENS, data_mod$TRT), nrow = nrow(data_mod), ncol = ncol(data_mod))

####

#LL_cox_ph(params_beta = c(1,1), params_cp = c(5, -999), data = data_mod )
data_nimble = list(zero = 0, data =data_mod, max_cp = max(data_mod[,1])*0.7)
n_cp <- 1
constants_nimble = list(n_cp = n_cp,   
                        n_param = n_cp +1,
                        df_nrow = nrow(data_mod),
                        df_ncol = ncol(data_mod),
                        n_param_covar = 3)


undebug(LL_cox_ph)


#runMCMC(mod, niter = 10, nchains = 3,)
#https://oliviergimenez.github.io/nimble-workshop/#25
model <- nimbleMCMC(mod, 
                    data = data_nimble,
                    constants = constants_nimble,
                    monitors = c("params_beta", "params_cp", "LL"),
                    niter = n.samp*10, 
                    nchains =2, 
                    nburnin = n.burnin*10)
model_output_compare2 <- model_res_final <- rbind(model[[1]], model[[2]])

plot(density(model_res_final[, "params_cp[1]"]))
plot(density(exp(model_res_final[, "params_beta[1]"])),xlim = c(0,2))
plot(density(exp(model_res_final[, "params_beta[2]"])), xlim = c(0,2))
plot(density(exp(model_res_final[, "params_beta[3]"])), xlim = c(0,2))

df_params <- rbind(quantile(model_res_final[, "params_cp[1]"], probs = c(0.025,0.5,0.975)),
      quantile(exp(model_res_final[, "params_beta[1]"]), probs = c(0.025,0.5,0.975)),
      quantile(exp(model_res_final[, "params_beta[2]"]), probs = c(0.025,0.5,0.975)),
      quantile(exp(model_res_final[, "params_beta[3]"]), probs = c(0.025,0.5,0.975)))

df_params_final <- cbind(Param = c("Change-point", "HR 1", "HR 2", "HR Age"),
                as.data.frame(df_params))

write.csv(df_params_final,paste0(pathway, "summary_output_df_cox_partial.csv"))

plot(density(model_res_final[, "LL"]))
plot(density(model_res_final[, "params_cp[1]"]))
hist(model_res_final[, "params_cp[1]"])

res.df <- data.frame(changepoint = c(model_output_compare[,"cp[2]"],model_output_compare2[, "params_cp[1]"]),
           Model = rep(c("Weibull", "Cox"), c(nrow(model_output_compare),nrow(model_output_compare2))))


mu <- res.df %>% group_by(Model) %>% summarise(grp.mean=mean(changepoint))

ggplot(res.df, aes(x = changepoint,fill = Model))+
  geom_density(alpha=0.4)+
  xlab("Time (Years)")+
  ggtitle("Posterior Distribution of Change-point")+
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
  #           linetype="dashed")+
  theme_bw()
  
           
ggsave(paste0(pathway, "Changepoint-Posterior-Distribution.png"), width = 10, height = 5, units = 'in')

                       
                       
# Waning an applied example.

path_KEYNOTE <- paste0(pathway,"/KEYNOTE-24/KEYNOTE-024/")

TA531_OS_Pembro_df <- read.xlsx(paste0(path_KEYNOTE,"Pembro","/Pseudo_IPD.xlsx"), header = T, sheetIndex = 1)[,-1]
TA531_OS_Pembro_df$arm <- 1#"Pembro"

TA531_OS_Chemo_df <- read.xlsx(paste0(path_KEYNOTE,"Chemo","/Pseudo_IPD.xlsx"), header = T, sheetIndex = 1)[,-1]
TA531_OS_Chemo_df$arm <- 0#"Chemo"

TA531_OS_all_df <- rbind(TA531_OS_Pembro_df, TA531_OS_Chemo_df)
TA531_OS_km <- survfit(Surv(time,status)~as.factor(arm), data = TA531_OS_all_df )




# data_jags <- list()
# data_jags$N <- nrow(TA531_OS_all_df)
# data_jags$time <- TA531_OS_all_df$time
# data_jags$status <- TA531_OS_all_df$status
# data_jags$MAXX <- max(TA531_OS_all_df$time)
# data_jags$N_CP <- 1
# data_jags$ncovar <- 2
# data_jags$index_pred <- c(1,3) #needs to be a vector 
# 
# 
# #data_jags$N_CP <-1
# data_jags$ncovar_fixed <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix
# data_jags$ncovar_cp <- 2 # will always be at least 2 even if thery are none.. just add a dummy matrix
# 
# data_jags$t_pred <- seq(0,100, by = 0.2)
# 
# #X_mat_fixed is just a dummy matrix
# #PH Model with Changepoint
# 
# data_jags$X_mat_fixed <- matrix(rep(0, data_jags$N*2),  ncol = 2)
# data_jags$X_mat_trt <- matrix(c(rep(0, data_jags$N), TA531_OS_all_df$arm),  ncol = 2)
# #independent parameters 
# 
# 
# #Change-inits
# inits <- function(){
#   list(unif = runif(data_jags$N_CP,max = data_jags$MAXX/data_jags$N_CP)
#        #,
#        #beta_cp_0 = c(0),
#        #beta_cp = rnorm(data_jags$N_CP+1)
#   )
# }
# 
# #data_jags$t_pred <- 0:20
# mod <- "wei"
# mod.cp <- runjags::run.jags(
#   model = get(paste0("jags.piecewise.",mod,"_chng_wane")),
#   #model = jags.piecewise.expo_wane,
#   data = data_jags,
#   n.chains = 2,
#   monitor = c("cp", "beta_cp", "total_loglik", "cum_haz", "lower_int", "upper_int", "beta_cp_anc", "loglik"), 
#   sample=n.samp, 
#   thin = n.thin, 
#   burnin = n.burnin,
#   inits = inits, 
#   method ='rjparallel')
# 
# add.summary(mod.cp, "cp")
# param_names <- colnames(mod.cp$mcmc[[1]]) #rownames(expo.mod_1chng$BUGSoutput$summary)
# mcmc_output <- rbind(mod.cp$mcmc[[1]],mod.cp$mcmc[[2]])
# 
# 
# loglik_mat <- mcmc_output[,grep("loglik",param_names)]
# loglik_mat <- loglik_mat[,-grep("total_loglik",colnames(loglik_mat))]
# 
# waic <- loo::waic(loglik_mat)[["estimates"]][3, 1]
# pml <- -2 * sum(log(nrow(loglik_mat)/colSums(1/exp(loglik_mat))))
# 
# St_mat <- exp(-mcmc_output[,grep("cum_haz[",param_names, fixed = T)])
# St_mat_quantile <-apply(St_mat, 2, quantile, probs = c(0.025,0.5,0.975))
# 
# plot(TA531_OS_km, xlim = c(0, max(data_jags$t_pred)))
# points(x = TA531_OS_all_df$time[TA531_OS_all_df$arm == 0], y = St_mat_quantile[2,TA531_OS_all_df$arm == 0], col = "red")
# points(x = TA531_OS_all_df$time[TA531_OS_all_df$arm == 1], y = St_mat_quantile[2,TA531_OS_all_df$arm == 1], col = "blue")
# 
# t_pred <- seq(0, 50, by = 0.1)
# HR_mat <- St_pred2 <- St_pred <- matrix(ncol = length(t_pred), nrow = nrow(mcmc_output))
# 
# if(weib){
#   incomplete_gamma_upper<- function(s,x){
#     exp(lgamma(s))*(1-pgamma(x,s,1))
#   }
#   #Write it out in integral calc
#   cum_haz_seg <- function(t_0,t_1,a, m, initial_HR, lambda_wane  ){
#     constant <- (initial_HR-1)*exp(t_0*lambda_wane)
#     int_1 <- (t_0^a)/a-(incomplete_gamma_upper(a,lambda_wane*t_0)*(t_0^a)*constant)/((lambda_wane*t_0)^a)
#     int_2 <- (t_1^a)/a-(incomplete_gamma_upper(a,lambda_wane*t_1)*(t_1^a)*constant)/((lambda_wane*t_1)^a)
#     return((int_2-int_1)*a*m)
#   }
#   
#   
# }else{
#   cum_haz_seg <- function(t_0,t_1, m, initial_HR, lambda_wane , a = NULL ){
#     lower_int <- m*(t_0-((initial_HR-1)*exp(-lambda_wane*(t_0-t_0))/lambda_wane))
#     upper_int <- m*(t_1-((initial_HR-1)*exp(-lambda_wane*(t_1-t_0))/lambda_wane))
#     #lambda*(t-((r-1)*e^(-lambda_2*(t-q)))/lambda_2) same a int cal
#     #t0
#     #Not equivalent to setting t_0 equal to zero!!
#     return((upper_int-lower_int))
#   }
#   
# }
# 
# #Fix plot
# for(i in 1:nrow(mcmc_output)){
#   
#   m1_1 <- exp(mcmc_output[i,"beta_cp[1,1]"] + mcmc_output[i,"beta_cp[1,2]"] )
#   m1_2 <- exp(mcmc_output[i,"beta_cp[1,1]"])
#   
#   if(weib == FALSE){
#     a1 <- 1
#   }else{
#     a1 <- exp(mcmc_output[i,"beta_cp_anc[1,1]"])
#   }
#   
#   St_pred_cp <- flexsurv::pweibullPH(mcmc_output[i,"cp[2]"], scale = m1_1, shape = a1, lower.tail = F, log = F)
#   St_pred2_cp <- flexsurv::pweibullPH(mcmc_output[i,"cp[2]"], scale = m1_2, shape = a1, lower.tail = F, log = F)
#   
#   
#   for(j in 1:length(t_pred)){
#     if(t_pred[j] <= mcmc_output[i,"cp[2]"]){
#       St_pred[i,j] <- flexsurv::pweibullPH(t_pred[j], scale = m1_1, shape = a1, lower.tail = F, log = F)
#       St_pred2[i,j] <- flexsurv::pweibullPH(t_pred[j], scale = m1_2, shape = a1, lower.tail = F, log = F)
#       HR_mat[i,j] <-  exp(mcmc_output[i,"beta_cp[1,2]"])
#     }else{
#       
#       m2 <- exp(mcmc_output[i,"beta_cp[2,1]"])
#       
#       if(weib == FALSE){
#         a2 <- 1
#       }else{
#         a2 <- exp(mcmc_output[i,"beta_cp_anc[2,1]"])
#       }
#       
#       
#       X_eval <- t_pred[j] -mcmc_output[i,"cp[2]"]
#       
#       initial_HR <-  exp(mcmc_output[i,"beta_cp[1,2]"])
#       lambda_wane <- exp(mcmc_output[i,"beta_cp[2,2]"])
#       HR_final <-1-(1-initial_HR)*exp(-lambda_wane*X_eval)
#       HR_mat[i,j] <-  HR_final
#       
#       t_0 <- mcmc_output[i,"cp[2]"]
#       t_1 <- t_pred[j]
#       
#       St_pred[i,j] <- exp(-(-log(St_pred_cp) +cum_haz_seg(t_0 = t_0, t_1 = t_1,a = a2, m = m2, initial_HR = initial_HR, lambda_wane = lambda_wane)))
#       St_pred2[i,j] <- exp(-(-log(St_pred2_cp) +flexsurv:::HweibullPH(t_1, scale = m2, shape = a2)-flexsurv:::HweibullPH(t_0, scale = m2, shape = a2)))
#       
#       
#       
#     }
#   }
#   
#   
# }
# mean_cp <- mean(mcmc_output[,grep("cp[2]",param_names,fixed = T)])
# 
# png(paste0(pathway, "Converge-Haz-Weibull-model-KEYNOTE-024.png"), width = 10, height = 5, units = 'in', res = 300)
# plot(TA531_OS_km, xlim = c(0, max(data_jags$t_pred)), xlab = "Time (Years)", ylab = "S(t)", 
#      main = "Weibull 1 change-point model (Converge Haz hazards)",
#      sub=paste0("WAIC - ", round(waic, digits = 2),": Mean change-point - ", round(mean_cp, digits= 2) ))
# lines(y = colMeans(St_pred), t_pred, col = "blue")
# lines(y = colMeans(St_pred2), t_pred, col = "red")
# dev.off()
# 
# quantile_HR <- apply(HR_mat, 2, quantile, probs = c(0.025, 0.5,0.975))
# png(paste0(pathway, "Converge-HR-Weibull-model-KEYNOTE-024.png"), width = 10, height = 5, units = 'in', res = 300)
# plot(x = 1:length(t_pred), y =quantile_HR[2,] ,xlim = c(0, 50), ylim = c(0,1.5), type = "l", 
#      xlab = "Time (Months)", ylab = "HR")
# lines(x = 1:length(t_pred), y =quantile_HR[1,],lty= 2)
# lines(x = 1:length(t_pred), y =quantile_HR[3,], lty = 2)
# dev.off()

weibull_ph <- read.csv(paste0(pathway,"summary_output_df_weib_cox.csv")) %>% mutate_if(is.numeric, round, digits = 2)
cox_ph <- read.csv(paste0(pathway,"summary_output_df_cox_partial.csv")) %>% mutate_if(is.numeric, round, digits = 2)
print(xtable::xtable(cox_ph[,-1]), include.rownames=FALSE)
print(xtable::xtable(weibull_ph[,-1]), include.rownames=FALSE)


