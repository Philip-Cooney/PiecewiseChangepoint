
path <- "C:/Users/phili/Desktop/Hazard-Changepoint-main/"
setwd(path)
source(paste0(path,"Functions Gibbs Sampler.R"))
source(paste0(path,"Functions Collapsing.R"))
library(loo)



#Stanford
library(survival)
stanford2 <- survival::stanford2
stanford2$time <- stanford2$time/365
surv.stanford <- survfit(Surv(time,status)~1,data = stanford2)
plot(surv.stanford)

dist.res  <- flexsurvreg(Surv(time,status)~1,data = stanford2, dist = "gengamma")

cum_haz_stan <- survminer::ggsurvplot(surv.stanford, risk.table = TRUE, fun = "cumhaz",break.time.by =1)
ggsave(file = "cum_haz_stan.png", print(cum_haz_stan))
ggsave("cum_haz_stan.png")



test <- collapsing.model(n.iter = 20750, df = stanford2,n.chains = 1,  burn_in = 750,
                         max_predict = 10,lambda.prior = 1, alpha.hyper =  1,beta.hyper1 = 1, beta.hyper2 = 1, obs_true = stanford2)

test[["output"]][["2"]]

stanford_test <-stanford2
for(i in 1:nrow(stanford_test)){
  if(stanford_test$time[i] >= 2){
    stanford_test$time[i] <- 2
    stanford_test$status[i] <-0
  }
}


# Returns the log likelihood of the model fitted to the observed data (df) vs the actual data obs_true
test <- collapsing.model(n.iter = 20750, df = stanford_test,n.chains = 1,  burn_in = 750,
                         max_predict = 10,lambda.prior = 1, alpha.hyper =  1,beta.hyper1 = 1, beta.hyper2 = 1, obs_true = stanford2, interval = 0.02)


test$plot_Sur_all

expo <- "model{
  for(i in 1:N){
    is.censored[i]~dinterval(t[i],t.cen[i])
    t[i] ~ dexp(lambda)
    Like[i] <- ifelse(is.censored[i], 1- pexp(t.cen[i],lambda), dexp(t[i], lambda))
    invLik[i] <- 1/Like[i]
    
  }
  for(i in 1:length(t_pred)){
    St_pred[i] <- 1- pexp(t_pred[i],lambda)
  }
  lambda ~ dgamma(1,1)
  
}"
df <- stanford_test[,c("time","status")]
data_new <- list()
df_jags <- df
df_jags$t <- df$time

tinits1<-df_jags$t + max(df$time)
is.na(tinits1)<-df_jags$status==1
tinits2<-tinits1 + 5

is.na(df_jags$t)<-df_jags$status==0
df_jags$is.censored<-1-df_jags$status
df_jags$t.cen<-df_jags$time+df_jags$status


modelinits <- list(list(t = tinits1),
                   list(t = tinits2))
data_jags <- list(N = nrow(df_jags),
                  t.cen = df_jags$t.cen,
                  is.censored = df_jags$is.censored,
                  t = df_jags$t)
data_jags$t_pred <- seq(0, 10, by = 0.02)


expo.jags<-R2jags::jags(model.file = textConnection(expo),
                        data=data_jags,
                        inits=modelinits,
                        n.chains=2,
                        n.iter = 10000,
                        parameters.to.save = c("Like", "lambda","invLik","St_pred"))

Like.sims.expo <- expo.jags$BUGSoutput[["sims.matrix"]][,grep("Like",colnames(expo.jags$BUGSoutput[["sims.matrix"]]))]
LL.dens.expo <- density(rowSums(log(Like.sims.expo)))
PML.expo <- 1/expo.jags$BUGSoutput[["summary"]][grep("invLik",rownames(expo.jags$BUGSoutput[["summary"]])),1]
#Same as approach PML.expo above
PML.expo <- nrow(Like.sims.expo)/colSums(1/Like.sims.expo)
Surv.expo <- expo.jags$BUGSoutput[["summary"]][grep("St_pred",rownames(expo.jags$BUGSoutput[["summary"]])),1]

# Weibull Model


weibull <- "model{
  for(i in 1:N){
is.censored[i]~dinterval(t[i],t.cen[i])
t[i] ~ dweib(v,lambda)
Like[i] <- ifelse(is.censored[i], 1- pweib(t.cen[i],v,lambda), dweib(t[i],v, lambda))
invLik[i] <- 1/Like[i]
    
  }
for(i in 1:length(t_pred)){
St_pred[i] <- 1- pweib(t_pred[i],v,lambda)
}
lambda ~ dgamma(1,1)
v ~ dgamma(1,1)
}"

weibull.jags<-R2jags::jags(model.file = textConnection(weibull),
                           data=data_jags,
                           inits=modelinits,
                           n.chains=2,
                           n.iter = 10000,
                           parameters.to.save = c("lambda","v", "Like", "invLik","St_pred"))

PML.weib <- 1/weibull.jags$BUGSoutput[["summary"]][grep("invLik",rownames(weibull.jags$BUGSoutput[["summary"]])),1]
Like.sims.weib <- weibull.jags$BUGSoutput[["sims.matrix"]][,grep("Like",colnames(weibull.jags$BUGSoutput[["sims.matrix"]]))]
LL.dens.weibull <- density(rowSums(log(Like.sims.weib)))
Surv.weibull <- weibull.jags$BUGSoutput[["summary"]][grep("St_pred",rownames(weibull.jags$BUGSoutput[["summary"]])),1]




gamma.jags <- "model{
  for(i in 1:N){
is.censored[i]~dinterval(t[i],t.cen[i])
t[i] ~ dgamma(shape,lambda)
Like[i] <- ifelse(is.censored[i], 1- pgamma(t.cen[i],shape,lambda), dgamma(t[i],shape, lambda))
invLik[i] <- 1/Like[i]
  }
 for(i in 1:length(t_pred)){
    St_pred[i] <- 1- pgamma(t_pred[i],shape,lambda)
  }

lambda ~ dgamma(1,1)
shape ~ dgamma(1,1)
}"

gamma.mod <-R2mod::jags(model.file = textConnection(gamma.jags),
                        data=data_jags,
                        inits=modelinits,
                        n.chains=2,
                        n.iter = 10000,
                        parameters.to.save = c("lambda","shape","Like", "invLik","St_pred"))


PML.gamma <- 1/gamma.mod$BUGSoutput[["summary"]][grep("invLik",rownames(gamma.mod$BUGSoutput[["summary"]])),1]
Like.sims.gamma <- gamma.mod$BUGSoutput[["sims.matrix"]][,grep("Like",colnames(gamma.mod$BUGSoutput[["sims.matrix"]]))]
LL.dens.gamma <- density(rowSums(log(Like.sims.gamma)))
Surv.gamma <- gamma.mod$BUGSoutput[["summary"]][grep("St_pred",rownames(gamma.mod$BUGSoutput[["summary"]])),1]


lnorm.jags <- "model{
  for(i in 1:N){
is.censored[i]~dinterval(t[i],t.cen[i])
t[i] ~ dlnorm(mu,tau)
Like[i] <- ifelse(is.censored[i], 1- plnorm(t.cen[i],mu,tau), dlnorm(t[i],mu, tau))
invLik[i] <- 1/Like[i]
  }
 for(i in 1:length(t_pred)){
    St_pred[i] <- 1- plnorm(t_pred[i],mu,tau)
  }
mu ~ dnorm(0,0.1)
sd ~ dunif(0,5)
tau <- pow(sd,-2)
}"

lnorm.mod <-R2jags::jags(model.file = textConnection(lnorm.jags),
                         data=data_jags,
                         inits=modelinits,
                         n.chains=2,
                         n.iter = 10000,
                         parameters.to.save = c("mu","sd", "Like", "invLik","St_pred"))



PML.lnorm <- 1/lnorm.mod$BUGSoutput[["summary"]][grep("invLik",rownames(lnorm.mod$BUGSoutput[["summary"]])),1]
Like.sims.lnorm <- lnorm.mod$BUGSoutput[["sims.matrix"]][,grep("Like",colnames(lnorm.mod$BUGSoutput[["sims.matrix"]]))]
LL.dens.lnorm <- density(rowSums(log(Like.sims.lnorm)))
Surv.lnorm <- lnorm.mod$BUGSoutput[["summary"]][grep("St_pred",rownames(lnorm.mod$BUGSoutput[["summary"]])),1]



llogis.jags <- "
model{
for(i in 1:N){
is.censored[i]~dinterval(t.log[i],t.cen.log[i])
t.log[i] ~ dlogis(mu,tau)
Like[i] <- ifelse(is.censored[i], 1/(1 + pow(exp(t.cen.log[i])/beta, alpha)), 
          (alpha/beta)*pow(exp(t.log[i])/beta, alpha-1)/pow(1 + pow(exp(t.log[i])/beta,alpha),2))
invLik[i] <- 1/Like[i]

}
 for(i in 1:length(t_pred)){
    St_pred[i] <- 1/(1 + pow(t_pred[i]/beta, alpha))
  }



mu ~ dnorm(0,0.1)
scale ~ dgamma(1,1)
tau <- pow(scale,-1) # Inverse of scale which is beta on the log-logistic dist
beta <- exp(mu)
alpha <- tau

}"

data_jags_llogis <- data_jags
data_jags_llogis$t.log <- log(data_jags$t)
data_jags_llogis$t.cen.log <- log(data_jags$t.cen)

llogis.mod <-R2jags::jags(model.file = textConnection(llogis.jags),
                          data=data_jags_llogis,
                          inits=modelinits,
                          n.chains=2,
                          n.iter = 10000,
                          parameters.to.save = c("alpha","beta","Like", "invLik","St_pred"))

PML.llogis <- 1/llogis.mod$BUGSoutput[["summary"]][grep("invLik",rownames(llogis.mod$BUGSoutput[["summary"]])),1]
Like.sims.llogis<- llogis.mod$BUGSoutput[["sims.matrix"]][,grep("Like",colnames(llogis.mod$BUGSoutput[["sims.matrix"]]))]
LL.dens.llogis <- density(rowSums(log(Like.sims.llogis)))
Surv.llogis <- llogis.mod$BUGSoutput[["summary"]][grep("St_pred",rownames(llogis.mod$BUGSoutput[["summary"]])),1]



gompertz.jags <- "
data{
for(i in 1:N){
zero[i] <- 0}
}

model{

C <- 10000
for(i in 1:N){

logHaz[i] <- (log(b)+ a*time[i])*status[i]  
logSurv[i] <- (-b/a)*(exp(a*time[i])-1)

LL[i] <- logHaz[i]+ logSurv[i]
Like[i] <- exp(LL[i])

invLik[i] <-pow(Like[i],-1)


zero[i] ~ dpois(zero.mean[i])
zero.mean[i] <- -logHaz[i]-logSurv[i] + C
}

for(i in 1:length(t_pred)){
   St_pred[i] <- exp((-b/a)*(exp(a*t_pred[i])-1))
  }


a ~ dnorm(0,0.01)
b ~ dnorm(0,0.01)
#b ~ dgamma(0,5)
}"

data_surv <- list()

data_surv$time <- df$time
data_surv$status <- df$status
data_surv$N <- nrow(df)
data_surv$t_pred <- data_jags$t_pred

model_gomp <- flexsurvreg(Surv(time,status)~1, data = stanford_test, dist = "gompertz")
gompertz.mod <-R2jags::jags(model.file = textConnection(gompertz.jags),
                            data=data_surv,
                            inits=list(list(a = 1, b = 2),
                                       list(a = 2, b = 1)),
                            n.chains=2,
                            n.iter = 10000,
                            parameters.to.save = c("a","b", "Like", "invLik","St_pred"))



PML.gomp <- 1/gompertz.mod$BUGSoutput[["summary"]][grep("invLik",rownames(gompertz.mod$BUGSoutput[["summary"]])),1]
Like.sims.gomp<- gompertz.mod$BUGSoutput[["sims.matrix"]][,grep("Like",colnames(gompertz.mod$BUGSoutput[["sims.matrix"]]))]
LL.dens.gomp <- density(rowSums(log(Like.sims.gomp)))

Surv.gomp <- gompertz.mod$BUGSoutput[["summary"]][grep("St_pred",rownames(gompertz.mod$BUGSoutput[["summary"]])),1]



gen.gamma.jags <- "model{
  for(i in 1:N){
is.censored[i]~dinterval(t[i],t.cen[i])
t[i] ~ dgen.gamma(r,lambda,b)
Like[i] <- ifelse(is.censored[i], 1- pgen.gamma(t.cen[i],r,lambda,b), dgen.gamma(t[i],r,lambda,b))
invLik[i] <- 1/Like[i]

  }

 for(i in 1:length(t_pred)){
    St_pred[i] <- 1- pgen.gamma(t_pred[i],r,lambda,b)
  }
r ~ dgamma(1,1)
lambda ~ dgamma(1,1)
b ~ dgamma(1,1)
}"

gen.gamma.mod <-R2jags::jags(model.file = textConnection(gen.gamma.jags),
                             data=data_jags,
                             inits=modelinits,
                             n.chains=2,
                             n.iter = 10000,
                             parameters.to.save = c("r","lambda","b", "Like", "invLik","St_pred"))



PML.gen.gamma <- 1/gen.gamma.mod$BUGSoutput[["summary"]][grep("invLik",rownames(gen.gamma.mod$BUGSoutput[["summary"]])),1]
Like.sims.gen.gamma <- gen.gamma.mod$BUGSoutput[["sims.matrix"]][,grep("Like",colnames(gen.gamma.mod$BUGSoutput[["sims.matrix"]]))]
LL.dens.gen.gamma <- density(rowSums(log(Like.sims.gen.gamma)))
Surv.gen.gamma <- gen.gamma.mod$BUGSoutput[["summary"]][grep("St_pred",rownames(gen.gamma.mod$BUGSoutput[["summary"]])),1]



## Get AUC
mean.Surv_df_all <- test$mean.Surv_df_all
t <- mean.Surv_df_all$time

true_km <- survfit(Surv(time,status)~1,data = stanford2)
true_km.df <- data.frame(cbind(true_km[[c("time")]],true_km[[c("surv")]], true_km[[c("upper")]],true_km[[c("lower")]]))
colnames(true_km.df) <- c("time", "survival", "upper", "lower")

surv_km_time <- rep(NA, length(t))
for(i in 1:length(t)){
  if(t[i] < true_km.df$time[1]){
    surv_km_time[i] <- 1
  }else if(t[i] > true_km.df$time[nrow(true_km.df)]){
    surv_km_time[i] <- 0
  }else{
    surv_km_time[i] <- true_km.df$survival[max(which(true_km.df$time <= t[i]))]
  }
}

AUC_true <- integrate.xy(t, surv_km_time)
t2 <- t[which(t <= true_km.df$time[nrow(true_km.df)])]
t_pred <- data_jags$t_pred

AUC = c(AUC_true, integrate.xy(t2, mean.Surv_df_all$survival[1:length(t2)]),
        integrate.xy(t_pred, Surv.expo),
        integrate.xy(t_pred, Surv.weibull),
        integrate.xy(t_pred, Surv.llogis),
        integrate.xy(t_pred, Surv.lnorm),
        integrate.xy(t_pred, Surv.gomp),
        integrate.xy(t_pred, Surv.gen.gamma))

AUC_diff = c(NA, integrate.xy(t2, abs(mean.Surv_df_all$survival[1:length(t2)]-surv_km_time[1:length(t2)])),
             integrate.xy(t_pred, abs(Surv.expo - surv_km_time[1:length(t_pred)])),
             integrate.xy(t_pred, abs(Surv.weibull - surv_km_time[1:length(t_pred)])),
             integrate.xy(t_pred, abs(Surv.llogis - surv_km_time[1:length(t_pred)])),
             integrate.xy(t_pred, abs(Surv.lnorm - surv_km_time[1:length(t_pred)])),
             integrate.xy(t_pred, abs(Surv.gomp - surv_km_time[1:length(t_pred)])),
             integrate.xy(t_pred, abs(Surv.gen.gamma - surv_km_time[1:length(t_pred)]))) 




df$enter <- 0
df <- df[order(df$status),]

changepoint <- test$changepoint[, -ncol(test$changepoint)]
lambda <- test$lambda
#Collapsing.Model$prob.changepoint

# separate out the changepoint models



ind.expo <- function(time,status, lambda){
  if(status ==0){
    return(pexp(time,rate = lambda, lower.tail = F, log.p = TRUE ))
  }else{
    return(dexp(time,rate = lambda, log = TRUE ))
  }
}



changepoint.num <- apply(changepoint,1, function(x){length(na.omit(x))})
df.plot <- NULL
haz.plt <- NULL
lamda.all <- NULL 
time.seq <- seq(0,max(df$time), by = 0.01)
lambda_res_final <-lambda_res_curr<- NULL
changepoint.table <- table(changepoint.num)
indiv.log.lik.final <- indiv.log.lik <-NULL
lambda.list <- list()

plot_on = FALSE 


for(i in seq_along(unique(changepoint.num))){
  #print(i)
  #Get Log Lik
  
  index <- unique(changepoint.num)[order(unique(changepoint.num))][i]
  
  if(changepoint.table[i] <2){
    next
  }
  
  if(index == 0){
    
    lambda_curr <- lambda[which(changepoint.num == index),1:(index+1)]
    lambda_res_final <- matrix(rep(lambda_curr, each = length(time.seq)), 
                               nrow = length(lambda_curr), ncol =  length(time.seq), byrow = T) 
    
    indiv.log.lik.final <- matrix(NA, nrow = length(lambda_curr), ncol = nrow(df))
    
    for(q in 1:nrow(df)){
      indiv.log.lik.final[,q] <- sapply(lambda_curr, FUN = ind.expo, time = df$time[q],status = df$status[q] )
    }
    if(plot_on == TRUE){
      df.changepoint<- data.frame(timepoints = rep(c(0, max(df$time)), times = length(lambda_curr)), 
                                  hazards = rep(lambda_curr, each = 2),
                                  id = rep(1:length(lambda_curr), each = 2))
    }
    
  }else{
    
    changepoint_curr <- changepoint[which(changepoint.num == index),1:index]
    lambda_curr <- lambda[which(changepoint.num == index),1:(index+1)]
    indiv.log.lik <- matrix(NA, nrow = nrow(lambda_curr), ncol = nrow(df))
    
    for(x in 1:nrow(lambda_curr)){
      indiv.log.lik[x,] <- piecewise_loglik.indiv(df, as.numeric(data.frame(changepoint_curr)[x,]), lambda_curr[x,])
      
    }
    
    indiv.log.lik.final <- rbind(indiv.log.lik.final,indiv.log.lik)
    
    
    if(plot_on == TRUE){
      lambda_res_curr <- NULL
      changepoint_curr_samp <- cbind(changepoint_curr, Inf)
      indiv.log.lik <- matrix(NA, nrow = length(lambda_curr), ncol = nrow(df))
      
      for(j in 1:length(time.seq)){
        index.lambda<- apply(changepoint_curr_samp, 1, function(x){which.min( time.seq[j]>x)})
        lambda_res_curr <- cbind(lambda_res_curr,lambda_curr[cbind(1:nrow(changepoint_curr_samp),index.lambda)])
      }
      lambda.list[[i]] <- lambda_res_curr
      print(paste0("NROW lambda final",nrow(lambda_res_final), "NROW lambda current",nrow(lambda_res_curr )))
      lambda_res_final <- rbind(lambda_res_final,lambda_res_curr)
      print(nrow(lambda_res_final))
      change_vec <- t(changepoint_curr)
      change_vec2 <- data.frame(rbind(0,change_vec, max(df$time)))
      lambda_vec <- t(lambda_curr)
      lambda_vec2 <- data.frame(rbind(lambda_vec,lambda_curr[,index+1]))
      
      if(length(changepoint.table[i-1]) != 0){
        index.start <- sum(changepoint.table[1:(i-1)])+1 
      }else{
        index.start <- 1
      }  
      
      df.changepoint<- data.frame(timepoints = unlist(change_vec2), 
                                  hazards = unlist(data.frame(lambda_vec2)),
                                  id = rep(index.start:(index.start+ncol(change_vec)-1), each = (index+2)))
    }
    
  }
  if(plot_on == TRUE){
    
    df.plot <- rbind(df.plot,df.changepoint)
  }
}

if(plot_on == TRUE){
  
  lambda_res_final_df <- data.frame(t(apply(lambda_res_final,2,  FUN = quantile, probs = c(0.05, 0.5, 0.95))), time = time.seq)
  
  
  
  
  max.num.post <- 500
  if(max.num.post < nrow(changepoint)){
    post_id <-  sample(1:nrow(changepoint), size = max.num.post)
    df.plot.final<- dplyr::filter(df.plot,id %in% post_id)
    
  }
  
  interval <- 0.1
  plot_haz <- ggplot(df.plot.final, aes(timepoints, hazards))+ 
    geom_step(aes(group = id), linetype = "dashed", alpha = 0.05, colour = "red")+
    geom_line(data = lambda_res_final_df, aes(time, X50.), size = 1.5)+
    geom_line(data = lambda_res_final_df, aes(time, X5.),linetype = "dotted", size = 1.25)+
    geom_line(data = lambda_res_final_df, aes(time, X95.), linetype = "dotted",size = 1.25)+
    ylim(c(0,1.25))+
    scale_x_continuous(breaks = seq(0, max(df$time), by = interval))+
    theme_classic()
  
  plot_haz
}



indiv.lik.final <-exp(indiv.log.lik.final) 
PML <-1/indiv.lik.final
i <- which.max(apply(PML, 2, var))

y <- seq_along(PML[,i])/cumsum(PML[,i]) 
plot(y = y, x = 1:length(y))
#PML changes with the model chosen the way the datapoints are load
#makes the issue look worse
#... as long as enough sims are run it should be fine
plot(PML[,i])
# If we randomize it this issue isn't problematic.

PML <- PML[sample(nrow(PML),replace = FALSE ),]
i <- which.max(apply(PML, 2, var))

y <- seq_along(PML[,i])/cumsum(PML[,i]) 
plot(y = y, x = 1:length(y))


PML.final <- nrow(PML)/colSums(PML)




log_PML_example <- c(sum(log(PML.final)), 
                 sum(log(PML.expo)),
                 sum(log(PML.weib)),
                 sum(log(PML.llogis)),
                 sum(log(PML.lnorm)),
                 sum(log(PML.gomp)),
                 sum(log(PML.gen.gamma)))

WAIC_example <- c(waic(indiv.log.lik.final)$estimates[3,1],
                  waic(log(Like.sims.expo))$estimates[3,1],
                  waic(log(Like.sims.weib))$estimates[3,1],
                  waic(log(Like.sims.llogis))$estimates[3,1],
                  waic(log(Like.sims.lnorm))$estimates[3,1],
                  waic(log(Like.sims.gomp))$estimates[3,1],
                  waic(log(Like.sims.gen.gamma))$estimates[3,1])

cbind(c(NA,round(-2*log_PML_example, digits = 2)),
      c(NA,round(WAIC_example, digits = 2)),
      round(AUC, digits = 2),
      round(AUC_diff, digits = 2))
cbind(round(-2*log_PML_example, digits = 2),round(WAIC_example, digits = 2))
#plot_Sur_all <- test$plot_Sur_all
# plot_Sur_all[["layers"]][[1]] <- NULL
# plot_Sur_all[["layers"]][[2]] <- NULL
# plot_Sur_all[["layers"]][[4]] <- NULL
# plot_Sur_all[["layers"]][[4]] <- NULL
# 
# plot_Sur_all[["layers"]][[3]][["aes_params"]][["colour"]] <- "black"
# plot_Sur_all

df_surv_expo <- data.frame(Surv.expo,t_pred)
df_surv_weib <- data.frame(Surv.weibull,t_pred)
df_surv_llogis <- data.frame(Surv.llogis,t_pred)
df_surv_lnorm <- data.frame(Surv.lnorm,t_pred)
df_surv_gomp <- data.frame(Surv.gomp,t_pred)
df_surv_gen.gamma <- data.frame(Surv.gen.gamma,t_pred)



colors <- c("KM curve" = "black",
            "Piecewise Expo" = "purple",
            "Exponential" = "red",
            "Weibull"= "cyan",
            "Log-Logistic"= "blue",
            "Log-Normal" = "pink",
            "Gompertz" = "green",
            "Gen. Gamma" = "orange")
#https://stackoverflow.com/questions/13701347/force-the-origin-to-start-at-0

plot_Sur_all <- ggplot(data = mean.Surv_df_all, aes(x = time, y = survival))+
  geom_line(aes(col = "Piecewise Expo"))+
  geom_line(df_surv_expo,mapping =aes(y = Surv.expo, x = t_pred,col = "Exponential"))+
  geom_line(df_surv_weib,mapping =aes(y = Surv.weibull, x = t_pred, col = "Weibull"))+
  geom_line(df_surv_llogis,mapping =aes(y = Surv.llogis, x = t_pred, col = "Log-Logistic"))+
  geom_line(df_surv_lnorm,mapping =aes(y = Surv.lnorm, x = t_pred, col = "Log-Normal"))+
  geom_line(df_surv_gomp,mapping =aes(y = Surv.gomp, x = t_pred, col = "Gompertz"))+
  geom_line(df_surv_gen.gamma,mapping =aes(y = Surv.gen.gamma, x = t_pred, col = "Gen. Gamma"))+
geom_step(data = true_km.df, aes(x = time, y = survival, col = "KM curve"))+
  scale_color_manual(values = colors)+
  labs(x = "Time",
       y = "Survival",
       color = "Survival Functions")+
  scale_y_continuous(expand = c(0, 0),breaks = seq(0, 1, by = 0.1), limits = c(0,1))+
  scale_x_continuous(expand = c(0, 0),breaks = seq(0, 10, by = 1), limits = c(0,NA))+
  theme_classic()+
  theme(legend.position = c(0.9, 0.8))+
  geom_vline(xintercept=c(2), linetype="dotted")


ggsave("Stanford full.png")


partial_km <- survfit(Surv(time,status)~1,data = stanford_test)
partial_km.df <- data.frame(cbind(partial_km[[c("time")]],partial_km[[c("surv")]], partial_km[[c("upper")]],partial_km[[c("lower")]]))
colnames(partial_km.df) <- c("time", "survival", "upper", "lower")


plot_Sur_all+
  geom_step(data = partial_km.df, aes(x = time, y = survival), col = "black")+
  xlim(c(0, 2.0))
ggsave("Stanford partial.png")


plot_Sur_all +theme(legend.position="bottom", legend.box = "horizontal")

