# 1. Glioblastoma Example ----------------------------------------------------
## 1.1 Load data and Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("curatedBreastData")
library("survminer")

#https://bioconnector.github.io/workshops/r-survival.html

# Load the bioconductor installer.
# Try http:// if https:// doesn't work.

#BiocManager::install("RTCGA", lib = "C:/Users/phili/Documents/R/win-library/4.0")
# Install the main RTCGA package
library("RTCGA")
#BiocManager::install("RTCGA.clinical")
# Create the clinical data


#https://bioconnector.github.io/workshops/r-survival.html#tcga
#Survival data.
library("RTCGA.clinical")
library("survival")
library("bshazard")
library("muhaz")
library("survHE")
library("sfsmisc")
devtools::load_all()


clinfinal <- survivalTCGA(GBM.clinical,
                          extract.cols="admin.disease_code")

clinfinal$time <- clinfinal$times#/365
clinfinal[which(clinfinal$time ==0),"time"] <- 0.001
clinfinal$status <- clinfinal$patient.vital_status

surv.fit_glio <- survfit(Surv(time, status)~1, data=clinfinal)

ggsurvplot(surv.fit_glio)


## 1.2 Fit Collapsing Model ----

Collapsing_Model <- collapsing.model(clinfinal,
                                     n.iter = 20750,
                                     burn_in = 750,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1/365)



## 1.3 Get the various Hazard Functions ----
clinfinal <- clinfinal %>% mutate(time_years = time/365,
                                        enter = 0)

bshazard.fit <- bshazard(formula = Surv(enter, time_years,
                                        status) ~ 1,
                         data = clinfinal,
                         nbin = length(unique(clinfinal[which(clinfinal$status == 1),"time"])))


plot(bshazard.fit)

bshazard.df  <- data.frame(time = bshazard.fit[["time"]],
                           hazard = bshazard.fit[["hazard"]],
                           upper.ci = bshazard.fit[["upper.ci"]],
                           lower.ci = bshazard.fit[["lower.ci"]])



harzard.plt.mod <- function(time.vec, cens.vec, Outcome = "Survival",
                            lrg.haz.int = 1,  bw.method = "local"){

  #Plot the hazards
  max.time <- max(time.vec)
  result.hazard.pe.lrg <- pehaz(time.vec, cens.vec, width= lrg.haz.int, max.time=max.time)

  result.smooth <- muhaz(time.vec, cens.vec, b.cor="left",
                         max.time=max.time, bw.method = bw.method)

  #Plot the Survival function

  plot(result.hazard.pe.lrg,  col="black", main = paste0("Hazards for ", Outcome))
  #result.hazard.pe.sml <- pehaz(time.vec, cens.vec, width=sml.haz.int.inv, max.time=max.time)
  #lines(result.hazard.pe.sml,  col="blue")

  lines(result.smooth, col = "red", lty= 2)
  legend("topright", legend = c(paste0(lrg.haz.int," year hazard step function"),
                                #paste0(sml.haz.int," month hazard step function"),
                                "Smoothed hazard function"),
         col = c("black",
                 #"blue",
                 "red"),
         lty = c(1,
                 #2,
                 2), cex = 0.8)

  return(result.hazard.pe.lrg)
}


hazard.PFS.mod <-harzard.plt.mod(time.vec = clinfinal$time_years ,
                                 cens.vec =  clinfinal$status ,
                                 Outcome = "Death", lrg.haz.int = 1)

df.haz <- data.frame(time = hazard.PFS.mod$Cuts[-1]-1,
                     hazard = c(hazard.PFS.mod$Hazard))


View(Collapsing_Model)
haz_plot_collapsing <- plot(Collapsing_Model, type = "hazard")
haz_df <- haz_plot_collapsing$data %>% mutate(timepoints = timepoints/365,
                                              hazards = hazards*365)
haz_df_summary <- haz_df %>% group_by(timepoints) %>% summarize(hazards.mean = mean(hazards),
                                                               hazards.975 = quantile(hazards, 0.975),
                                                               hazards.025 = quantile(hazards, 0.025))

df.summary[nrow(df.summary),1] <-
path_latex <- "C:/Users/phili/OneDrive/PhD/Latex Folder/"
ggplot(bshazard.df[which(bshazard.df$time <= 8.5),], aes(x = time, y = hazard))+
  geom_line(colour = "blue", size = 1.3,linetype = "longdash")+
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.2)+
  geom_step(data = df.haz[which(df.haz$time <= 8.5),],  aes(time, hazard), size = 1.1,linetype = "dotdash")+
  scale_x_continuous(breaks=0:9)+
  xlab("Time")+
  ylab("Hazard")+
  ggtitle("Hazard Functions for Gioblastoma data")+
  geom_step(data = haz_df ,aes(timepoints, hazards,group = id), linetype = "dashed", alpha = 0.03, colour = "red")+
  geom_line( aes(timepoints, hazards.mean ),data = haz_df_summary, colour = "purple", size = 1.3)+
  scale_y_continuous(breaks = seq(from = 0, to = 1.7, by =0.1))+
  coord_cartesian(ylim = c(0,1.7),xlim = c(0,8.5))+
  theme_light()
#Non-parametric Hazards
ggsave(paste0(path_latex,"Hazards Gioblastoma.png"), width = 7, height = 5)




# 2. Stanford Example ----------------------------------------------------

## 2.1 Load data and fit collapsing model and standard Parametric models -----

stanford2 <- survival::stanford2
stanford2$time <- stanford2$time
surv.stanford <- survfit(Surv(time,status)~1,data = stanford2)

dist.res  <- flexsurvreg(Surv(time,status)~1,data = stanford2, dist = "gengamma")

cum_haz_stan <- survminer::ggsurvplot(surv.stanford, risk.table = TRUE, fun = "cumhaz",break.time.by =1)
#ggsave(file = "cum_haz_stan.png", print(cum_haz_stan))
#ggsave("cum_haz_stan.png")

time_censor <- 2#365 #2years

stanford_test <- stanford2 %>% mutate(time = time/365) %>% rename(time_event = time) %>%
                               mutate(time = ifelse(time_event <= time_censor, time_event, time_censor),
                                      status = ifelse(time_event <= time_censor, status, 0))

Collapsing_Model <- collapsing.model(stanford_test,
                                     n.iter = 20750,
                                     burn_in = 750,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1)



mod_comp <-compare.surv.mods(Collapsing_Model,
                             max_predict = 10,
                             plot.best =2,
                             chng.num = 2,
                             n.iter.jags = 10000,
                             n.thin.jags = 1,
                             n.burnin.jags = 1000)


Piecewise_LL<- get.loglik(Collapsing_Model)
index_2 <- apply(Collapsing_Model$k.stacked, 1, function(x){length(na.omit(x))}) ==2

#2 changepoint + 3 hazards
num.param_piece <- mean(2*2 +1)
num.param <- c(1,2,2,2,2,2,3)

AIC_vec <- BIC_vec <- rep(NA, 9)

#Brooks Approach
# total_piecewiseLL <- mean(rowSums(Piecewise_LL[index_2,]))
# AIC_vec[8] <- -2*total_piecewiseLL + 2*num.param[1]
# BIC_vec[8] <- -2*total_piecewiseLL + num.param[1]*log(sum(stanford_test$status))

#Raftery Approach
LL_max <- mean(rowSums(Piecewise_LL[index_2,])) +var(rowSums(Piecewise_LL[index_2,]))
AIC_vec[8] <- -2*LL_max + 2*num.param_piece
BIC_vec[8] <- -2*LL_max + num.param_piece*log(sum(stanford_test$status))


mod.flexsurv <- c("exp", "weibullPH","gamma",  "lnorm", "llogis","gompertz", "gengamma.orig")

for(i in 1:length(mod_comp[["jag.models"]])){
  index <- grep("total_LLik",rownames(mod_comp[["jag.models"]][[i]][["BUGSoutput"]][["summary"]]))
#Brooks approach
# minus_2LL<- -2*mod_comp[["jag.models"]][[i]][["BUGSoutput"]][["summary"]][index,1]
#  AIC_vec[i] <- minus_2LL + 2*num.param[i]
#  BIC_vec[i] <- minus_2LL + num.param[i]*log(sum(stanford_test$status))

#Raftery Approach
  LL_max <- mod_comp[["jag.models"]][[i]][["BUGSoutput"]][["summary"]][index,1] + (mod_comp[["jag.models"]][[i]][["BUGSoutput"]][["summary"]][index,2])^2

  # MLE_true <- flexsurv::flexsurvreg(formula=Surv(time,status)~1,
  #                                   data=stanford_test,
  #                                   dist = mod.flexsurv[i])
  #print(c(LL_max, MLE_true$loglik))

  AIC_vec[i] <- -2*LL_max + 2*num.param[i]
  BIC_vec[i] <- -2*LL_max + num.param[i]*log(sum(stanford_test$status))
}


## 2.2 Fit Royston Parmer models -----

param_expert_vague <- list()
param_expert_vague[[1]] <- data.frame(dist = "beta", wi = 1, param1 = 1, param2 = 1, param2 = NA)

rps.1  <- expertsurv::fit.models.expert(formula=Surv(time,status)~1,data=stanford_test,
                                        distr=c("rps"),
                                        k = 1,
                                        method="hmc",
                                        iter = 5000,
                                        warmup = 200,
                                        opinion_type = "survival",
                                        times_expert = 1,
                                        param_expert = param_expert_vague)


LL_max <- mean(rstan::extract(rps.1$models$`Royston-Parmar`)[["lp__"]]) + var(rstan::extract(rps.1$models$`Royston-Parmar`)[["lp__"]])

parm_rps.1 <-3
AIC_vec[9] <- -2*LL_max + 2*parm_rps.1
BIC_vec[9] <- -2*LL_max + parm_rps.1*log(sum(stanford_test$status))

#flexsurv::flexsurvspline(formula=Surv(time,status)~1, data=stanford_test, k = 1)
#Brooks approach

#AIC_vec[9] <- rps.1$model.fitting$aic
#BIC_vec[9] <- mean(rstan::extract(rps.1$models$`Royston-Parmar`)[["lp__"]]) +parm_rps.1*log(sum(stanford_test$status))

#rps.1$model.fitting$aic - 2*3 + 3*log(length(stanford_test$status))


# 1 Knot fits best

# rps.2  <- expertsurv::fit.models.expert(formula=Surv(time,status)~1,data=stanford_test,
#                                         distr=c("rps"),
#                                         k = 2,
#                                         method="hmc",
#                                         iter = 5000,
#                                         warmup = 200,
#                                         opinion_type = "survival",
#                                         times_expert = 1,
#                                         param_expert = param_expert_vague)


mod_comp2 <- rbind(mod_comp$mod.comp,c(NA,NA,NA))
mod_comp2[nrow(mod_comp2),1] <- "Royston-Parmer"
mod_comp2[nrow(mod_comp2),2] <- rps.1[["model.fitting"]]$pml
mod_comp2[nrow(mod_comp2),3] <- rps.1[["model.fitting"]]$waic

## 2.3 Add RP to Survival Curve -----


max_predict <- 10
psa_rps.1 <- make.surv(fit = rps.1, nsim = 1000, t = seq(0, 10, by = max_predict/100))


rps_St_df <- data.frame(time = psa_rps.1$mat[[1]][,1] %>% pull(),
                          St_rps.1 = rowMeans(psa_rps.1$mat[[1]][,-1])
                        #,St_rps.2 = rowMeans(psa_rps.2$mat[[1]][,-1])
                        )

col_vec <-    c(
  "black",
  "purple",
  "red",
  "cyan",
  "brown",
  "blue",
  "pink",
  "green",
  "orange")

col_name <- c(
  "KM curve",
  "Piecewise Expo",
  "Exponential" ,
  "Weibull" ,
  "Gamma",
  "Log-Logistic",
  "Log-Normal" ,
  "Gompertz" ,
  "Gen. Gamma"
)

colors <- col_vec
names(colors) <- col_name

col_vec_incl <- colors[names(colors) %in% c("KM curve","Piecewise Expo","Log-Logistic","Log-Normal")]
col_vec_incl <- c(col_vec_incl, "Royston-Parmer" = "red")

plot_all <- mod_comp$plot_Surv_all+
  geom_line(data = rps_St_df, aes(x = time, y = St_rps.1, colour = "Royston-Parmer"), inherit.aes = F)+
  scale_color_manual(values = col_vec_incl)


add_km(plot_all, df = stanford2 %>% mutate(time = time/365))

ggsave(filename = "Stanford Survival.png",width = 8, height = 4.24)

#Frequentist (In case Bayesian appproach doesn't fit)
# sp1 <- flexsurvspline(Surv(time, status) ~ 1, data = stanford_test, k = 1,
#                       scale = "hazard")
# sp1_flexsurv <- data.frame(summary(sp1, t = seq(0, 10, by = 0.01)))


St_piecewise <- get_Surv(Collapsing_Model,chng.num = 2, max_predict = 10)
St_piecewise_mean <- rowMeans(St_piecewise)
t <- seq(0, 10, 10/100)

true_km <- survfit(Surv(time,status)~1,data = stanford2 %>% mutate(time = time/365))
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

#survival:::survmean(true_km, rmean=10)[[1]]["*rmean"]

mod_names <- c("expo", "weibull", "gamma", "llogis", "lnorm", "gomp", "gen.gamma")

for(i in 1:length(mod_names)){

  temp_mod <- mod_comp$jag.models[[i]]$BUGSoutput[["summary"]]
  assign(paste0("Surv.",mod_names[i]),temp_mod[grep("St_pred",rownames(temp_mod)),1])


}


AUC = c(
        integrate.xy(t, Surv.expo),
        integrate.xy(t, Surv.weibull),
        integrate.xy(t, Surv.gamma),
        integrate.xy(t, Surv.lnorm),
        integrate.xy(t, Surv.llogis),
        integrate.xy(t, Surv.gomp),
        integrate.xy(t, Surv.gen.gamma),
        integrate.xy(t, St_piecewise_mean),
        integrate.xy(t, rps_St_df[,"St_rps.1"]),AUC_true)

AUC_diff = c(
  integrate.xy(t, abs(Surv.expo - surv_km_time[1:length(t)])),
  integrate.xy(t, abs(Surv.weibull - surv_km_time[1:length(t)])),
  integrate.xy(t, abs(Surv.gamma - surv_km_time[1:length(t)])),
  integrate.xy(t, abs(Surv.lnorm - surv_km_time[1:length(t)])),
  integrate.xy(t, abs(Surv.llogis - surv_km_time[1:length(t)])),
  integrate.xy(t, abs(Surv.gomp - surv_km_time[1:length(t)])),
  integrate.xy(t, abs(Surv.gen.gamma - surv_km_time[1:length(t)])),
  integrate.xy(t, abs(St_piecewise_mean-surv_km_time[1:length(t)])),
  integrate.xy(t, abs(rps_St_df[,"St_rps.1"] - surv_km_time[1:length(t)])),
  NA)

mod_compfinal <- rbind(cbind(mod_comp2,AIC_vec,BIC_vec),c(NA,NA,NA))
mod_compfinal$Model[nrow(mod_compfinal)] <- "True Observations"
mod_compfinal <- cbind(mod_compfinal, AUC, AUC_diff)

mod_compfinal<- mod_compfinal[order(mod_compfinal$AUC_diff),] %>% mutate_if(is.numeric, round, digits = 2)

print(xtable::xtable(mod_compfinal, type = "latex"), include.rownames=FALSE)

df_IC <- round(data.frame(AIC = mod_compfinal$AIC_vec, BIC =  mod_compfinal$BIC_vec), 2)
df_IC <- cbind("&",df_IC$AIC, "&", df_IC$BIC )
paste(df_IC[,1],df_IC[,2],df_IC[,3],df_IC[,4] )
# 3. Leukemia Example ----------------------------------------------------

library("openxlsx")
leukemia.df <- read.xlsx(paste0("C:\\Users\\phili\\OneDrive\\PhD\\R Codes\\Changepoint Analysis\\","Achar Leukemia Data.xlsx"),1)
leukemia.df <- leukemia.df[which(leukemia.df[,1] !=182),]


## 3.1 Fit Collapsing Model ------------



Collapsing_Model_leukemia <- collapsing.model(leukemia.df,
                                     n.iter = 20750,
                                     burn_in = 750,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1/365)

## 3.2 Fit Gibbs Samplers ------------

gibbs_1chng <- PiecewiseChangepoint:::gibbs.sampler(df = leukemia.df,
                                                    num.breaks = 1, # Number of Changepoints
                                                    n.iter = 5000,
                                                    burn_in = 100,
                                                    num.chains = 2,
                                                    n.thin = 1,
                                                    alpha.hyper = 1, #Hyperparameter for Marginal Likelihood
                                                    beta.hyper = 365, #Hyperparameter for Marginal Likelihood
                                                    MLE = FALSE)
plt_surv1 <- PiecewiseChangepoint:::plot.changepoint(gibbs_1chng, max_predict = 2000)


gibbs_2chng <- PiecewiseChangepoint:::gibbs.sampler(df = leukemia.df,
                          num.breaks = 2, # Number of Changepoints
                          n.iter = 5000,
                          burn_in = 100,
                          num.chains = 2,
                          n.thin = 1,
                          alpha.hyper = 1, #Hyperparameter for Marginal Likelihood
                          beta.hyper = 365, #Hyperparameter for Marginal Likelihood
                          MLE = FALSE)
  #If TRUE, considers on Uniform Prior
plt_surv2 <- PiecewiseChangepoint:::plot.changepoint(gibbs_2chng, max_predict = 2000)

ggarrange(plt_surv1+xlab("Time"),
          plt_surv2+xlab("Time")+ylab(""))
ggsave(paste0(path_latex,"Leukemia_both.png"),width=12,height =4)





# 4 Comparison of the Collapsing approach with Marginal Likelihood evaluation ------------

set.seed(123)
n_obs =100
n_events_req=100
max_time =  2

rate = c(0.75,0.25)
t_change =1

df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                   num.breaks = length(t_change),rate = rate ,
                   t_change = t_change, max_time = max_time)



Collapsing_Model <- collapsing.model(df,
                                     n.iter = 50000,
                                     burn_in = 750,
                                     n.chains = 3,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1)
#remove.packages("PiecewiseChangepoint")

library("Brobdingnag")
library("combinat")

# generating combinations of the
# alphabets taking 2 at a time

beta.hyper <- 1
time_diffs<- PiecewiseChangepoint:::df_recast(df)
#n.sims <- 100000
n.sims <- 10000
beta_vals <- 400
marg_vec2 <- marg_vec <- c()
for(k in 1:3){

  k_change <- k
  n_events <- nrow(time_diffs)
  combs <- combn(1:n_events,k_change)
  probs <-apply(combs,2,function(x){PiecewiseChangepoint:::pos_prob(x, n_events, FALSE)})


  #PiecewiseChangepoint:::pos_prob(c(2,4,8), 10, FALSE)

  index.drp <- which(probs == 0)
  probs <- probs[-index.drp]
  combs <- combs[,-index.drp, drop = F]
  cum_prob <- cumsum(probs)
  max(cum_prob)
  unif_sims <- runif(n.sims)
  k_prior_loc <- sapply(unif_sims, function(x){min(which(cum_prob >= x))})


  # marg_vec <- c(marg_vec,log(mean(exp(apply(combs[,k_prior_loc,drop = F], 2,
  #                    function(x){PiecewiseChangepoint::margin_lik_calc_log(time_diffs,x,1,365)+
  #                        log(PiecewiseChangepoint:::pos_prob(x, n_events, FALSE))})))))

  # marg_vec2 <- c(marg_vec2,log(mean(exp(apply(combs[,k_prior_loc,drop = F], 2,
  #                                           function(x){PiecewiseChangepoint::margin_lik_calc_log(time_diffs,x,1,500)})))))



  beta_array <- array(dim = c(dim(combs)[1]+1, dim(combs)[2],beta_vals))

  for(i in 1:beta_vals){
    beta_array[,,i] <- rbind(combs, rgamma(dim(combs)[2],1, 1/beta.hyper))

  }


  res <- matrix(nrow = n.sims, ncol = beta_vals)

  for(j in 1:beta_vals){

    res[,j] <-  apply(beta_array[,k_prior_loc,j,drop = F],2,
                      function(x){PiecewiseChangepoint:::margin_lik_calc_log(time_diffs,x[1:k],1,x[k+1])})

  }

  brob_res <- exp(as.brob(as.vector(res)))
  marg_vec2 <- c(marg_vec2, sum(brob_res)/length(brob_res))

}



#NAIVE ESTIMATION
# lambda_naive  <- rgamma(100000, 1, 500)
# log(mean((lambda_naive^sum(df$status))*exp(-lambda_naive*sum(df$time))))

res_no_chng <- c()


for(i in 1:beta_vals){
  res_no_chng <- c(res_no_chng, PiecewiseChangepoint:::margin_lik_calc_log(time_diffs,NA,1,rgamma(1, 1, 1/beta.hyper)))
}

brob_res_no_chng <- exp(as.brob(as.vector(res_no_chng)))
marg_no_chng <- sum(brob_res_no_chng)/length(brob_res_no_chng)
marg.like <- c(marg_no_chng,marg_vec2)

#exp(marg_vec2[1]-margin_lik_calc_log(time_diffs,NA,1,500))

log_marg <- unlist(lapply(marg.like,function(x){log(x)}))

options(scipen = 100)
log_marg2 <- log_marg+ dpois(0:(length(log_marg)-1), 1, log = T)
round(exp(diff(log_marg2)), digits = 4)

Collapsing_Model$prob.changepoint[2]/Collapsing_Model$prob.changepoint[1]
Collapsing_Model$prob.changepoint[3]/Collapsing_Model$prob.changepoint[2]
Collapsing_Model$prob.changepoint[4]/Collapsing_Model$prob.changepoint[3]


# We show that we can evaluate the marginal likelihood by brute force and obtain similar Bayes Factors although much quicker


# 5 Convergence assessment of Collapsing Models ------------




RJMCMC_conver_plt <- function(k,lambda, batch.size = 50, test.index){

  n.chains <- dim(k)[3]
  complete.batch <- floor(nrow(k[,,1])/batch.size)
  nrow.array <- complete.batch*batch.size
  lambda <- lambda[1:nrow.array,,]
  k <- k[1:nrow.array,,]
  chain.id <- rep(1:n.chains, each = nrow.array)
  batch.id <- rep(1:complete.batch, each = batch.size )

  k.stacked.converge <- k[,,1]
  lambda.stacked.converge <- lambda[,,1]
  #changepoint.stacked.converge <- changepoint[,,1]

  for(i in 2:n.chains){
    k.stacked.converge  <-rbind(k.stacked.converge,k[,,i])
    lambda.stacked.converge  <-rbind(lambda.stacked.converge,lambda[,,i])
    #changepoint.stacked.converge <-rbind(changepoint.stacked.converge,changepoint[,,i])
  }

  converge.array <-   array(NA,dim = c(nrow.array*n.chains,5))
  converge.array[,1] <- apply(k.stacked.converge,1,
                              FUN= index.loc,index = test.index)

  for(i in 1:nrow(converge.array)){
    converge.array[i,2] <- lambda.stacked.converge[i,converge.array[i,1]]
  }
  num.changepoints<- apply(k.stacked.converge,1,function(x){length(na.omit(x))})
  converge.array[,3] <- num.changepoints
  converge.array[,4] <- batch.id
  converge.array[,5] <- chain.id
  converge.array <- data.frame(converge.array)
  colnames(converge.array) <- c("Location","Hazard", "NumChangepoint", "Batch","Chain")


  PSRF1_1 <- rep(NA,complete.batch)
  PSRF1_2 <- rep(NA,complete.batch)

  PSRF2_1 <- rep(NA,complete.batch)
  PSRF2_2 <- rep(NA,complete.batch)

  MPSRF1 <- rep(NA,complete.batch)
  MPSRF2  <- rep(NA,complete.batch)
  eigen1 <- rep(NA,complete.batch)
  eigen2 <- rep(NA,complete.batch)
  eigen3 <- rep(NA,complete.batch)
  eigen4 <- rep(NA,complete.batch)


  for(i in 1:complete.batch){

    temp.converge.array <- converge.array %>% dplyr::filter(Batch %in% c(1:i))
    T.val <- batch.size*i
    M.val <- length(unique(temp.converge.array$NumChangepoint))

    eq.11 <- temp.converge.array %>% group_by(Chain,NumChangepoint) %>%
      dplyr::summarise_all(mean,na.rm =T) %>% dplyr::select(-Batch)
    names(eq.11) <- c("Chain","NumChangepoint","Location_avg", "Hazard_avg")

    eq.12 <- temp.converge.array %>% group_by(Chain) %>% dplyr::summarise_all(mean,na.rm =T)
    eq.12 <- eq.12[,c(1:3)]
    names(eq.12) <- c("Chain","Location_avg", "Hazard_avg")

    eq.13 <- temp.converge.array %>% group_by(NumChangepoint) %>%
      dplyr::summarise_all(mean,na.rm =T)
    eq.13 <- eq.13[,c(1:3)]
    names(eq.13) <- c("NumChangepoint","Location_avg", "Hazard_avg")

    eq.14  <- colMeans(temp.converge.array) #Equation 14 in Castelloe et al
    eq.14 <- eq.14[c(1,2)]

    #Multivariate & univariate
    eq.37 <- var(temp.converge.array[,c(1,2)])
    eq.15 <- diag(var(temp.converge.array[,c(1,2)]))

    eq.38 <- temp.converge.array[,-c(3,4)] %>% group_by(Chain) %>% dplyr::left_join(eq.12)%>%
      mutate(diff_haz = Hazard-Hazard_avg,
             diff_loc = Location-Location_avg)

    eq.16 <- colSums((eq.38[,c("diff_loc","diff_haz")])^2)/(n.chains*(T.val -1))

    eq.38 <- matrix(rowSums(apply(eq.38[,c("diff_loc","diff_haz")],1,FUN = function(x){x%*%t(x)}))/
                      (n.chains*(T.val -1)),2,2)

    eq.39 <- temp.converge.array[,-c(4,5)]  %>% dplyr::left_join(eq.13, by = "NumChangepoint") %>%
      mutate(diff_haz = Hazard-Hazard_avg,
             diff_loc = Location-Location_avg)

    eq.17 <- colSums((eq.39[,c("diff_loc","diff_haz")]^2))/(n.chains*T.val -M.val)

    eq.39 <- matrix(rowSums(apply(eq.39[,c("diff_loc","diff_haz")],1,
                                  FUN = function(x){x%*%t(x)}))/(n.chains*T.val -M.val),2,2)

    eq.40 <- temp.converge.array[,-c(4)]  %>% dplyr::left_join(eq.11, by = c("Chain","NumChangepoint")) %>%
      mutate(diff_haz = Hazard-Hazard_avg,
             diff_loc = Location-Location_avg)

    eq.18 <- colSums((eq.40[,c("diff_loc","diff_haz")]^2))/(n.chains*(T.val -M.val))

    eq.40 <- matrix(rowSums(apply(eq.40[,c("diff_loc","diff_haz")],1,FUN = function(x){x%*%t(x)}))/(n.chains*(T.val -M.val)),2,2)

    PSRF1_1[i] <- eq.15[1]/eq.16[1]
    PSRF1_2[i] <- eq.15[2]/eq.16[2]

    PSRF2_1[i] <- eq.17[1]/eq.18[1]
    PSRF2_2[i] <- eq.17[2]/eq.18[2]

    MPSRF1[i] <- max(eigen(solve(eq.38)%*%eq.37)$values)
    MPSRF2[i]  <- max(eigen(solve(eq.40)%*%eq.39)$values)
    eigen1[i] <- max(eigen(eq.37)$values)
    eigen2[i] <- max(eigen(eq.38)$values)
    eigen3[i] <- max(eigen(eq.39)$values)
    eigen4[i] <- max(eigen(eq.40)$values)

  }

  plot(MPSRF1, typ = "l")
  lines(PSRF1_1, col = 2)
  lines(PSRF1_2, col = 3)

  plot(MPSRF2, typ = "l")
  lines(PSRF2_1, col = 2)
  lines(PSRF2_2, col = 3)

  plot(eigen1, typ = "l")
  lines(eigen2, col = 2)

  plot(eigen3, typ = "l")
  lines(eigen4, col = 2)

}

RJMCMC_conver_plt2 <- function(k, batch.size = 50){

  n.chains <- dim(k)[3]
  complete.batch <- floor(nrow(k[,,1])/batch.size)

  nrow.array <- complete.batch*batch.size
  k <- k[1:nrow.array,,]
  chain.id <- rep(1:n.chains, each = nrow.array)
  batch.id <- rep(1:complete.batch, each = batch.size )

  k.stacked.converge <- k[,,1]

  for(i in 2:n.chains){
    k.stacked.converge  <-rbind(k.stacked.converge,k[,,i])
    #changepoint.stacked.converge <-rbind(changepoint.stacked.converge,changepoint[,,i])
  }

  num.changepoints <-unlist(apply(k.stacked.converge,1, function(x){length(na.omit(x))}))

  converge.array <- array(NA,dim = c(nrow.array*n.chains,3))
  converge.array[,1] <- num.changepoints
  converge.array[,2] <- batch.id
  converge.array[,3] <- chain.id
  converge.array <- data.frame(converge.array)
  colnames(converge.array) <- c( "NumChangepoint", "Batch","Chain")
  PSRF1_1 <- rep(NA,complete.batch)


  for(i in 1:complete.batch){

    temp.converge.array <- converge.array %>% dplyr::filter(Batch %in% c(1:i))
    T.val <- batch.size*i
    M.val <- unique(temp.converge.array$NumChangepoint)

    eq.11 <- temp.converge.array %>% group_by(Chain,NumChangepoint) %>%
      dplyr::summarise_all(mean,na.rm =T) %>% dplyr::select(-Batch)
    names(eq.11) <- c("Chain","NumChangepoint")

    eq.15 <- var(temp.converge.array[,1])

    eq.16 <- temp.converge.array%>% group_by(Chain) %>% dplyr::summarize(NumChangepoint_avg = mean(NumChangepoint))
    eq.16 <- temp.converge.array %>% dplyr::left_join(eq.16, by = "Chain") %>% mutate(diff_NumChange = (NumChangepoint-NumChangepoint_avg)^2)
    eq.16 <- sum(eq.16[,"diff_NumChange"])/(n.chains*(T.val -1))

    PSRF1_1[i] <- eq.15/eq.16


  }

  tiff("PSRF1_1.png",  units="in", width=6, height=4,res=80)

  plot(PSRF1_1, type = "l", ylab = "PSRF", xlab = "Batch Number")

  dev.off()

}



k<- Collapsing_Model[["k"]]
lambda <- Collapsing_Model[["lambda"]]

#have to unstack the lambda

seq(from = 1, to = nrow(lambda), by = dim(k)[1] )

upper <- dim(k)[1]*(1:dim(k)[3])
lower <- head(c(1, 1+ upper), n= -1)
dim_lamb <- dim(k)
dim_lamb[2] <- dim_lamb[2]+1
lambda_array <- array(NA, dim = dim_lamb)

for(i in 1:dim(k)[3]){
  lambda_array[,,i] <- lambda[lower[i]:upper[i],]

}

batch.size <- 5000
test.index <- Mode(k[,1,1])

undebug(RJMCMC_conver_plt)
#Have to comment back in the code to define a lambda_array, for some reason this only works with a weibull model.
RJMCMC_conver_plt(k,lambda = lambda_array, batch.size = batch.size, test.index = test.index)
RJMCMC_conver_plt2(k, batch.size = batch.size)



# 6 NICE plot of Changepoint Log-Likelihood  ------------


library("purrr")
library("pch")
library("plotly")


url.path <- "http://merlot.stat.uconn.edu/~mhchen/survbook/dataset/e1690.missing.dat"
E1690.dat <- read.delim(url(url.path), header = T, sep="", skip=12, as.is=TRUE)
#Drop PFS events with time equal zero
E1690.dat <- E1690.dat[-which(E1690.dat$FAILTIME ==0),]

#Convert to the correct notation for survival objects
E1690.dat[which(E1690.dat$FAILCENS == 1),"FAILCENS"] <-0
E1690.dat[which(E1690.dat$FAILCENS == 2),"FAILCENS"] <-1
E1690.dat[which(E1690.dat$SURVCENS == 1),"SURVCENS"] <-0
E1690.dat[which(E1690.dat$SURVCENS == 2),"SURVCENS"] <-1

source("C:/Users/phili/OneDrive/PhD/R packages/PiecewiseChangepoint/PiecewiseChangepoint/Publication Examples/Backup PhD/Analysis functions.R")
trt.df <- E1690.dat[which( E1690.dat$TRT == 2 & E1690.dat$STUDY == "1684"),]


result <- grid.search.piecewise.pchreg(min.break = 0.2,
                                       max.break = 7,
                                       grid.width = 0.25,
                                       num.breaks = 2,
                                       min.break.width = 1,
                                       time = trt.df$FAILTIME,
                                       status = trt.df$FAILCENS)
View(result)



## Add in the simulation studies


library(ggplot2)
library(ggpubr)
library(stringr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
sim.path <- "C:\\Users\\phili\\OneDrive\\PhD\\R codes\\Changepoint Analysis\\Simulation Study\\Simulation Collapsing\\One Changepoint\\"
model.names <- c("decreasing_haz_lrg_","decreasing_haz_sml_","increasing_haz_lrg_","increasing_haz_sml_")
model.size <- c("200","500","1000")

#sim.path <- "C:/Users/phili/OneDrive/PhD/Simulation Study/Simulation Collapsing/"
#model.names <- c("bathtub_haz","decreasing_haz","increasing_haz","invertbath_haz")
#model.size <- c("300","500","1000")

model.com <- rep(NA, length(model.names)*length(model.size))
cnt <- 0
for(i in 1:length(model.names)){
  for(j in 1:length(model.size)){
    cnt <- cnt +1
    model.com[cnt] <- paste0(model.names[i],model.size[j])

  }
}
library(stringr)

n <- 1000

df.prob <- data.frame(probability= rep(NA,12),
                      model= rep(NA,12),
                      n =  rep(NA,12))

dict.names <- cbind(c("Increasing Small", "Increasing Large", "Decreasing Small", "Decreasing Large"),
                    c("increasing_haz_sml_","increasing_haz_lrg_", "decreasing_haz_sml_","decreasing_haz_lrg_"))

for(i in 1:length(model.com)){
  model.name <- model.com[i]
  sim.size <- as.numeric(gsub("[^0-9.]", "",  model.name))
  assign(model.name,loadRData(paste0(sim.path,model.com[i],".RData")))
  assign("current.model",loadRData(paste0(sim.path,model.com[i],".RData")))
  txt.string <- str_remove(model.name,paste0(sim.size))#str_remove(model.name,paste0(sim.size))#

  df.prob[i,1] <- current.model$prob.correct*100 #
  df.prob[i,2] <- dict.names[which(dict.names[,2]==txt.string),1]
  df.prob[i,3] <- sim.size

  for(j in 1:1){
    shape <- ((current.model$haz[j,1])^2)/(current.model$haz[j,2]^2)
    rate <- (current.model$haz[j,1])/(current.model$haz[j,2]^2)
    assign(paste0("df.res",j),data.frame(hazard = rgamma(n, shape = shape, rate = rate),
                                         changepoint = rnorm(n, current.model$chng[j,1], current.model$chng[j,2]),
                                         n = as.factor(sim.size), interval= as.factor(j), model = dict.names[which(dict.names[,2]==txt.string),1]))



  }
  shape <- ((current.model$haz[2,1])^2)/(current.model$haz[2,2]^2)
  rate <- (current.model$haz[2,1])/(current.model$haz[2,2]^2)

  df.res2 <-data.frame(hazard = rgamma(n, shape = shape, rate = rate),
                       changepoint = NA,
                       n = as.factor(sim.size), interval = as.factor(3),
                       model =  dict.names[which(dict.names[,2]==txt.string),1])


  assign(paste0("df.res_all",txt.string,sim.size),rbind(df.res1,df.res2))


  rm(df.res1,df.res2,txt.string)

}

for(i in 1:nrow(dict.names)){

  assign(paste0("df.res_final",dict.names[i,2]),rbind(get(paste0("df.res_all",dict.names[i,2],200)),
                                                      get(paste0("df.res_all",dict.names[i,2],500)),
                                                      get(paste0("df.res_all",dict.names[i,2],1000))))

}

df.prob$n <- as.factor(df.prob$n)
#df.prob$model <- as.factor(df.prob$model)

df.prob$model <- factor(df.prob$model, levels = c("Increasing Small", "Increasing Large", "Decreasing Small", "Decreasing Large"))
prob.plot <- ggplot(df.prob, aes(fill=n, y=probability, x=model)) +
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(breaks = seq(0,100,10))+
  xlab("")+
  ylab("Probability")+
  theme_bw()+
  theme(legend.position = "none", plot.margin = unit(c(1.5,0.5,0.5,0.5), "lines"))
true.haz <- data.frame(model = dict.names[,2],
                       interval1  = c(0.25,0.75,0.75,0.2),
                       interval2  = c(0.5,0.5,0.2,0.75),
                       interval3  = c(0.75,0.25,0.75,0.2))


for(i in 1:nrow(dict.names)){

  df.plt.temp <- get(paste0("df.res_final",dict.names[i,2]))
  df.plt.temp$interval <- as.factor(df.plt.temp$interval)
  df.plt.temp$model <- as.factor(df.plt.temp$model)
  df.plt.temp$n <- as.factor(df.plt.temp$n)

  #index <- which(dict.names[i,2]==true.haz[,1])
  #print(index)
  gg.plt.temp <-  ggplot(df.plt.temp, aes(x=interval, y=hazard, fill=n)) +
    geom_boxplot(position=position_dodge(1), outlier.shape = NA)+
    #geom_segment(aes(x=2/3,xend=4/3,y=true.haz[index,1+1],yend=true.haz[index,1+1]), colour = "black", linetype = 1, size = 1, show.legend   = FALSE)+
    #geom_segment(aes(x=5/3,xend=7/3,y=true.haz[index,2+1],yend=true.haz[index,2+1]), colour = "orange", linetype = 1, size = 1, show.legend   = FALSE) +
    #geom_segment(aes(x=8/3,xend=10/3,y=true.haz[index,3+1],yend=true.haz[index,3+1]), colour = "purple", linetype = 1, size = 1, show.legend   = FALSE)+
    xlab("")+
    scale_y_continuous(breaks = seq(0,1.2,.2))+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5, linetype = "solid"))+
    coord_cartesian(ylim = c(0,1.2),expand = FALSE)
  #print(c(true.haz[index,1+1],true.haz[index,2+1],true.haz[index,3+1]))

  gg.plt.temp
  assign(paste0("p_hazard_",dict.names[i,2]),gg.plt.temp)
  rm(gg.plt.temp)

  gg.plt.temp2  <- ggplot(df.plt.temp[!is.na(df.plt.temp$changepoint),], aes(x=interval, y=changepoint, fill=n)) +
    geom_boxplot(position=position_dodge(1), outlier.shape = NA)+
    #geom_segment(aes(x=2/3,xend=4/3,y=0.5,yend=0.5), colour = "black", linetype = 1, size = 1, show.legend  = FALSE)+
    #geom_segment(aes(x=5/3,xend=7/3,y=1,yend=1), colour = "orange", linetype = 1, size = 1, show.legend   = FALSE)+
    xlab("")+
    ylab("")+
    scale_y_continuous(breaks = seq(.2,1.4,.2))+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    coord_cartesian(ylim = c(.2,1.4),expand = FALSE)

  assign(paste0("p_change_",dict.names[i,2]),gg.plt.temp2)
  rm(gg.plt.temp2)

}

plt.output <- ggarrange(
  ggarrange(p_hazard_increasing_haz_sml_ +ylab("Increasing Small"),
            p_change_increasing_haz_sml_,
            p_hazard_increasing_haz_lrg_ +ylab("Increasing Large"),
            p_change_increasing_haz_lrg_,
            p_hazard_decreasing_haz_sml_ +ylab("Decreasing Small"),
            p_change_decreasing_haz_sml_,
            p_hazard_decreasing_haz_lrg_ +ylab("Decreasing Large"),
            p_change_decreasing_haz_lrg_,
            nrow = 4, ncol = 2, common.legend = TRUE, legend = "bottom"), prob.plot+theme(axis.text=element_text(size=8)))
#Could annotate figure at bottom

ggsave(annotate_figure(plt.output, fig.lab.pos = c("top.left"),fig.lab.face= "bold",fig.lab =  c("Hazard across intervals         Changepoints                           Probability of selecting correct model")) ,
       filename = "Collapsing Model One Changepoint.png", width = 8.5, height =6)



# 7 Plot of Model results from simulation study ------------
#C:\Users\phili\OneDrive\PhD\R codes\Changepoint Analysis\Simulation Study


library(ggplot2)
library(ggpubr)
library(stringr)

## Need to update for  Chapell analysis
## Include mean difference from true survival function as per Chapell

sim.study <- function(n_obs,n_events_req,rate,t_change, max_time =2,
                      n.sims, sims = 10750,burn_in = 750){

  if(is.na(t_change)){
    num.breaks = 0
  }else{
    num.breaks <- length(t_change)

  }


  n.events <- rep(NA, n.sims)
  avg.Haz <- array(NA, dim = c(n.sims, 7))
  avg.change <- array(NA, dim = c(n.sims, 6))

  for(i in 1:n.sims){
    print(paste0("Sim ",i," of ",n.sims))
    Sys.sleep(1 / n.sims)

    n_events <- 0

    while(n_events  != req_obs){
      df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                         num.breaks = num.breaks,rate = rate ,
                         t_change = t_change, max_time = max_time)
      n_events <-  df$status
    }

    n.events[i] <- sum(df$status == 1)
    #tryCatch({
    mod <- RJMC.piece.compact_rcpp(sims, df,n.chains = 1, burn_in = burn_in )
    avg.Haz[i,1:length(mod[[1]])] <- mod[[1]]
    avg.change[i,1:length(mod[[2]])] <- mod[[2]]
    #}, error = function(e){print("Nuts")})
  }

  chang.vals <-apply(avg.change,1, function(x){length(na.omit(x))})
  uniq.chang <- unique(chang.vals)[order(unique(chang.vals))]
  final.chng <- data.frame(mean = rep(NA,length(uniq.chang)),sd = rep(NA,length(uniq.chang)), nsims = rep(NA,length(uniq.chang)))
  for(i in 1:length(uniq.chang)){
    final.chng[i,1]  <- mean(avg.change[which(chang.vals == uniq.chang[i]),uniq.chang[i]], na.rm =T)
    final.chng[i,2] <- sd(avg.change[which(chang.vals == uniq.chang[i]),uniq.chang[i]], na.rm =T)
    final.chng[i,3] <- table(chang.vals)[i]
  }

  row.names(final.chng) <- uniq.chang
  prob.correct <- mean(chang.vals == num.breaks)

  if(length(which(chang.vals == num.breaks)) >1){

    if(num.breaks == 0 ){
      se.chng <- mean.chng <- NA
      mean.haz <- mean(avg.Haz[which(chang.vals == num.breaks),
                               1:(num.breaks+1)], na.rm = T)
      se.haz <- sd(avg.Haz[which(chang.vals == num.breaks),
                           1:(num.breaks+1)], na.rm = T)
    }else if(num.breaks == 1 ){

      mean.chng <- mean(avg.change[which(chang.vals == num.breaks),
                                   1:num.breaks], na.rm =T)
      se.chng <- sd(avg.change[which(chang.vals == num.breaks),1:num.breaks],na.rm = T)

      mean.haz <- colMeans(avg.Haz[which(chang.vals == num.breaks),
                                   1:(num.breaks+1)], na.rm = T)

      se.haz <- apply(avg.Haz[which(chang.vals == num.breaks),
                              1:(num.breaks+1)],2,FUN = sd, na.rm = T)


    }else{

      mean.chng <- colMeans(avg.change[which(chang.vals == num.breaks),
                                       1:num.breaks], na.rm =T)
      se.chng <- apply(avg.change[which(chang.vals == num.breaks),1:num.breaks]
                       ,2,FUN= sd,na.rm = T)

      mean.haz <- colMeans(avg.Haz[which(chang.vals == num.breaks),
                                   1:(num.breaks+1)], na.rm = T)

      se.haz <- apply(avg.Haz[which(chang.vals == num.breaks),
                              1:(num.breaks+1)],2,FUN = sd, na.rm = T)
    }

  }else{

    mean.chng<- NA
    se.chng <- NA
    mean.haz <- NA
    se.haz <- NA

  }

  return(list(final.chng = final.chng,prob.correct = prob.correct,
              chng = data.frame(mean.chng,se.chng),
              haz = data.frame(mean.haz,se.haz),
              prob.mod = chang.vals))
}


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
sim.path <- "C:\\Users\\phili\\OneDrive\\PhD\\R codes\\Changepoint Analysis\\Simulation Study\\Simulation Collapsing\\One Changepoint\\"
model.names <- c("decreasing_haz_lrg_","decreasing_haz_sml_","increasing_haz_lrg_","increasing_haz_sml_")
model.size <- c("200","500","1000")

#sim.path <- "C:/Users/phili/OneDrive/PhD/Simulation Study/Simulation Collapsing/"
#model.names <- c("bathtub_haz","decreasing_haz","increasing_haz","invertbath_haz")
#model.size <- c("300","500","1000")

model.com <- rep(NA, length(model.names)*length(model.size))
cnt <- 0
for(i in 1:length(model.names)){
  for(j in 1:length(model.size)){
    cnt <- cnt +1
    model.com[cnt] <- paste0(model.names[i],model.size[j])

  }
}
library(stringr)

n <- 1000

df.prob <- data.frame(probability= rep(NA,12),
                      model= rep(NA,12),
                      n =  rep(NA,12))

dict.names <- cbind(c("Increasing Small", "Increasing Large", "Decreasing Small", "Decreasing Large"),
                    c("increasing_haz_sml_","increasing_haz_lrg_", "decreasing_haz_sml_","decreasing_haz_lrg_"))

for(i in 1:length(model.com)){
  model.name <- model.com[i]
  sim.size <- as.numeric(gsub("[^0-9.]", "",  model.name))
  assign(model.name,loadRData(paste0(sim.path,model.com[i],".RData")))
  assign("current.model",loadRData(paste0(sim.path,model.com[i],".RData")))
  txt.string <- str_remove(model.name,paste0(sim.size))#str_remove(model.name,paste0(sim.size))#

  df.prob[i,1] <- current.model$prob.correct*100 #
  df.prob[i,2] <- dict.names[which(dict.names[,2]==txt.string),1]
  df.prob[i,3] <- sim.size

  for(j in 1:1){
    shape <- ((current.model$haz[j,1])^2)/(current.model$haz[j,2]^2)
    rate <- (current.model$haz[j,1])/(current.model$haz[j,2]^2)
    assign(paste0("df.res",j),data.frame(hazard = rgamma(n, shape = shape, rate = rate),
                                         changepoint = rnorm(n, current.model$chng[j,1], current.model$chng[j,2]),
                                         n = as.factor(sim.size), interval= as.factor(j), model = dict.names[which(dict.names[,2]==txt.string),1]))



  }
  shape <- ((current.model$haz[2,1])^2)/(current.model$haz[2,2]^2)
  rate <- (current.model$haz[2,1])/(current.model$haz[2,2]^2)

  df.res2 <-data.frame(hazard = rgamma(n, shape = shape, rate = rate),
                       changepoint = NA,
                       n = as.factor(sim.size), interval = as.factor(3),
                       model =  dict.names[which(dict.names[,2]==txt.string),1])


  assign(paste0("df.res_all",txt.string,sim.size),rbind(df.res1,df.res2))


  rm(df.res1,df.res2,txt.string)

}

for(i in 1:nrow(dict.names)){

  assign(paste0("df.res_final",dict.names[i,2]),rbind(get(paste0("df.res_all",dict.names[i,2],200)),
                                                      get(paste0("df.res_all",dict.names[i,2],500)),
                                                      get(paste0("df.res_all",dict.names[i,2],1000))))

}

df.prob$n <- as.factor(df.prob$n)
#df.prob$model <- as.factor(df.prob$model)

df.prob$model <- factor(df.prob$model, levels = c("Increasing Small", "Increasing Large", "Decreasing Small", "Decreasing Large"))
prob.plot <- ggplot(df.prob, aes(fill=n, y=probability, x=model)) +
  geom_bar(position="dodge", stat="identity")+
  scale_y_continuous(breaks = seq(0,100,10))+
  xlab("")+
  ylab("Probability")+
  theme_bw()+
  theme(legend.position = "none", plot.margin = unit(c(1.5,0.5,0.5,0.5), "lines"))
true.haz <- data.frame(model = dict.names[,2],
                       interval1  = c(0.25,0.75,0.75,0.2),
                       interval2  = c(0.5,0.5,0.2,0.75),
                       interval3  = c(0.75,0.25,0.75,0.2))


for(i in 1:nrow(dict.names)){

  df.plt.temp <- get(paste0("df.res_final",dict.names[i,2]))
  df.plt.temp$interval <- as.factor(df.plt.temp$interval)
  df.plt.temp$model <- as.factor(df.plt.temp$model)
  df.plt.temp$n <- as.factor(df.plt.temp$n)

  #index <- which(dict.names[i,2]==true.haz[,1])
  #print(index)
  gg.plt.temp <-  ggplot(df.plt.temp, aes(x=interval, y=hazard, fill=n)) +
    geom_boxplot(position=position_dodge(1), outlier.shape = NA)+
    #geom_segment(aes(x=2/3,xend=4/3,y=true.haz[index,1+1],yend=true.haz[index,1+1]), colour = "black", linetype = 1, size = 1, show.legend   = FALSE)+
    #geom_segment(aes(x=5/3,xend=7/3,y=true.haz[index,2+1],yend=true.haz[index,2+1]), colour = "orange", linetype = 1, size = 1, show.legend   = FALSE) +
    #geom_segment(aes(x=8/3,xend=10/3,y=true.haz[index,3+1],yend=true.haz[index,3+1]), colour = "purple", linetype = 1, size = 1, show.legend   = FALSE)+
    xlab("")+
    scale_y_continuous(breaks = seq(0,1.2,.2))+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5, linetype = "solid"))+
    coord_cartesian(ylim = c(0,1.2),expand = FALSE)
  #print(c(true.haz[index,1+1],true.haz[index,2+1],true.haz[index,3+1]))

  gg.plt.temp
  assign(paste0("p_hazard_",dict.names[i,2]),gg.plt.temp)
  rm(gg.plt.temp)

  gg.plt.temp2  <- ggplot(df.plt.temp[!is.na(df.plt.temp$changepoint),], aes(x=interval, y=changepoint, fill=n)) +
    geom_boxplot(position=position_dodge(1), outlier.shape = NA)+
    #geom_segment(aes(x=2/3,xend=4/3,y=0.5,yend=0.5), colour = "black", linetype = 1, size = 1, show.legend  = FALSE)+
    #geom_segment(aes(x=5/3,xend=7/3,y=1,yend=1), colour = "orange", linetype = 1, size = 1, show.legend   = FALSE)+
    xlab("")+
    ylab("")+
    scale_y_continuous(breaks = seq(.2,1.4,.2))+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    coord_cartesian(ylim = c(.2,1.4),expand = FALSE)

  assign(paste0("p_change_",dict.names[i,2]),gg.plt.temp2)
  rm(gg.plt.temp2)

}

plt.output <- ggarrange(
  ggarrange(p_hazard_increasing_haz_sml_ +ylab("Increasing Small"),
            p_change_increasing_haz_sml_,
            p_hazard_increasing_haz_lrg_ +ylab("Increasing Large"),
            p_change_increasing_haz_lrg_,
            p_hazard_decreasing_haz_sml_ +ylab("Decreasing Small"),
            p_change_decreasing_haz_sml_,
            p_hazard_decreasing_haz_lrg_ +ylab("Decreasing Large"),
            p_change_decreasing_haz_lrg_,
            nrow = 4, ncol = 2, common.legend = TRUE, legend = "bottom"), prob.plot+theme(axis.text=element_text(size=8)))
#Could annotate figure at bottom

ggsave(annotate_figure(plt.output, fig.lab.pos = c("top.left"),fig.lab.face= "bold",fig.lab =  c("Hazard across intervals         Changepoints                           Probability of selecting correct model")) ,
       filename = "Collapsing Model One Changepoint.png", width = 8.5, height =6)



# 8. Read in Simulation Study Results for Publication -----


#Simulation
path <-"C:\\Users\\phili\\OneDrive\\PhD\\R codes\\Changepoint Analysis\\Biometrics publication\\Simulation\\"


rate1 <- c(0.5,0.75)
rate2 <- c(0.25,0.75)
rate3 <- c(0.75,0.5)
rate4 <- c(0.75,0.25)

rate5 <- c(0.25,0.5,0.75)
rate6 <- c(0.75,0.5,0.25)
rate7 <- c(0.75,0.2,0.75)
rate8 <- c(0.2,0.75,0.2)

time1 <- c(0.5,1)

censor_vec <- c(0,0.33)
sample_vec <- c(300,500,1000)
rate_list <- list()

for(i in 1:8){
  current_rate <-get(paste0("rate",i))
  rate_list[[i]] <- current_rate
}

res_vec <- array(NA, dim = c(length(censor_vec)*length(sample_vec)*length(rate_list),2))
res_vec <- array(NA, dim = c(1*length(sample_vec)*length(rate_list),7))
res_vec <- data.frame(res_vec)
f <- 1
for(i in 1:length(rate_list)){
  for(j in 1:length(sample_vec)){
    for(q in 1:1 ){
      #for(q in 1:length(censor_vec) ){
      current_name <- paste0("Gibbs_rate_",paste0(rate_list[[i]],collapse = "_"),"_size_",sample_vec[j],"_censor_",censor_vec[q])
      res_vec[f,1] <- current_name
      res <- readRDS(file = paste0(path,"Gibbs Simulations 2021\\",current_name,".rds"))
      res_vec[f,2]<- as.numeric(res$prob.correct)
      res_vec[f,3] <- as.numeric(res$chng[1,1])
      res_vec[f,4] <- as.numeric(res$chng[1,2])
      if(nrow(res$chng)>1){
        res_vec[f,5] <- as.numeric(res$chng[2,1])
        res_vec[f,6] <- as.numeric(res$chng[2,2])
      }
      res_vec[f,7] <- i
      f <- f +1
    }
  }
}

library(dplyr)
colnames(res_vec) <- c("Name","Prob","Tau1","SE1","Tau2","SE2","Model")


res_vec <- as_tibble(res_vec) %>% mutate_if(is.numeric,round, digits = 2)
res_vec$Prob <- res_vec$Prob*100
res_vec$Comb1 <- paste0(res_vec$Tau1," (",res_vec$SE1,")")
result <- aggregate(Comb1 ~ Model, data = res_vec, paste, collapse = " & ")
result$Comb1 <- paste0(result$Comb1," & 0.5 ")
res_vec$Comb2 <- paste0(res_vec$Tau2," (",res_vec$SE2,")")
result2 <- aggregate(Comb2 ~ Model, data = res_vec, paste, collapse = " & ")
result2$Comb2 <- paste0(result2$Comb2," & 1.0 ")
result3 <- aggregate(Prob ~ Model, data = res_vec, paste, collapse = " & ")


#No Changepoint

rate_list <- c(0.25, 0.5,0.75)
sample_vec <- c(100,200)
censor_vec <- c(0,0.5)
res_vec2 <- res_vec <- data.frame(array(NA, dim = c(1*length(sample_vec)*length(rate_list),2)))
f <- 1
for(i in 1:length(rate_list)){
  for(j in 1:length(sample_vec)){
    for(q in 1:length(censor_vec) ){

      current_name <- paste0("Gibbs_rate_",paste0(rate_list[[i]],collapse = "_"),"_size_",sample_vec[j],"_censor_",censor_vec[q])
      current_name2 <- paste0("Collapsing_rate_",paste0(rate_list[[i]],collapse = "_"),"_size_",sample_vec[j],"_censor_",censor_vec[q])
      res <- readRDS(file = paste0(path,"No Changepoint//",current_name,".rds"))

      res_vec[f,1] <-current_name
      res_vec[f,2] <-res$prob.correct

      res2 <- readRDS(file = paste0(path,"No Changepoint//",current_name2,".rds"))
      res_vec2[f,1] <-current_name2
      res_vec2[f,2] <-mean(res2$prob.mod == 0)

      f <- f +1
    }
  }
}




res_vec$X2 <- round(res_vec$X2*100,0)
res_vec2$X2 <- round(res_vec2$X2*100,0)
res_vec <- res_vec[-grep("size_200_censor_0$", res_vec$X1),]
res_vec2 <- res_vec2[-grep("size_200_censor_0$", res_vec2$X1),]

cbind(paste0(res_vec$X2," & ",res_vec2$X2))




#Simulation Read in Collapsing


path <-"C:\\Users\\phili\\OneDrive\\PhD\\R codes\\Changepoint Analysis\\Biometrics publication\\Simulation\\"


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


rate1 <- c(0.5,0.75)
rate2 <- c(0.25,0.75)
rate3 <- c(0.75,0.5)
rate4 <- c(0.75,0.25)

rate5 <- c(0.25,0.5,0.75)
rate6 <- c(0.75,0.5,0.25)
rate7 <- c(0.75,0.2,0.75)
rate8 <- c(0.2,0.75,0.2)

time1 <- c(0.5,1)

censor_vec <- c(0,0.33)
sample_vec <- c(300,500,1000)
rate_list <- list()

for(i in 1:8){
  current_rate <-get(paste0("rate",i))
  rate_list[[i]] <- current_rate
}

res_vec <- array(NA, dim = c(length(censor_vec)*length(sample_vec)*length(rate_list),2))
res_vec <- array(NA, dim = c(1*length(sample_vec)*length(rate_list),7))
res_vec <- data.frame(res_vec)
f <- 1
for(i in 1:length(rate_list)){
  for(j in 1:length(sample_vec)){
    for(q in 1:1 ){
      #for(q in 1:length(censor_vec) ){
      current_name <- paste0("Collapsing_rate_",paste0(rate_list[[i]],collapse = "_"),"_size_",sample_vec[j],"_censor_",censor_vec[q])
      res_vec[f,1] <- current_name

      #if(j ==3 ){
      #res <- loadRData(fileName = paste0(path,"Collapsing Simulations without 1000\\",current_name,".RData"))
      #}else{
      #res <- readRDS(file = paste0(path,"Collapsing Simulations without 1000\\",current_name,".rds"))
      #}
      res <- readRDS(file = paste0(path,"Collapsing Simulations\\",current_name,".rds"))


      res_vec[f,2]<- as.numeric(res$prob.correct)
      res_vec[f,3] <- as.numeric(res$chng[1,1])
      res_vec[f,4] <- as.numeric(res$chng[1,2])
      if(nrow(res$chng)>1){
        res_vec[f,5] <- as.numeric(res$chng[2,1])
        res_vec[f,6] <- as.numeric(res$chng[2,2])
      }
      res_vec[f,7] <- i
      f <- f +1
    }
  }
}




colnames(res_vec) <- c("Name","Prob","Tau1","SE1","Tau2","SE2","Model")


res_vec <- as_tibble(res_vec) %>% mutate_if(is.numeric,round, digits = 2)
res_vec$Prob <- res_vec$Prob*100
res_vec$Comb1 <- paste0(res_vec$Tau1," (",res_vec$SE1,")")
result <- aggregate(Comb1 ~ Model, data = res_vec, paste, collapse = " & ")
#result$Comb1 <- paste0(result$Comb1," & 0.5 ")
res_vec$Comb2 <- paste0(res_vec$Tau2," (",res_vec$SE2,")")
result2 <- aggregate(Comb2 ~ Model, data = res_vec, paste, collapse = " & ")
#result2$Comb2 <- paste0(result2$Comb2," & 1.0 ")
result3 <- aggregate(Prob ~ Model, data = res_vec, paste, collapse = " & ")


#Or just write a latex table




