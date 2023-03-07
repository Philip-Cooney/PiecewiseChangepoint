
devtools::load_all()
set.seed(125)
n_obs =500
n_events_req=500
max_time =  2

rate = c(0.75,0.25)
t_change =1

df <- PiecewiseChangepoint::gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                   num.breaks = length(t_change),rate = rate ,
                   t_change = t_change, max_time = max_time)
debug(PiecewiseChangepoint::rpwexp)
library(survival)
plot(survfit(Surv(time,status)~1, data = df))

sum(df$time[which(df$time < 1)])/sum(df$status[which(df$time < 1 & df$status ==1)])
recasted_df <- PiecewiseChangepoint::df_recast(df)


sum(recasted_df[which(as.numeric(rownames(recasted_df)) < 1), 2])/sum(recasted_df[which(as.numeric(rownames(recasted_df)) < 1), "time_diff_event"])

sum(recasted_df[which(as.numeric(rownames(recasted_df)) >1), 2])/sum(recasted_df[which(as.numeric(rownames(recasted_df)) > 1), "time_diff_event"])


Collapsing_Model <- PiecewiseChangepoint::collapsing.model(df,
                                     n.iter = 10000,
                                     burn_in = 5,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1)
PiecewiseChangepoint:::rpwexp
plot(Collapsing_Model)

St_df <- get_Surv(Collapsing_Model, max_predict = 10)
str(St_df)

St_mu <- rowMeans(St_df)
as.numeric(names(St_mu))

plot(survfit(Surv(df$time_event,rep(1,nrow(df)))~1))
lines(y = St_mu, x = as.numeric(names(St_mu)))
lines(y = data.frame(psa$S[[1]])[,"S"], x = data.frame(psa$S[[1]])[,"t"])





data.frame(psa$S[[1]])[,"S"]
m5 <- fit.models(Surv(time, status) ~ 1, data = df, distr = "rps", k = 2, method = "hmc", iter = 10000, warmup = 1000)

plot(m5, add.km = T)

psa <- make.surv(fit = m5, nsim = 1000, t = seq(0, 10, by = 0.01))
psa.plot(psa, add_km = T)


mod_comp <-compare.surv.mods(Collapsing_Model)

mod_comp$mod.comp

add_km(mod_comp$plot_Surv_all, data.frame(time = df$time_event, status = 1), colour = "black")


mod_comp$plot_Surv_all


library(flexsurv)
sp1 <- flexsurvspline(Surv(time, status) ~ 1, data = df, k = 1,
                       scale = "hazard")

data_jags_spline <- list()
k <- 2
knots <- quantile(log(df$time[df$status == 1]),
                  seq(0, 1, length = k + 2))
B <- flexsurv::basis(knots, log(df$time))
DB <- flexsurv::dbasis(knots, log(df$time))
max_predict <- 10
interval <- 0.5
data_jags_spline$t_pred <- seq(0.0001, max_predict, by = interval)

B_pred <- flexsurv::basis(knots, log(data_jags_spline$t_pred))
data_jags_spline$B_pred <- B_pred
data_jags_spline$B <- B
data_jags_spline$DB <- DB
data_jags_spline$t <- df$time
data_jags_spline$N <- nrow(df)
data_jags_spline$K <- 2
data_jags_spline$status <- df$status

flex_spline <- flexsurvspline(Surv(time, status) ~ 1, data = df, k = k,
                      scale = "hazard")

data_jags_spline$gamma_prec_mat <- solve(flex_spline$cov)
data_jags_spline$gamma_mu_vec <- flex_spline$res[,1]



spline.mod <- "
data{
  for(i in 1:N){
    zeros[i] <- 0
  }
}

model{
  C <- 1000

  eta <- B%*%gamma
  eta_prime <- DB%*%gamma
  for(i in 1:N){
    log_lik[i] = status[i]*(-log(t[i]) + log(eta_prime[i]) + eta[i]) - exp(eta[i])
    Like[i] <- exp(log_lik[i])
    zeros[i] ~ dpois(-log_lik[i]+C)
  }


# Poor convergence; supply priors from MLE
  #for(i in 1:(K+2)){
  #  gamma[i] ~ dnorm(0,0.001)
  #}

  gamma ~ dmt(gamma_mu_vec, gamma_prec_mat,3)


  eta_pred = B_pred%*%gamma;

  for(i in 1:length(t_pred)){
    St_pred[i] = exp(-exp(eta_pred[i]));
  }

}
"


rcs.mod <- R2jags::jags(
  model.file = textConnection(spline.mod),
  data = data_jags_spline,
  #inits = inits_list("weibull", n.chains),
  n.chains = 2,
  n.iter = 100000,
  n.thin = 10,
  parameters.to.save = c("gamma")
  #,
  #inits = inits_spline()
  )
library(ggmcmc)
jagsfit.mcmc <- as.mcmc(rcs.mod)
S <- ggs(jagsfit.mcmc)
ggmcmc(S, file="model_simple-diag.pdf")


