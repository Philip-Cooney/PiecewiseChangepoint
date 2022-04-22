
library("PiecewiseChangepoint")


set.seed(123)
n_obs =100
n_events_req=100
max_time =  2

rate = c(0.75,0.25)
t_change =1

df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                   num.breaks = length(t_change),rate = rate ,
                   t_change = t_change, max_time = max_time)


survfit <- survfit(Surv(time,status)~1, data = df)

summary(survfit)

png("Cumulative hazard.png")
plot(y = summary(survfit)[["cumhaz"]], x =  summary(survfit)[["time"]],
     xlab = "time", ylab = "Cumulative Hazard")
dev.off()

Collapsing_Model <- collapsing.model(df,
                                     n.iter = 20750,
                                     burn_in = 750,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1)



Collapsing_Model_uniform <- collapsing.model(df,
                                     n.iter = 20750,
                                     burn_in = 750,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1,
                                     MLE =TRUE)


library(progress)
Gibbs_sampler <- gibbs.sampler(df,
                          num.breaks =1, # Number of Changepoints
                          n.iter =10000,
                          burn_in = 100,
                          num.chains = 2,
                          n.thin = 1,
                          alpha.hyper = 1, #Hyperparameter for Marginal Likelihood
                          beta.hyper = 1, #Hyperparameter for Marginal Likelihood
                          MLE = FALSE) #If TRUE, considers on Uniform Prior){

colMeans(Gibbs_sampler$changepoint)

Gibbs_sampler_uniform <- gibbs.sampler(df,
                               num.breaks =1, # Number of Changepoints
                               n.iter =10000,
                               burn_in = 100,
                               num.chains = 2,
                               n.thin = 1,
                               alpha.hyper = 1, #Hyperparameter for Marginal Likelihood
                               beta.hyper = 1, #Hyperparameter for Marginal Likelihood
                               MLE = T) #If TRUE, considers on Uniform Prior){

colMeans(Gibbs_sampler_uniform$changepoint)

