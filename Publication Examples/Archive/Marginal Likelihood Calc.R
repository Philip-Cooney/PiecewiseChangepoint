#library(devtools)

#devtools::load_all()

library("PiecewiseChangepoint")


# set.seed(123)
# n_obs =100
# n_events_req=100
# max_time =  2
#
# rate = c(0.75,0.25)
# t_change =1
#
# df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
#                    num.breaks = length(t_change),rate = rate ,
#                    t_change = t_change, max_time = max_time)


library("xlsx")

leukemia.df <- read.xlsx(paste0("C:\\Users\\phili\\OneDrive\\PhD\\R Codes\\Changepoint Analysis\\","Achar Leukemia Data.xlsx"),1)
df <- leukemia.df <- leukemia.df[which(leukemia.df[,1] !=182),]
#gibbs<- PiecewiseChangepoint::gibbs.sampler(df, 1, 1000, num.chains = 2)

Collapsing_Model <- collapsing.model(df,
                                     n.iter = 50000,
                                     burn_in = 750,
                                     n.chains = 3,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1/365)
#remove.packages("PiecewiseChangepoint")

library("Brobdingnag")

library("combinat")

# generating combinations of the
# alphabets taking 2 at a time

beta.hyper <- 365
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

# exp(log_marg[3]-log_marg[2]-log(2))
# exp(log_marg[2]-log_marg[1])
#
# exp(log_marg[4]-log_marg[3]- (dpois(c(2),1,log = T)-dpois(c(3),1,log = T)) )
# exp(-1.342470)
#
# 39.6/17.2
#
# dpois(c(2),1,log = T)-  dpois(c(3),1,log = T)
# dpois(c(2),1,log = F)/ dpois(c(3),1,log = F)
#
# 34.8/1.6
#
# log(marg.like[[1]])
#
# unlist(marg.like)
# exp(log(marg_vec2[[2]]-marg_vec2[1]-log(2))
# sum(df$status)/sum(df$time)
#
#
#
# dpois(1,lambda = 1)/dpois(2,lambda = 1)
#
#
# apply(combs[,k_prior_loc,drop = F], 2,
#                                           function(x){PiecewiseChangepoint::margin_lik_calc_log(time_diffs,x,1,365)+
#                                               log(PiecewiseChangepoint:::pos_prob(x, n_events, FALSE))})
#
# log(PiecewiseChangepoint:::pos_prob(combs[,1], n_events, FALSE))
#
# apply(combs[,k_prior_loc,drop = F], 2,
#       function(x){PiecewiseChangepoint::margin_lik_calc_log(time_diffs,x,1,365)})
#
# df2 <- df
# df2$time <- df2$time/365
# time_diffs2 <- df_recast(df2)
# margin_lik_calc_log(time_diffs2,NA,1,1)
#
# log(exp(margin_lik_calc_log(time_diffs2,NA,1,1)), base = 10)
#
# library("Brobdingnag")
# as.brob(exp(margin_lik_calc_log(df_recast(leukemia.df),NA, alpha = 1,beta =365)))
#
#
# # library(abind)
# #
# # res <- abind::abind(gibbs$k,gibbs$lambda.array, along = 2)
# # mcmc.list_res <- list()
# # for(c in 1:dim(res)[3]){
# #   mcmc.list_res[[paste0("Chain ",c)]] <- as.mcmc(res[,,c])
# #
# # }
# # as.mcmc.list(mcmc.list_res)
#
#
# Collapsing_Model <- collapsing.model(df,
#                                      n.iter = 20750,
#                                      burn_in = 750,
#                                      n.chains = 2,
#                                      alpha.hyper = 1,
#                                      beta.hyper1 = 1,
#                                      beta.hyper2 = 1/365)


# num.breaks = 1 # Number of Changepoints
# n.iter =500
# burn_in = 0
# num.chains = 2
# n.thin = 1
# alpha.hyper = 1 #Hyperparameter for Marginal Likelihood
# beta.hyper = 1 #Hyperparameter for Marginal Likelihood
# MLE = FALSE
#
#
# View(Collapsing_Model)
