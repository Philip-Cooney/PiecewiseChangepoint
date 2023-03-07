#path <- "C:/Users/phili/OneDrive/PhD/R codes/Changepoint Analysis/"
path <- "C:/Users/phili/OneDrive/PhD/R codes/Changepoint Analysis/GitHub/"


path_biometrics <-"C:\\Users\\phili\\OneDrive\\PhD\\R codes\\Changepoint Analysis\\Biometrics publication\\"

source(paste0(path, "Functions Gibbs Sampler.R"))
source(paste0(path, "Functions Collapsing.R"))


#Change in Changed December has no impact on the results
#source(paste0(path, "Functions Gibbs Markdown changed December.R"))
#source(paste0(path, "Functions Reversible Jump chains final hyper.R"))

#Leukemia dataset 
leukemia.df <- read.xlsx(paste0("C:\\Users\\phili\\OneDrive\\PhD\\R Codes\\Changepoint Analysis\\","Achar Leukemia Data.xlsx"),1)
leukemia.df <- leukemia.df[which(leukemia.df[,1] !=182),]

#res_final1 <- gibbs.sampler(leukemia.df, 1, 10100, 100,1,alpha.hyper = 1, beta.hyper = 365, max_predict = 2000, interval = 50)
#res_final2 <- gibbs.sampler(leukemia.df, 2, 10100, 1000,1,alpha.hyper = 1, beta.hyper = 365, max_predict = 2000, interval = 50)
#res_final3 <- gibbs.sampler(leukemia.df, 3, 10100, 1000,1,alpha.hyper = 1, beta.hyper = 365,max_predict = 5*365, interval = 50)
#res_final4 <- gibbs.sampler(leukemia.df, 4, 10100, 1000,1,alpha.hyper = 1, beta.hyper = 365,max_predict = 5*365, interval = 50)
#res_final5 <- gibbs.sampler(leukemia.df, 5, 10100, 1000,1,alpha.hyper = 1, beta.hyper = 365,max_predict = 5*365, interval = 50)
#res_final6 <- gibbs.sampler(leukemia.df, 6, 10100, 1000,1,alpha.hyper = 1, beta.hyper = 365,max_predict = 5*365, interval = 50)

#saveRDS(res_final1, paste0(path_biometrics,"Leukemia Example\\Gibbs 1.rds"))
#saveRDS(res_final2, paste0(path_biometrics,"Leukemia Example\\Gibbs 2.rds"))
#saveRDS(res_final3, paste0(path_biometrics,"Leukemia Example\\Gibbs 3.rds"))
#saveRDS(res_final4, paste0(path_biometrics,"Leukemia Example\\Gibbs 4.rds"))
#saveRDS(res_final5, paste0(path_biometrics,"Leukemia Example\\Gibbs 5.rds"))
#saveRDS(res_final6, paste0(path_biometrics,"Leukemia Example\\Gibbs 6.rds"))


# marg_res_days <- c(as.brob(exp(margin_lik_calc_log(df_recast(leukemia.df),NA, alpha = 1,beta =365))),
#                     res_final1[["mean.marg.lik.log"]],
#                     res_final2[["mean.marg.lik.log"]],
#                     res_final3[["mean.marg.lik.log"]],
#                     res_final4[["mean.marg.lik.log"]],
#                     res_final5[["mean.marg.lik.log"]],
#                    res_final6[["mean.marg.lik.log"]])


# gibbs.model.selc(marg_res_days)
# 
# marg_res_days <- exp(sapply(marg_res_days, log))#*dpois(0:5,lambda = 1)
# marg_res_days2 <- sapply(marg_res_days, log10)#*dpois(0:5,lambda = 1)
# round(marg_res_days2,2)

# round(marg_res_days[-1]/marg_res_days[-length(marg_res_days)],2)

#I tweaked the RJMCMC_Cpp_function_v5.cpp to set the prior to be uniform

leuk.RJMC2 <- collapsing.model(n.iter = 20750, df = leukemia.df,n.chains = 1,  burn_in = 750,
                          alpha.hyper = 1, beta.hyper1 = 1, beta.hyper2 = 1/365)

leuk.RJMC2$prob.changepoint
 # leuk.RJMC2 <- collapsing.model(n.iter = 20750, df = leukemia.df,n.chains = 1,  burn_in = 750,
 #                          alpha.hyper = 1, beta.hyper1 = 3650, beta.hyper2 = 10)
 # 
 # leuk.RJMC2$prob.changepoint
 # 
leuk.RJMC2[["changepoint"]]

# No issues with different timescales
 # time_diffs <- df_recast(leukemia.df)
 # margin_lik_calc_log(time_diffs, 10, 1, 365)- margin_lik_calc_log(time_diffs, NA, 1, 365)
 # 
  leukemia.df2 <- leukemia.df
  leukemia.df2$time <- leukemia.df$time/365
 # time_diffs2 <- df_recast(leukemia.df2)
 # margin_lik_calc_log(time_diffs2, 10, 1, 1)- margin_lik_calc_log(time_diffs2, NA, 1, 1)
 # 
 # 
 # 
  leuk.RJMC2_time <- collapsing.model(n.iter = 20750, df = leukemia.df2,n.chains = 1,  burn_in = 750,
                                 alpha.hyper = 1, beta.hyper1 = 1, beta.hyper2 = 1)
 # 
  leuk.RJMC2_time$prob.changepoint
 # 
 
leuk.RJMC2 <- readRDS(paste0(path_biometrics,"Leukemia Example\\Collapsing.rds"))

round(leuk.RJMC2[["prob.changepoint"]],3)*100
summary(leuk.RJMC2[["beta.array"]])

options(scipen=999)

leuk.RJMC2[["prob.changepoint"]][-1]/leuk.RJMC2[["prob.changepoint"]][-length(leuk.RJMC2[["prob.changepoint"]])]
round(leuk.RJMC2[["prob.changepoint"]],3)*100


# saveRDS(leuk.RJMC2, paste0(path_biometrics,"Leukemia Example\\Collapsing.rds"))


#Compare the changepoint models

#Collapsing
chang_num_vec <- apply(leuk.RJMC2[["k.stacked"]],1,function(x){length(na.omit(x))})
table(chang_num_vec)

changepoint_num <- 1
#summary(leuk.RJMC2[["changepoint"]][which(chang_num_vec == changepoint_num),])
#summary(leuk.RJMC2[["lambda"]][which(chang_num_vec == changepoint_num),])

Mode(leuk.RJMC2[["changepoint"]][which(chang_num_vec == changepoint_num),1])

apply(leuk.RJMC2[["changepoint"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)
apply(leuk.RJMC2[["changepoint"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)
apply(leuk.RJMC2[["changepoint"]][which(chang_num_vec == changepoint_num),],2,quantile,na.rm =T, probs = c(0.025,0.5,0.975))

apply(leuk.RJMC2[["lambda"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)
apply(leuk.RJMC2[["lambda"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)
apply(leuk.RJMC2[["lambda"]][which(chang_num_vec == changepoint_num),],2,quantile,na.rm =T, probs = c(0.025,0.5,0.975))

apply(leuk.RJMC2[["lambda.df_mu"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)
apply(leuk.RJMC2[["lambda.df_mu"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)
apply(leuk.RJMC2[["lambda.df_mu"]][which(chang_num_vec == changepoint_num),],2,quantile,na.rm =T, probs = c(0.025,0.5,0.975))

apply(sqrt(leuk.RJMC2[["lambda.df_var"]][which(chang_num_vec == changepoint_num),]),2,mean,na.rm =T)
apply(leuk.RJMC2[["lambda"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)


#Gibbs

res_final1<- readRDS(paste0(path_biometrics,"Leukemia Example\\Gibbs 1.rds"))

summary(res_final1[["stacked_df"]])

apply(res_final1[["stacked_df"]],2,mean,na.rm =T)
apply(res_final1[["stacked_df"]],2,sd,na.rm =T)
apply(res_final1[["stacked_df"]],2,quantile,na.rm =T, probs = c(0.025,0.5,0.975))

Lambda_column_mu <- c(apply(res_final1[["stacked_df"]],2,mean,na.rm =T)[2:3], #Gibbs Sampler
                      na.omit(apply(leuk.RJMC2[["lambda"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)),#Uncollapsed Collapsing model
                      na.omit(apply(leuk.RJMC2[["lambda.df_mu"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)) #Collapsing Post-hoc Expect
                      )
index_1 <- seq(from = 1, to =length(Lambda_column_mu), by = 2)
index_2 <- seq(from = 2, to =length(Lambda_column_mu), by = 2)

Lambda_column_mu1 <- Lambda_column_mu[index_1]
Lambda_column_mu2 <- Lambda_column_mu[index_2]

Lambda_column_sd <- c(apply(res_final1[["stacked_df"]],2,sd,na.rm =T)[2:3], #Gibbs Sampler
                      na.omit(apply(leuk.RJMC2[["lambda"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)),#Uncollapsed Collapsing model
                      na.omit(apply(sqrt(leuk.RJMC2[["lambda.df_var"]][which(chang_num_vec == changepoint_num),]),2,mean,na.rm =T)) #Collapsing Post-hoc Expect
)

Lambda_column_sd1 <- Lambda_column_sd[index_1]
Lambda_column_sd2 <- Lambda_column_sd[index_2]


Changepoint_column <- c(apply(res_final1[["stacked_df"]],2,mean,na.rm =T)[1],
                        na.omit(apply(leuk.RJMC2[["changepoint"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)),
                        NA)


Changepoint_sd <- c(apply(res_final1[["stacked_df"]],2,sd,na.rm =T)[1],
                        na.omit(apply(leuk.RJMC2[["changepoint"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)),
                        NA)

results.vec_leuk <- cbind(Lambda_column_mu1, Lambda_column_sd1, Lambda_column_mu2, Lambda_column_sd2,Changepoint_column,Changepoint_sd)

results.vec_leuk <- apply(results.vec_leuk,2, round,digits = 5)


rownames(results.vec_leuk) <- c("Gibbs", "Collapsing - Uncollapsed", "Collapsing - Post Hoc")

colnames(results.vec_leuk) <- c("Lambda 1", "sd", "Lambda 2", "sd", "Changepoint", "sd")

#Flatter priors can be obtainted increasing beta.hyper 2 keeping beta.hyper1 in the same proportion




path_latex <- "C:/Users/phili/OneDrive/PhD/Latex Folder/"
res_final1[["plot_Surv"]]
ggsave(paste0(path_latex,"Leukemia_1_change.png"),width = 6,height = 4)
res_final2[["plot_Surv"]]
ggsave(paste0(path_latex,"Leukemia_2_change.png"),width=6,height =4)
res_final2[["plot_Surv"]]
res_final2[[plot]]


ggarrange(res_final1[["plot_Surv"]]+xlab("Time"),
          res_final2[["plot_Surv"]]+xlab("Time")+ylab(""))
ggsave(paste0(path_latex,"Leukemia_both.png"),width=12,height =4)

round(t(apply(res_final1[["stacked_df"]],2,quantile, probs = c(0.025,0.5,0.975))),5)
round(t(apply(res_final1[["stacked_df"]],2,quantile, probs = c(0.025,0.5,0.975)))[-1,]*365.25,3)


round(t(apply(res_final1[["stacked_df"]],2,quantile, probs = c(0.05,0.5,0.95))))
num.changepoints <- apply(leuk.RJMC2[["changepoint"]],1,function(x){length(na.omit(x))})
quantile(leuk.RJMC2[["changepoint"]][which(num.changepoints==1),1], probs = c(0.05,0.5,0.95))


round(t(apply(res_final2[["stacked_df"]],2,quantile, probs = c(0.025,0.5,0.975))),5)

apply(leuk.RJMC2[["changepoint"]][which(num.changepoints==2),1:2],2,quantile, probs = c(0.025,0.5,0.975))

#round(t(apply(res_final5[["stacked_df"]],2,quantile, probs = c(0.025,0.5,0.975))),5)


