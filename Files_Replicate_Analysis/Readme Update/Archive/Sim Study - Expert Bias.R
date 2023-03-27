
list.of.packages <- need<-c("overlapping", "tidyr", "expertsurv", "xlsx", "dplyr","flexsurv") #needed libraries

res <- lapply(list.of.packages, require, character.only = TRUE)

not.loaded <-   list.of.packages[which(sapply(res, unlist) ==F)]
not.loaded2 <- lapply(not.loaded, require, character.only = TRUE)
not.installed <-   not.loaded2[which(sapply(not.loaded2, unlist) ==F)]
#load the packages
if(length(not.installed)) install.packages(not.loaded)
if(length(not.installed)) lapply(not.loaded, require, character.only = TRUE)

path_expertsurv <- "~/expertsurv - Sim Study/"
source(paste0(path_expertsurv,"Flexsurv functions v2.R"))


shape_vec <- c(0.75, 1, 1.25)
scale_vec <- seq(0.1,1, by = 0.2)
n_vec <- c(30,50,100)#samp_vec[samp]
cut_off <- c(2)
expert_time_1 <- c(4,NA)
expert_time_2 <- c(10,NA)
expert_mean1 <- c(0.1, 0.3)
expert_mean2 <- c(0.1, 0.05)
sd <- c(0.1,0.05,0.025)
prop.correct <- c(.75,1,1.25)
n.knots <- 1
beta_dist =F
time_horizon <- 15
time_seq <- seq(0,time_horizon, by = 0.1)


expert_time_all <- as.numeric(na.omit(unique(c(expert_time_1,expert_time_2))))
n.boots <- 500

path_output <- gen_pathway(paste0(path_expertsurv,"Bias Expert/"))


l <- list(param_1 = shape_vec, param_2 = scale_vec, samp_size = n_vec, cut_off = cut_off,
          expert_time_1 = expert_time_1,expert_time_2 =expert_time_2, 
          sd = sd, prop.correct= prop.correct)
sim_comb <- expand.grid(l)

surv_vals <- t(apply(sim_comb, 1, function(x){flexsurv::pweibullPH(na.omit(c(expert_time_1, expert_time_2)), 
                                                                   shape = x["param_1"], scale = x["param_2"],
                                                                   lower.tail = F)}))
surv_vals <- data.frame(surv_vals)
names(surv_vals) <- c("true_surv1","true_surv2" )

sim_comb2 <- cbind(sim_comb,surv_vals)

sim_comb2 <- sim_comb2 %>% data.frame() %>% filter(true_surv2 > 0.05, true_surv2 < 0.5, !is.na(expert_time_1) | !is.na(expert_time_2))



index <- which(is.na(sim_comb2$expert_time_1) & sim_comb2$expert_time_2 == 10&
                 sim_comb2$param_1==0.75 & sim_comb2$param_2 == 0.5 & sim_comb2$samp_size ==30 & sim_comb2$cut_off  ==2&sim_comb2$sd == 0.1)
# sim_comb2[index,]
# i <- index
index <- sum(which(is.na(sim_comb2$expert_time_1) & sim_comb2$expert_time_2 == 10) < 135)



sim_comb2$expert_mean1 <- round(sim_comb2$true_surv1*sim_comb2$prop.correct,2)
sim_comb2$expert_mean2 <- round(sim_comb2$true_surv2*sim_comb2$prop.correct,2)

mod_all <- c("exp","weibullPH","lnorm","gamma","gompertz","llogis","lognormal","rps") 
mod_selc <- "weibullPH"
RMSE_res_mat <-ISSE_res_mat <- matrix(nrow = nrow(sim_comb2), ncol = 22)
colnames(RMSE_res_mat) <-colnames(ISSE_res_mat) <- c(paste0("Expert ", c("2.5%", "50%", "97.5%")),
                                                     paste0("MLE ", c("2.5%", "50%", "97.5%")),
                                                     "Mean Diff", "SD diff", "N",
                                                     paste0("RSME Ratio ", c("2.5%", "50%", "97.5%", "SD")),
                                                     paste0("Expert - RSt Ratio", c("2.5%", "50%", "97.5%")),
                                                     paste0("MLE - RSt Ratio ", c("2.5%", "50%", "97.5%")),
                                                     paste0("RSt - Expert/MLE Ratio ", c("2.5%", "50%", "97.5%")))

ISSE_res_mat <- ISSE_res_mat[,1:13]
for(i in 1:nrow(sim_comb2)){
  print(paste0("Sim ",i))
  cut_off_sim <- sim_comb2[i, "cut_off"]
  samp_size_sim <-  sim_comb2[i, "samp_size"]
  param_1_sim  <-   sim_comb2[i, "param_1"] 
  param_2_sim  <-   sim_comb2[i, "param_2"] 
  expert_mult  <-   sim_comb2[i, "prop.correct"] 
  
  True_RSt <- flexsurv:::rmst_weibullPH(time_horizon,shape = param_1_sim, scale = param_2_sim)
  True_St <- flexsurv:::pweibullPH(time_seq,shape = param_1_sim, scale = param_2_sim,lower.tail = F)
  
  time_expert <- as.numeric(na.omit(as.numeric(sim_comb2[i, grep("expert_time_", colnames(sim_comb2))])))
  mu_expert <- as.numeric(na.omit(as.numeric(sim_comb2[i, grep("expert_mean", colnames(sim_comb2))])))
  sd_expert <- as.numeric(rep(sim_comb2[i, grep("sd", colnames(sim_comb2))], length(mu_expert)))
  
  index_selc <- which( expert_time_all%in% time_expert)
  
  mu_expert <-  mu_expert[index_selc]
  sd_expert <-sd_expert[index_selc]
  
  a_beta <- as.numeric(((mu_expert*(1-mu_expert))/pow(sd_expert,2) -1)*mu_expert)
  b_beta <-as.numeric(((mu_expert*(1-mu_expert))/pow(sd_expert,2) -1)*(1-mu_expert))
  
  
  mu_lnorm <- log(mu_expert/sqrt(pow(sd_expert,2)/pow(mu_expert,2) +1 )) 
  sigma_lnorm <- sqrt(log(pow(sd_expert,2)/pow(mu_expert,2) +1 ))
  
  a_gamma <- pow(mu_expert,2)/pow(sd_expert,2)
  b_gamma <- pow(mu_expert,1)/pow(sd_expert,2)
  #plot(seq(0,1,by = 0.01), dbeta(seq(0,1,by = 0.01), a_beta[2], b_beta[2]))
  #plot(seq(0,1,by = 0.01), dlnorm(seq(0,1,by = 0.01), mu_lnorm[2], sigma_lnorm[2]))
  #plot(seq(0,1,by = 0.01), dgamma(seq(0,1,by = 0.01), a_gamma[2], b_gamma[2]))
  #plot(seq(0,1,by = 0.01), dnorm(seq(0,1,by = 0.01), mu_expert[2], sd_expert[2]))
  
 
  #plot(seq(0,1,by = 0.01), dtruncnorm(seq(0,1,by = 0.01),a = 0, b = 1, mu_expert[2], sd_expert[2]))
  
  expert_opinion<- param_expert <- list()
  
  for(j in 1:length(time_expert)){
    
    if(beta_dist){
      param_expert[[j]] <-  data.frame(dist = "beta", wi = 1, param1 = a_beta[j],
                                       param2 = b_beta[j], param3 = NA)
    }else{
      param_expert[[j]] <-  data.frame(dist = "norm", wi = 1, param1 = mu_expert[j],
                                       param2 = sd_expert[j], param3 = NA)
    }
    
    
  }
  expert_opinion$param_expert <- make_data_expert(param_expert, time_expert)
  expert_opinion$times <- time_expert
  expert_opinion$pool <- 1
  
  #folder_name <- paste0("Times- ",paste0( time_expert, collapse = ","), " N-Sim ",samp_size_sim , "/")
  
  plt_name <- paste0(paste0("Cut off - ",cut_off_sim), paste0(" Mean - ", paste0(mu_expert, collapse = ",")),
                     " SD - ",paste0( sd_expert, collapse = ","), paste0(" Shape - ",param_1_sim),
                     paste0(" Scale - ",param_2_sim), " Expert Multiplier - ", expert_mult ,collapse = "")
  #names_survHE <- survHE:::load_availables()
  #model_eval_expertsurv <- mod_all[mod]
  #model_eval_flesurv <- names(which(names_survHE$mle==model_eval_expertsurv))[1]
  
  RSt_MLE_err <- RSt_expert_err <- ISSE_expert <- ISSE_MLE <- RMSE_MLE<- RMSE_expert <- rep(NA,n.boots )
  Surv_mat_expert <- Surv_mat_mle <- matrix(NA, nrow = n.boots, ncol = length(time_seq))
  
  b <- 1
  
  while(b <= n.boots){
  
 # for(b in 1:n.boots){
    
    time_sim <- rweibullPH(samp_size_sim, shape = param_1_sim, scale = param_2_sim)
    data_use <- data.frame(time_event = time_sim) %>% mutate(time = ifelse(time_event > cut_off_sim,
                                                                           cut_off_sim, time_event),
                                                             status = ifelse(time_event > cut_off_sim, 0, 1))
    
    
    #for(q in 1:length(mod_all)){
    if(mod_selc == "rps"){
      fit_mle<-  try({flexsurvspline(formula = Surv(time, status) ~ 1, 
                                     data = data_use,k = n.knots, expert_opinion = NULL)},
                     silent = TRUE)
    }else{
      fit_mle<- try({flexsurvreg(formula = Surv(time, status) ~ 1, 
                                 data = data_use, dist = mod_selc)},
                    silent = TRUE)
    }
    
    #assign(paste0("mle.ests_",mod_selc),fit_expert_vague)
    
    
    if(class(fit_mle)=="try-error"){
      print("Could not fit MLE model")
      next
    }
    
    if(mod_selc == "rps"){
      fit_expert<-  try({flexsurvspline(formula = Surv(time, status) ~ 1, 
                                        data = data_use,k = n.knots, expert_opinion = expert_opinion)},
                        silent = TRUE)
      
    }else{
      fit_expert<- try({flexsurvreg(formula = Surv(time, status) ~ 1, 
                                    data = data_use, dist=mod_selc,expert_opinion = expert_opinion)},
                       silent = TRUE)
    }
    
    #plot(fit_expert, t =0:time_horizon, xlim = c(0,time_horizon))
    #plot(fit_mle, t =0:time_horizon, xlim = c(0,time_horizon), col = "blue")
    
    
    if(class(fit_expert)=="try-error"){
      print("Could not fit Expert model")
      next
    }
    
    #Doesn't matter if the integration is done seperately or at the mean of survival you get same E[R_St]
    # where R_St is restricted mean survival
    #sep_int <- apply(Surv.all.collapsing,2, function(x){sfsmisc::integrate.xy(x = time_seq,fx = x)})
    #mean(sep_int)
    
    
    
    RMSE_MLE_mod <- try({summary(fit_mle, type = "rmst", t = time_horizon)[[1]][1,"est"]}, silent =T)
    RMSE_expert_mod <- try({summary(fit_expert, type = "rmst", t = time_horizon)[[1]][1,"est"]}, silent =T)
    if(class(RMSE_expert_mod )=="try-error"|class(RMSE_MLE_mod)=="try-error"){
      print("Error Integration Failed")
      next
    }
    
    
    Surv_MLE_mod <-  summary(fit_mle, type = "survival", t = time_seq)[[1]][,"est"]
    ISSE_MLE[b] <- sfsmisc::integrate.xy(x = time_seq, fx = (True_St-Surv_MLE_mod)^2)
    RMSE_MLE[b] <- (RMSE_MLE_mod-True_RSt)^2
    RSt_MLE_err[b]<- abs(RMSE_MLE_mod-True_RSt)/True_RSt
    
    Surv_expert_mod <-  summary(fit_expert, type = "survival", t = time_seq)[[1]][,"est"]
    ISSE_expert[b] <- sfsmisc::integrate.xy(x = time_seq, fx = (True_St-Surv_expert_mod)^2)
    RMSE_expert[b] <- (RMSE_expert_mod-True_RSt)^2
    RSt_expert_err[b]<- abs(RMSE_expert_mod-True_RSt)/True_RSt

    
    
    Surv_mat_expert[b,] <- Surv_expert_mod
    Surv_mat_mle[b,] <- Surv_MLE_mod
    
    b <- b +1
  }
  
  folder_name <- paste0("Times- ",paste0( time_expert, collapse = ","), " N-Sim ",samp_size_sim , "/")
  
  png(gen_pathway(paste0(path_output,folder_name,plt_name,".png" )))
  
  plot(x = time_seq, apply(Surv_mat_expert,2, quantile, na.rm = T, probs = 0.5), type = "l", ylim = c(0,1),col = "purple",
       lwd=2.0, ylab = "S(t)", xlab = "Time")
  lines(x = time_seq,apply(Surv_mat_expert,2, quantile, na.rm = T, probs = 0.025), col = "purple", lty= 5)
  lines(x = time_seq,apply(Surv_mat_expert,2, quantile, na.rm = T, probs = 0.975), col = "purple", lty= 5)
  lines(x = time_seq,apply(Surv_mat_mle,2, quantile, na.rm = T, probs = 0.5), col = "red", lty= 1,lwd=2.0)
  lines(x = time_seq,apply(Surv_mat_mle,2, quantile, na.rm = T, probs = 0.025), col = "red", lty= 5)
  lines(x = time_seq,apply(Surv_mat_mle,2, quantile, na.rm = T, probs = 0.975), col = "red", lty= 5)
  lines(x = time_seq,y = True_St,col = "blue",lwd=2.0)
  legend("topright", legend=c("True Survival", "MLE", "Data + Expert"),
         col=c("blue", "red", "purple"), lty=1, cex=0.8)
  for(x in 1:length(time_expert)){
    if(beta_dist){
      quants <- qbeta(c(0.025,0.975),a_beta[x], b_beta[x])
    }else{
      quants <- qnorm(c(0.025,0.975),mu_expert[x], sd_expert[x])
    }
    segments(x0 = time_expert[x], y0 = quants[1], y1 = quants[2], lty= 2)
    
  }
  dev.off()
  
  
  RMSE_res_mat[i,1:3] <- quantile(RMSE_expert, na.rm = T, prob = c(0.025,.5,0.975))
  RMSE_res_mat[i,4:6] <- quantile(RMSE_MLE, na.rm = T, prob = c(0.025,.5,0.975))
  RMSE_res_mat[i,7:8] <- c(mean(RMSE_expert-RMSE_MLE, na.rm  = T), sd(RMSE_expert-RMSE_MLE, na.rm  = T))
  RMSE_res_mat[i,9] <- length(na.omit(RMSE_expert))
  RMSE_res_mat[i,10:12] <- quantile(RMSE_expert/RMSE_MLE, na.rm = T, prob = c(0.025,.5,0.975))
  RMSE_res_mat[i,13] <- sd(RMSE_expert/RMSE_MLE, na.rm = T)
  RMSE_res_mat[i,14:16] <- quantile(RSt_expert_err, na.rm = T, prob = c(0.025,.5,0.975))
  RMSE_res_mat[i,17:19] <- quantile(RSt_MLE_err, na.rm = T, prob = c(0.025,.5,0.975))
  RMSE_res_mat[i,20:22] <- quantile(RSt_expert_err/RSt_MLE_err, na.rm = T, prob = c(0.025,.5,0.975))
  
  
  ISSE_res_mat[i,1:3] <- quantile(ISSE_expert, na.rm = T, prob = c(0.025,.5,0.975))
  ISSE_res_mat[i,4:6] <- quantile(ISSE_MLE, na.rm = T, prob = c(0.025,.5,0.975))
  ISSE_res_mat[i,7:8] <- c(mean(ISSE_expert-ISSE_MLE, na.rm  = T), sd(ISSE_expert-ISSE_MLE, na.rm  = T))
  ISSE_res_mat[i,9] <- length(na.omit(ISSE_expert))
  ISSE_res_mat[i,10:12] <- quantile(ISSE_expert/ISSE_MLE, na.rm = T, prob = c(0.025,.5,0.975))
  ISSE_res_mat[i,13] <- sd(ISSE_expert/ISSE_MLE, na.rm = T)
  
  
}



write.xlsx(ISSE_res_mat,paste0(path_output,"ISSE_res.xlsx"))
write.xlsx(RMSE_res_mat,paste0(path_output,"RMSE_res_mat.xlsx"))


res_mod <- cbind(sim_comb2, Median.Ratio_err = RMSE_res_mat[,21],  Mean.diff_RMSE = RMSE_res_mat[,"Mean Diff"]) %>% mutate(time_points = paste0(expert_time_1, ",",expert_time_2)) %>% 
  group_by(time_points, prop.correct,sd) 
plot(res_mod$Mean.diff_RMSE,res_mod$Median.Ratio_err)
abline(v = 0, h = 1,col="red", lwd=3, lty=2)



res_mod_summary <- res_mod %>%summarize(mean_ratio = mean(Median.Ratio_err),
                             mean_diff = mean(Mean.diff_RMSE))
View(res_mod %>% filter(Median.Ratio_err > 1))


res_mod %>% filter()

png(paste0(path_output,"Density Diff - RMSE.png"))
plot(density(RMSE_expert-RMSE_MLE), main = "Density of difference between RMSE")
dev.off()

png(paste0(path_output,"Density Ratio - Abs RMST.png"))
plot(density(RSt_expert_err/RSt_MLE_err), main = "Density of ratio of the difference in RMST to true value")
dev.off()



write.xlsx(paste0(path_output, "res-final.xlsx"),res_mod_summary )

colnames(RMSE_res_mat)
