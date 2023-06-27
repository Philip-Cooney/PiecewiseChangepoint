
#============================================================================
# Preliminaries - load required packages
# Notes: If the following packages aren't already installed, then
#        they can be  installed from CRAN by typing, for example:
#        install.packages("package-name"). Note that this requires an
#        internet connection.
# For PiecewiseChangepoint package, either the package binary should be
# installed, or run devtools::install_github("Anon19820/PiecewiseChangepoint")
# (if devtools is installed)
#============================================================================

#Pathway to output the results -- Need to Change as Appropriate
pathway <- "~/SimStudy-Res/"
#All these Packages need to be installed for code to run
list.of.packages <-c("PiecewiseChangepoint", "flexsurv","xlsx", "dplyr")
#install.packages(list.of.packages)
library("PiecewiseChangepoint")
library("flexsurv")
library("xlsx")
library("dplyr")

sim.study <- function(n_obs,n_events_req,rate,t_change, max_time =2,
                      n.sims, sims = 10750,burn_in = 750, timescale, lambda, time_seq, right_wrong = F){

  param_true <- list()
  param_true$lambda <- matrix(rate, nrow = 1)
  param_true$changepoint <- matrix(t_change, nrow = 1)

  if(any(is.na(t_change))){
    num.breaks = 0
  }else{
    num.breaks <- length(t_change)

  }
  n.events <- rep(NA, n.sims)
  avg.Haz <- array(NA, dim = c(n.sims, 7))
  avg.change <- array(NA, dim = c(n.sims, 7))
  ISSE_vec <- RMSE_vec <- RSt_err <- matrix(ncol = 8, nrow =  n.sims)
  mod_param <- c("weibull", "lnorm", "gamma", "gompertz", "llogis", "exponential", "gengamma")

  for(i in 1:n.sims){
    print(paste0("Sim ",i," of ",n.sims))
    Sys.sleep(1 / n.sims)

     n_events <- 0


      df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                         num.breaks = num.breaks,rate = rate ,
                         t_change = t_change, max_time = max_time)

      df$status2 <- 1

      surv.obj.df  <- survfit(Surv(time_event ,status2)~1, data = df)
      True_St <- summary(surv.obj.df, t = time_seq)[["surv"]]
      time_seq_km <- summary(surv.obj.df, t = time_seq)[["time"]]
      True_RSt <- integrate.xy(x = time_seq_km, True_St)

      n_events <-  df$status

    n.events[i] <- sum(df$status == 1)

    mod <- collapsing.model(df,n.chains = 1,n.iter = sims, burn_in = burn_in, timescale = timescale )

    most_prob <-  as.numeric(names(mod$prob.changepoint[which.max(mod$prob.changepoint)]))

    index_prob <-  apply(mod$changepoint,1, function(x){length(na.omit(x))==most_prob})
    avg.Haz[i,] <- colMeans(mod[["lambda"]][index_prob,])
    avg.change[i,] <- colMeans(mod[["changepoint"]][index_prob,])
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


sim.study_right_wrong <- function(n_obs,n_events_req,rate,t_change, max_time =2,
                      n.sims, sims = 10750,burn_in = 750, timescale, lambda, time_seq){

  param_true <- list()
  param_true$lambda <- matrix(rate, nrow = 1)
  param_true$changepoint <- matrix(t_change, nrow = 1)

  if(any(is.na(t_change))){
    num.breaks = 0
  }else{
    num.breaks <- length(t_change)

  }
  n.events <- rep(NA, n.sims)
  avg.Haz <- array(NA, dim = c(n.sims, 7))
  avg.change <- array(NA, dim = c(n.sims, 7))
  ISSE_vec <- RMSE_vec <- RSt_err <- matrix(ncol = 10, nrow =  n.sims)
  mod_param <- c("weibull", "lnorm", "gamma", "gompertz", "llogis", "exponential", "gengamma")

  for(i in 1:n.sims){
    print(paste0("Sim ",i," of ",n.sims))
    Sys.sleep(1 / n.sims)

    n_events <- 0


    df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                       num.breaks = num.breaks,rate = rate ,
                       t_change = t_change, max_time = max_time)

    df$status2 <- 1

    surv.obj.df  <- survfit(Surv(time_event ,status2)~1, data = df)
    summary_km <- summary(surv.obj.df, t = time_seq)
    True_St <- summary_km[["surv"]]
    time_seq_km <- summary_km[["time"]]
      True_RSt <- integrate.xy(x = time_seq_km, True_St)

    n_events <-  df$status

    n.events[i] <- sum(df$status == 1)
   if(any(is.na(t_change))){
      num_breaks_correct <- NA
      num_breaks_wrong1 <- 1
      num_breaks_wrong2 <- 2
    }else{
      num_breaks_correct <- length(t_change)
      num_breaks_wrong1 <- length(t_change)-1
      num_breaks_wrong2 <- length(t_change)+1
    }

    if(!is.na(num_breaks_correct)){
      mod_correct <- gibbs.sampler(df, num.breaks = num_breaks_correct, n.iter =sims, burn_in = burn_in, num.chains = 2)
      mod_correct_St <- rowMeans(get_Surv(mod_correct,  time = time_seq_km))
    }else{
      mod_correct_St <- pexp(time_seq_km, sum(df$status)/sum(df$time), lower.tail = FALSE)
    }
    if(!is.na(num_breaks_wrong1)){
      mod_wrong1 <- gibbs.sampler(df, num.breaks = num_breaks_wrong1, n.iter =sims, burn_in = burn_in, num.chains = 2) # Too few
      mod_wrong1_St <- rowMeans(get_Surv(mod_wrong1, time = time_seq_km))

    }else{
      mod_wrong1_St <- pexp(time_seq_km, sum(df$status)/sum(df$time), lower.tail = FALSE)
    }

    if(!is.na(num_breaks_wrong2)){
      mod_wrong2 <- gibbs.sampler(df, num.breaks= num_breaks_wrong2, n.iter =sims, burn_in = burn_in, num.chains = 2) # Too many
      mod_wrong2_St <- rowMeans(get_Surv(mod_wrong2, time = time_seq_km))
    }else{
      mod_wrong2_St <- pexp(time_seq_km, sum(df$status)/sum(df$time), lower.tail = FALSE)
    }


    mod_correct_RSt <- integrate.xy(x = time_seq_km, mod_correct_St)
    mod_wrong1_RSt <- integrate.xy(x = time_seq_km, mod_wrong1_St)
    mod_wrong2_RSt <- integrate.xy(x = time_seq_km, mod_wrong2_St)

    ISSE_vec[i,1] <- integrate.xy(x = time_seq_km, fx = (True_St-mod_correct_St)^2)
    ISSE_vec[i,2] <- integrate.xy(x = time_seq_km, fx = (True_St-mod_wrong1_St)^2)
    ISSE_vec[i,3] <- integrate.xy(x = time_seq_km, fx = (True_St-mod_wrong2_St)^2)

    RMSE_vec[i,1] <- (mod_correct_RSt-True_RSt)^2
    RMSE_vec[i,2] <- (mod_wrong1_RSt-True_RSt)^2
    RMSE_vec[i,3] <- (mod_wrong2_RSt-True_RSt)^2

    RSt_err[i,1]<- abs(mod_correct_RSt-True_RSt)/True_RSt
    RSt_err[i,2]<- abs(mod_wrong1_RSt-True_RSt)/True_RSt
    RSt_err[i,3]<- abs(mod_wrong2_RSt-True_RSt)/True_RSt

    if(RSt_err[i,1] > 1 ){
      plot(time_seq_km, y = True_St)
      lines(time_seq_km,y = mod_correct_St, col = "red" )
    }

    for(b in 1:length(mod_param)){
      fit.parm <- flexsurvreg(formula = Surv(time, status) ~ 1, data = df, dist=mod_param[b])
      mod_St <- summary(fit.parm, t = time_seq_km)[[1]]$est
      mod_RSt <- integrate.xy(x = time_seq_km, fx = mod_St)
      ISSE_vec[i,3+b] <- integrate.xy(x = time_seq_km, fx = (True_St-mod_St)^2)
      RMSE_vec[i,3+b] <- (mod_RSt-True_RSt)^2
      RSt_err[i,3+b]<- abs(mod_RSt-True_RSt)/True_RSt
    }

    if(RSt_err[i,1] > 1 ){
      fit.parm <- flexsurvreg(formula = Surv(time, status) ~ 1, data = df, dist=mod_param[1])
      mod_St <- summary(fit.parm, t = time_seq_km)[[1]]$est

      plot(time_seq_km, y = True_St)
      lines(time_seq_km,y = mod_correct_St, col = "red" )
      lines(time_seq_km,y = mod_St, col = "blue" )
    }


  }
  if(any(RSt_err[,1] < 0)){
    RSt_err <- RSt_err[-which(RSt_err[,1] < 0),]
  }

 return(list( ISSE_vec = ISSE_vec,
              RMSE_vec = RMSE_vec,
              RSt_err = RSt_err,
              change_num = c(num_breaks_correct,num_breaks_wrong1,num_breaks_wrong2)))
}



#No Change-point

rate1<- c(0.25)
rate2<- c(0.5)
rate3<- c(0.75)
sample_vec <- c(100,200)
censor_vec <- c(0,0.5)
rate_list <- list()
for(i in 1:3){
  current_rate <-get(paste0("rate",i))
  rate_list[[i]] <- current_rate
}

sims_eval <- 500

#For results in Table 1 & Table 2 in Supplementary Appendix use function sim.study
#For results in Table 3 in Supplementary Appendix use function sim.study_right_wrong
for(q in 1:length(censor_vec) ){
  for(j in 1:length(sample_vec)){
    for(i in 1:length(rate_list)){

      current_name <- paste0("Collapsing_rate_",paste0(rate_list[[i]],collapse = "_"),"_size_",sample_vec[j],"_censor_",censor_vec[q])
      assign(current_name, sim.study_right_wrong(n_obs = sample_vec[j],
                                     n_events_req= round(sample_vec[j]*(1-censor_vec[q])),
                                     rate = rate_list[[i]],
                                     max_time = 2,
                                     t_change = NA,
                                     n.sims = sims_eval,
                                     sims =2500,
                                     burn_in = 750,
                                     timescale = "months",
                                     lambda = 1,
                                     time_seq = c(0:15)))
       saveRDS(get(current_name), file = paste0(pathway,current_name,".rds"))

    }
  }
}


#One and Two Change-point models


rate1 <- c(0.5,0.75)
rate2 <- c(0.25,0.75)
rate3 <- c(0.75,0.5)
rate4 <- c(0.75,0.25)

rate5 <- c(0.25,0.5,0.75)
rate6 <- c(0.75,0.5,0.25)
rate7 <- c(0.75,0.2,0.75)
rate8 <- c(0.2,0.75,0.2)

time1 <- c(0.5,1)

censor_vec <- c(0) #c(0,0.33)
sample_vec <- 300#c(300,500,1000)
rate_list <- list()

for(i in 1:8){
  current_rate <-get(paste0("rate",i))
  rate_list[[i]] <- current_rate
}


for(q in 1:length(censor_vec) ){
  for(j in 1:length(sample_vec)){
    for(i in 1:length(rate_list)){

      current_name <- paste0("Collapsing_rate_",paste0(rate_list[[i]],collapse = "_"),"_size_",sample_vec[j],"_censor_",censor_vec[q])
      assign(current_name, sim.study_right_wrong(n_obs = sample_vec[j],
                                     n_events_req= round(sample_vec[j]*(1-censor_vec[q])),
                                     rate = rate_list[[i]],
                                     max_time = 2,
                                     t_change = c(0.5,1)[1:(length(rate_list[[i]])-1)],
                                     n.sims = sims_eval,
                                     sims =2500,
                                     burn_in = 750,
                                     timescale = "months",
                                     lambda = 1,
                                     time_seq = c(0:15)))

      saveRDS(get(current_name), file = paste0(pathway,current_name,".rds"))

    }
  }
}



for(q in 1:length(censor_vec) ){
  for(j in 1:length(sample_vec)){
    for(i in 1:length(rate_list)){

      current_name <- paste0("Collapsing_rate_",paste0(rate_list[[i]],collapse = "_"),"_size_",sample_vec[j],"_censor_",censor_vec[q])
      print(colMeans(get(current_name)[["RSt_err"]][,1:3]))

    }
  }
}

## Plot results


#Comparision versus Chapel - See Chapple Simulation



