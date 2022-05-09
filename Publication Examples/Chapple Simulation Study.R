#devtools::load_all()


library("survival")
library("BayesReversePLLH")
library("PiecewiseChangepoint")
library("sfsmisc")

max_time =  2
rate1 <- c(0.25)
rate2 <- c(0.5)
rate3 <- c(0.75)

rate4 <- c(0.5,0.75)
rate5 <- c(0.25,0.75)
rate6 <- c(0.75,0.5)
rate7 <- c(0.75,0.25)

rate8 <- c(0.25,0.5,0.75)
rate9 <- c(0.75,0.5,0.25)
rate10 <- c(0.75,0.2,0.75)
rate11 <- c(0.2,0.75,0.2)


censor_vec <- c(0,0.33)
censor_vec <- c(0)
sample_vec <- c(100,300,500)
#sample_vec <- c(300,500,1000)
rate_list <- list()

for(i in 1:11){
  current_rate <-get(paste0("rate",i))
  rate_list[[i]] <- current_rate
}


pathway <- "~/Simulation Study 2022/Chapple Comparison/"


sim.study_comp <- function(n_obs,n_events_req,rate,t_change, max_time =2,
                           n.sims, sims = 10750,burn_in = 750, lambda.prior = 1,
                           n.chains = 1, time_horizon =25){


  time_seq <- seq(0, time_horizon, length.out = 1500)
  if(any(is.na(t_change))){
    num.breaks = 0
  }else{
    num.breaks <- length(t_change)
  }

  # 1 Generate holding variables ----

  selc_mod_Collasping <- selc_mod_Chapple <- n.events <-  rep(NA, n.sims)
  RMSE_Collapsing_time_horizon <- RMSE_Collapsing_restrict <- RMSE_Chapple_time_horizon <- RMSE_Chapple_restrict <- rep(NA, n.sims)
  ISSE_Collapsing_time_horizon <- ISSE_Collapsing_restrict <- ISSE_Chapple_time_horizon <- ISSE_Chapple_restrict <- rep(NA, n.sims)

  # 2 True Values ----

  ## 2.1 True Survival ----

  if(any(is.na(t_change))){
    True_Surv <- pexp(time_seq, rate, lower.tail = F)
  }else{
    True_Surv <- GetSurvPEH(x =time_seq,s = c(0,t_change,100), lam =  log(rate), J  = length(t_change))

  }
  True_Surv[1] <- 1


  ## 2.2 RSt (Restricted Mean Survival) ----

  True_RSt_restrict <- sfsmisc::integrate.xy(fx = True_Surv[time_seq<= max_time],
                                             x = time_seq[time_seq<= max_time])

  True_RSt_time_horizon <- sfsmisc::integrate.xy(fx = True_Surv,
                                                 x = time_seq)


  # 3 Run Simulations ----


  for(i in 1:n.sims){
    print(paste0("Sim ",i," of ",n.sims))


    df <- PiecewiseChangepoint::gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                                             num.breaks = length(t_change),rate = rate ,
                                             t_change = t_change, max_time = max_time)

    n.events[i] <- sum(df$status)


    Collapsing_Model <-    collapsing.model(df,
                                            n.iter = sims,
                                            burn_in = burn_in,
                                            n.chains = n.chains,
                                            alpha.hyper = 1,
                                            beta.hyper1 = 1,
                                            beta.hyper2 = 1,
                                            lambda.prior = lambda.prior)



    Chapple_Model <- BayesPiecewiseHazard(Y = df$time, I1 =df$status, Poi = lambda.prior,  B = sims*n.chains)
    for( j in 1:nrow(Chapple_Model[[1]])){ # Requires a large number to calculate survival for time_horizon
      Chapple_Model[[1]][j,which.max(Chapple_Model[[1]][j, ])] <- 100
    }

    Surv.all.Chapple  <- GetALLSurvPEH(time_seq, Chapple_Model)
    Surv.all.Chapple[,1] <- 1

    Surv.all.collapsing  <- PiecewiseChangepoint:::get_Surv(object = Collapsing_Model,time = time_seq)
    #Chapple's function is slower but gives the same results
    # Surv.all.collapsing  <- GetSurvPEH_collapse(time_seq, mod = mod)
    # Surv.all.collapsing[,1] <- 1



    ## 3.1 Calculate ISSE ----

    #Can use the rectangular integration from Chapple
    #Both approaches are the same for ISSE, however,
    #sfsmisc::integrate.xy is much more accurate for estimating restricted mean survival:

    #time_seq <- seq(0, max(Times), length.out = 1000)
    #dtheta <- mean(diff(time_seq))
    #sum(((True_Surv-colMeans(Surv.all.Chapple))^2)*dtheta)

    ISSE_Chapple_time_horizon[i] <- sfsmisc::integrate.xy(x = time_seq, fx = (True_Surv-colMeans(Surv.all.Chapple))^2)

    ISSE_Chapple_restrict[i] <- sfsmisc::integrate.xy(x = time_seq[time_seq<= max_time],
                                                      fx = ((True_Surv-colMeans(Surv.all.Chapple))^2)[time_seq<= max_time])


    ISSE_Collapsing_time_horizon[i] <- sfsmisc::integrate.xy(x = time_seq, fx = (True_Surv-rowMeans(Surv.all.collapsing))^2)
    ISSE_Collapsing_restrict[i] <- sfsmisc::integrate.xy(x = time_seq[time_seq<= max_time],
                                                         fx = ((True_Surv-rowMeans(Surv.all.collapsing))^2)[time_seq<= max_time])

    ## 3.2 Calculate RSME ----

    #Doesn't matter if the integration is done seperately or at the mean of survival you get same E[R_St]
    # where R_St is restricted mean survival
    #sep_int <- apply(Surv.all.collapsing,2, function(x){sfsmisc::integrate.xy(x = time_seq,fx = x)})
    #mean(sep_int)
    RMSE_Collapsing_time_horizon[i] <- (sfsmisc::integrate.xy(x = time_seq,
                                                              fx = rowMeans(Surv.all.collapsing))-True_RSt_time_horizon)^2

    RMSE_Collapsing_restrict[i] <- (sfsmisc::integrate.xy(x = time_seq[time_seq<= max_time],
                                                          fx = rowMeans(Surv.all.collapsing)[time_seq<= max_time]) - True_RSt_restrict)^2


    RMSE_Chapple_time_horizon[i] <- (sfsmisc::integrate.xy(x = time_seq,
                                                           fx = colMeans(Surv.all.Chapple))-True_RSt_time_horizon)^2

    RMSE_Chapple_restrict[i] <- (sfsmisc::integrate.xy(x = time_seq[time_seq<= max_time],
                                                       fx = colMeans(Surv.all.Chapple)[time_seq<= max_time])
                                 - True_RSt_restrict)^2


    # 4 Most Probable model----

    selc_mod_Collasping[i] <- as.numeric(names(which.max(Collapsing_Model$prob.changepoint)))
    selc_mod_Chapple[i] <- as.numeric(names(table(Chapple_Model[[3]])[which.max(table(Chapple_Model[[3]]))]))

    # 5 Survival Plot ----

    # Not used but provided for reference

    #https://stackoverflow.com/questions/43173044/how-to-compute-the-mean-survival-time

    # km <- survfit(Surv(time,status)~1, data = df)
    # sum(diff(c(0,km$time))*c(1,km$surv[1:(length(km$surv)-1)]))
    # survival:::survmean(true_km, rmean=max(Times_event))[[1]]["*rmean"]
    #
    #
    # plot(km)
    # lines(y = Survs, x = time_seq)
    # lines(y = colMeans(Surv.all.Chappel), x = time_seq, col = "red")
    # lines(y = colMeans(Surv.all.collapsing), x = time_seq, col = "blue")


  }

  #Highest posterior model
  prob.correct_Collapsing <- mean(selc_mod_Collasping == num.breaks)
  prob.correct_Chapple <- mean(selc_mod_Chapple == num.breaks)

  list_result <- list()

  list_result[["Collapsing"]] <- list(selc_mod = selc_mod_Collasping,
                                      prob.correct = prob.correct_Collapsing,
                                      RMSE_time_horizon = RMSE_Collapsing_time_horizon,
                                      RMSE_restrict = RMSE_Collapsing_restrict,
                                      ISSE_time_horizon = ISSE_Collapsing_time_horizon,
                                      ISSE_restrict = ISSE_Collapsing_time_horizon)


  list_result[["RJMCMC"]] <- list(selc_mod = selc_mod_Chapple,
                                  prob.correct = prob.correct_Chapple,
                                  RMSE_time_horizon = RMSE_Chapple_time_horizon,
                                  RMSE_restrict = RMSE_Chapple_restrict,
                                  ISSE_time_horizon = ISSE_Chapple_time_horizon,
                                  ISSE_restrict = ISSE_Chapple_time_horizon)

  list_result[["n_events"]] <- n.events

  return(list_result)

}

#ISsue Row index is out of bounds: [index=1; row extent=1].

#Occurs when lambda_df row is equal to 1
for(q in 1:length(censor_vec) ){
  for(j in 1:length(sample_vec)){
    for(i in 1:length(rate_list)){

      if(length(rate_list[[i]]) ==1){
        t_change_curr <- NA
      }else{
        t_change_curr  <- c(0.5,1)[1:(length(rate_list[[i]])-1)]
      }

      current_name <- paste0("Compare_RJMCMC_",paste0(rate_list[[i]],collapse = "_"),"_size_",sample_vec[j],"_censor_",censor_vec[q])
      print(current_name)

      assign(current_name, sim.study_comp(n_obs = sample_vec[j],
                                          n_events_req = round(sample_vec[j]*(1-censor_vec[q])),
                                          rate = rate_list[[i]],
                                          t_change = t_change_curr,
                                          max_time = max_time,
                                          n.sims = 100,
                                          sims = 1500,
                                          burn_in = 5,
                                          lambda.prior = 1))

      saveRDS(get(current_name), file = paste0(pathway,current_name,".rds"))

    }
  }
}

for(i in 1:length(rate_list)){
  dir.create(paste0(pathway, paste0("Compare_RJMCMC_",paste0(rate_list[[i]],collapse = "_")), "_plots"))
}



names_sq_error <- c("RMSE_time_horizon","RMSE_restrict","ISSE_time_horizon","ISSE_restrict")
df_correct <-  NULL
for(i in 1:length(rate_list)){

  folder_output <- paste0(pathway, paste0("Compare_RJMCMC_",paste0(rate_list[[i]],collapse = "_")), "_plots")
  df_mod <- df_output <- NULL

  for(q in 1:length(censor_vec) ){
    for(j in 1:length(sample_vec)){


      if(length(rate_list[[i]]) ==1){
        t_change_curr <- NA
      }else{
        t_change_curr  <- c(0.5,1)[1:(length(rate_list[[i]])-1)]
      }

      current_name <- paste0("Compare_RJMCMC_",paste0(rate_list[[i]],collapse = "_"),"_size_",sample_vec[j],"_censor_",censor_vec[q])

      scen_eval <- paste0("Censoring ",censor_vec[q], "; Size ", sample_vec[j]  )
      scen_eval_list <- readRDS(file = paste0(pathway,current_name,".rds"))

      df_correct_temp <-  data.frame(Scenario = paste0(scen_eval,paste0("; Hazard c(", paste0(rate_list[[i]], collapse = ","),")")))
      df_correct_temp$Collaping <- scen_eval_list[[1]]$prob.correct
      df_correct_temp$RJMCMC <- scen_eval_list[[2]]$prob.correct

      df_correct <- rbind(df_correct, df_correct_temp)
      for(x in 1:2){#Number of models
        df_temp <- data.frame(Scen = scen_eval, Model = names(scen_eval_list)[x])
        for(p in 1:4){#Number of measures
          df_temp <- cbind(df_temp,scen_eval_list[[x]][[names_sq_error[p]]])
        }
        names(df_temp) <- c("Scenario", "Model", names_sq_error)
        df_mod <- rbind(df_mod, df_temp)
      }
      df_output  <- rbind(df_output, df_mod)
    }
  }

  for( r in 1:4){
    index <- grep(names_sq_error[r], colnames(df_output) )
    plot_save <- ggplot(df_output, aes(x=Scenario, y=df_output[,index], fill=Model)) +
      geom_boxplot()+
      ylim(c(0,.1))+
      ylab(gsub("_", " ",names_sq_error[r]))+
      theme(axis.text.x = element_text(angle = 45, hjust=1))


    ggsave(paste0(folder_output,"/", names_sq_error[r],".png" ), height = 6, width = 8)

  }


}




# BACKUP----

## 1 Validate Survival Functions ----
# Ensuring both functions give the same survival


### 1.1 Multiple Change-points ----
ncp <- 2
cp_sim <-runif(ncp)

lambda_df <- matrix(runif(ncp+1), nrow = 1)
changepoint_df <- matrix(cp_sim[order(cp_sim)], nrow = 1)

object <- list(lambda = lambda_df,
               changepoint = changepoint_df)

abs(as.numeric(get_Surv(object, time = seq(0,10, by = 0.5))) -
      GetSurvPEH(x =seq(0,10, by= 0.5),s = c(0,changepoint_df,100), lam =  log(lambda_df), J  = length(cp_sim))) < 0.0001

### 1.2 No Change-points ----

lambda_df <- matrix(runif(ncp+1), nrow = 1)[,1, drop = F]
changepoint_df <- matrix(NA, nrow = 1)

object <- list(lambda = lambda_df,
               changepoint = changepoint_df)

as.numeric(get_Surv(object, time = seq(0,10, by = 0.5)))
pexp(seq(0,10, by = 0.5), rate = lambda_df, lower.tail = F)

## 2 Calculate Time Horizon ----
# Calculate appropriate get time-horizon for the simulations
# Ensuring that less than 1% of population is surviving

source(paste0("~/Simulation Study 2022/Chapple Comparison/Chapple Simulation Backup Functions.R"))
time_horizon <- 25
for(i in 1:length(rate_list)){
  if(length(rate_list[[i]]) ==1){
    t_change_curr <- NA
    surv_time <- pexp(time_horizon, rate =rate_list[[i]], lower.tail = F)
  }else{
    t_change_curr  <- c(0.5,1)[1:(length(rate_list[[i]])-1)]
    surv_time <-  GetSurvPEH(x =time_horizon,s = c(0,t_change_curr,100), lam =  log(rate_list[[i]]), J  = length(t_change_curr))
  }
  print(surv_time<0.01)

}


## 3 Simulating Events ----

#Approach used by Chapple can give a situation where many of the censorsed observations
# have time 0 otherwise they would have negative survival
# By selecting the correct pooling (so that observations accure togehter )
# and a high accrual rate this can be avoided, however, I prefer my original approach,
# as you fix the "length" of the study. The disadvantage with my approach is that you
# have censoring throughout the study and at end of follow up and therefore the percentage
# censored refers to the proportion who are censored throughout the observed study time.

PCT <- .50
n <- 200

Rate_piece = c(0.75, 0.2,0.75)
Rate_ACC <-mean(Rate_piece)*20/PCT
chng = c(0.5,1)
B = 1000
factor_pool = n/(20*PCT)
censor_time <- GetCensorTimePIECE(PCT = PCT,
                                  n = n,
                                  Rate_ACC =Rate_ACC,
                                  Rate_piece = Rate_piece,
                                  chng = chng,
                                  B = B,
                                  factor_pool)

nboots <- 5000
events <- rep(NA, nboots)
for(i in 1:nboots){

  ACC <- cumsum(rexp(n/factor_pool, Rate_ACC))
  ACC <- sample(ACC, n, replace = T)

  TIMES<- PiecewiseChangepoint:::rpwexp(n, lam  =Rate_piece, s = chng)
  df_times <- data.frame(TIMES, ACC)

  df_final <- df_times %>% mutate(TIME_EVENT = ifelse(TIMES+ACC > censor_time,
                                                      censor_time - ACC, TIMES),
                                  status = ifelse(TIMES+ACC > censor_time,
                                                  0,1))

  if(any(df_final$TIME_EVENT <0)){
    stop("Negative Survival times, need a larger accural rate")
  }
  events[i] <- sum(df_final$status)
}

mean(events)/n
hist(events)
plot(survfit(Surv(TIME_EVENT, status)~1, data = df_final))


## 4 Integration of the Survival function ----
# I'm not sure how Chapple calculates the Area under the curve, however,
# if we use the rectangular box method it is less accurate than the function
# Shown by testing versus the analytic value of the Restricted Mean Survival


rate <- c(0.75,0.2,0.75)
t_change <- c(0.5,1)
#sfsmisc::integrate.xy

#integral cal exp(-lambdax)

time_seq <- seq(0, time_horizon, length.out = 1500)
dtheta <- mean(diff(time_seq))

Survs <- GetSurvPEH(x =time_seq,s = c(0,t_change,100), lam =  log(rate), J  = length(t_change))
Survs[1] <- 1
#TRUE
#exp(-lambdax)
t_change_diff <- diff(c(0, t_change, time_horizon))

cum_haz <- 0
cum_haz_seq <- c(0,t_change_diff*rate)
St_thres <- exp(-cumsum(cum_haz_seq))
#don't need the last St_thres
int_surv<- rep(NA, length(rate))
for(i in 1:length(rate)){
  int_surv[i] <- (exp(-0*rate[i])/rate[i] - exp(-t_change_diff[i]*rate[i])/rate[i])*St_thres[i]
}

sum(int_surv)
sfsmisc::integrate.xy(x =time_seq, fx = Survs)
sum(Survs*dtheta)



