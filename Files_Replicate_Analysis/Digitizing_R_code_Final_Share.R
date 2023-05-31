#============================================================================
# Preliminaries - load required packages
# Notes: If the following packages aren't already installed, then 
#        they can be  installed from CRAN by typing, for example: 
#        install.packages("package-name"). Note that this requires an 
#        internet connection. Installing rstan for the first time requires
#        a few additional steps, see http://mc-stan.org/interfaces/rstan
#        for details. Similarly installing rjags requires additional steps, 
#        see https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781119287995.app1  
#
# For PiecewiseChangepoint package, either the package binary should be 
# installed, or run devtools::install_github("Anon19820/PiecewiseChangepoint")
# (if devtools is installed)
#============================================================================


list.of.packages <- c("PiecewiseChangepoint", "flexsurv","xlsx", "dplyr", "ggplot2", "sjstats", "rstan", "R2jags", "rjags")
install.packages(list.of.packages)
library("PiecewiseChangepoint")
library("flexsurv")
library("xlsx")
library("dplyr")
library("ggplot2")
library("sjstats")
library("rstan")
library("R2jags")

#Helper function
compare_dco <- function (all_surv_mods, old_dco, new_dco, km_risk = 0.1, col_new_dco = "grey"){
  result.km_old <- survfit(Surv(time, status) ~ 1, data = old_dco)
  if (!is.null(km_risk)) {
    max_time <- result.km_old$time[max(which(result.km_old$n.risk/result.km_old$n >= 
                                               km_risk))]
  }
  result.km <- survfit(Surv(time, status) ~ 1, data = new_dco)
  km.data <- data.frame(cbind(result.km[[c("time")]], result.km[[c("surv")]], 
                              result.km[[c("upper")]], result.km[[c("lower")]]))
  colnames(km.data) <- c("time", "survival", "upper", "lower")
  if (is.null(km_risk)) {
    km.data <- km.data %>% dplyr::filter(time >= max(old_dco$time))
  }
  else {
    km.data <- km.data %>% dplyr::filter(time >= max_time)
  }
  all_surv_mods$plot_Surv_all + geom_step(data = km.data, aes(x = time, 
                                                              y = survival),
                                          colour = col_new_dco, inherit.aes = F) + 
    geom_step(data = km.data,aes(x = time, y = upper), colour = col_new_dco, linetype = "dashed",inherit.aes = F) +
    geom_step(data = km.data, aes(x = time,y = lower), colour = col_new_dco, linetype = "dashed", inherit.aes = F)
}

 

# Create a pathway which will import files and export all the results
#All the required files exist in the "Files_Replicate_Analysis" folder (if downloaded from Github).
path<-"~/Files_Replicate_Analysis/"

# 1 Read Data and generate pseudo-IPD  ----

#Read in the Conditional Death Probability for the General Population
#https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/mortalityratesqxbysingleyearofage
Conditional_Death_df <- read.xlsx(paste0(path, "Conditional_Death_UK.xlsx"),1)
time_horizon <- 100 # Max age at which survival probabilities is considered

folder_name <- excel_file <-c(
"TA396_D+T_OS_Initial_V",
"TA396_D+T_OS_Initial_D",
"TA396_D+T_PFS_Initial_V",
"TA396_D+T_PFS_Initial_D",
"TA396_V_PFS_Initial",
"TA396_D_PFS_Initial",
"TA428_PEM_OS_Initial",
"TA269_V_PFS_Initial",
"TA269_D_PFS_Initial",
"TA269_V_OS_Initial",
"TA269_D_OS_Initial",
"TA428_PEM_OS_Update",
"TA268_GP100_OS_Initial",
"TA268_GP100_PFS_Initial",
"TA268_Ipi_PFS_Initial",
"TA268_Ipi_OS_Initial",
"TA347_Doce_OS_Initial",
"TA347_N_OS_Initial",
"TA347_Doce_PFS_Initial",
"TA347_N_PFS_Initial",
"TA396_D+T_OS_Update",
"TA396_D+T_PFS_Update",
"TA447_PEM_OS_Initial",
"TA447_PEM_OS_Update")

for(i in 1:length(excel_file)){
  survival_df <- read.xlsx(paste0(path, folder_name[i],"/",excel_file[i],".xlsx"),
                           "survival")
  nrisk_df <- read.xlsx(paste0(path, folder_name[i],"/",excel_file[i], ".xlsx"),
                        "nrisk")
  
  write.table(survival_df, paste0(path, folder_name[i],"/survival.txt"), row.names = F)
  
  write.table(nrisk_df, paste0(path, folder_name[i],"/nrisk.txt"), row.names = F)
  
  digitise(paste0(path,folder_name[i],"/survival.txt"),
           paste0(path,folder_name[i],"/nrisk.txt"),
           km_output = paste0(path,folder_name[i], "/KMdata.txt"),
           ipd_output = paste0(path,folder_name[i], "/IPDdata.txt"))
  
  digitized_IPD <- read.table(paste0(path,folder_name[i],"/IPDdata.txt"), header = T)
  digitized_IPD <- data.frame(apply(digitized_IPD,2, as.numeric))
  digitized_IPD$arm <- folder_name[i]
  assign(folder_name[i], digitized_IPD)
  
  km.fit <- survfit(Surv(time, status)~1, data.frame(digitized_IPD))
  
  # plot(km.fit, xaxt = "n", yaxt = "n")
  # axis(1, at=0:45, labels=0:45)
  # axis(2, at=seq(0,1, by = 0.1), labels=seq(0,1, by = 0.1))
  #
  median_km <-
    round(survival:::survmean(km.fit, scale = 1, 1)$matrix[["median"]],digits = 2)
  
  plt_time <- (round(max(digitized_IPD[,1]),0)+1)
  
  #Plot the Kaplan Meier curve
  png(paste0(path, folder_name[i],"/KM curve.png"), width = 960, height = 480)
  plot(km.fit, main = paste0("KM curve : ", excel_file[i]),
       sub = paste0("Median : ",median_km),
       xlab = "Months",
       ylab = "S(t)",
       xaxt = "n",
       yaxt = "n")
  axis(1, at=0:plt_time, labels=0:plt_time)
  axis(2, at=seq(0,1, by = 0.1), labels=seq(0,1, by = 0.1))
  points(survival_df$time, y = survival_df$survival, col = "red")
  grid(nx = NULL, ny = NULL,
       lty = 2,      # Grid line type
       col = "gray", # Grid line color
       lwd = 2)
  rm(survival_df, nrisk_df)
  dev.off()
}

n.iter <-20750
burn_in <- 750
n.chains <- 2

n.iter.jags <- 5000
n.burnin.jags <- n.iter.jags/10
seed.val <- 123

# 2 TA268 ----

## 2.1 OS ----

TA268_GP100_OS_model <- collapsing.model(TA268_GP100_OS_Initial,
                                         n.iter = n.iter,
                                         burn_in = burn_in,
                                         n.chains = n.chains,
                                         timescale = "years",
                                         seed.val = seed.val)

TA268_Ipi_OS_model <- collapsing.model(TA268_Ipi_OS_Initial,
                                       n.iter = n.iter,
                                       burn_in = burn_in,
                                       n.chains = n.chains,
                                       timescale = "years",
                                       seed.val = seed.val)

## 2.2 PFS ----


TA268_Ipi_PFS_model <- collapsing.model(TA268_Ipi_PFS_Initial,
                                        n.iter = n.iter,
                                        burn_in = burn_in,
                                        n.chains = n.chains,
                                        timescale = "years",
                                        max.num.breaks = 10,
                                        seed.val = seed.val)

TA268_GP100_PFS_model <- collapsing.model(TA268_GP100_PFS_Initial,
                                          n.iter = n.iter,
                                          burn_in = burn_in,
                                          n.chains = n.chains,
                                          timescale = "years",
                                          seed.val = seed.val)




# 2 TA269 ----

## 2.1 PFS ----

TA269_V_PFS_model <- collapsing.model(TA269_V_PFS_Initial,
                                      n.iter = n.iter,
                                      burn_in = burn_in,
                                      n.chains = n.chains,
                                      timescale = "months",
                                      seed.val = seed.val)

TA269_D_PFS_model <- collapsing.model(TA269_D_PFS_Initial,
                                      n.iter = n.iter,
                                      burn_in = burn_in,
                                      n.chains = n.chains,
                                      timescale = "months",
                                      seed.val = seed.val)


# 3 TA347 ----

## 3.1 OS ----
TA347_Doce_OS_model <- collapsing.model(TA347_Doce_OS_Initial,
                                        n.iter = n.iter,
                                        burn_in = burn_in,
                                        n.chains = n.chains,
                                        timescale = "months",
                                        seed.val = seed.val)


TA347_N_OS_model <- collapsing.model(TA347_N_OS_Initial,
                                     n.iter = n.iter,
                                     burn_in = burn_in,
                                     n.chains = n.chains,
                                     timescale = "months",
                                     seed.val = seed.val)
#plot(TA347_N_OS_model, max_predict = 50, chng.num = 1)
#plot(TA347_N_OS_model, max_predict = 50, chng.num = 0)

## 3.2 PFS ----

TA347_Doce_PFS_model <- collapsing.model(TA347_Doce_PFS_Initial,
                                         n.iter = n.iter,
                                         burn_in = burn_in,
                                         n.chains = n.chains,
                                         timescale = "days",
                                         seed.val = seed.val)


TA347_N_PFS_model <- collapsing.model(TA347_N_PFS_Initial,
                                      n.iter = n.iter,
                                      burn_in = burn_in,
                                      n.chains = n.chains,
                                      timescale = "days",
                                      seed.val = seed.val)

# 4 TA396 ----

## 4.1 OS ----

### Dabrafenib + Trametinib (DT) ----

#We combine Dabrafenib + Trametinib from both COMBI-D and COMBI-V
TA396_DT_OS_Initial_df <- rbind(`TA396_D+T_OS_Initial_D`,`TA396_D+T_OS_Initial_V`)

#Fit the changepoint model
TA396_DT_OS_Initial_model <- collapsing.model(TA396_DT_OS_Initial_df,
                                      n.iter = n.iter, # Change to 20750
                                      burn_in = burn_in,
                                      n.chains = n.chains,
                                      timescale = "months",
                                      seed.val = seed.val)

# Plot the posterior distribution of the changepoint
TA396_DT_OS_Initial_posterior_chng <- plot.pos_changepoint(TA396_DT_OS_Initial_model,
                                                                      breaks = seq(from = 0, to = 25, by = 5))+theme_light()

max_prob_chng <-   as.numeric(names(which.max(TA396_DT_OS_Initial_model$prob.changepoint)))

# Compare to other survival models and adjust for General Population Mortaility
age_baseline_TA396 <- round(54.475, digits = 0)
prop_male_TA396 <- 0.567

Conditional_Death_df_temp <- Conditional_Death_df
Conditional_Death_df_temp[, "mix_prob"] <- Conditional_Death_df_temp[,2]*prop_male_TA396 + Conditional_Death_df_temp[,3]*(1-prop_male_TA396)
Conditional_Death_df_temp <- Conditional_Death_df_temp %>% filter(age >= age_baseline_TA396 & age <= time_horizon)





#Data (or Model Cycle length) is in months so need to adjust Conditional Death Probabilities
#If Survival data/Cycle Length is in: 
# - Years time_factor is = 1
# - Months time_factor is 12
# - Weeks time_factor is round(365.25/7, digits = 0)

#We need to convert probability conditional death to a rate and
# divide to the appropriate cycle length 
#See Fleurence et al "Rates and Probabilities in Economic Modelling."https://doi.org/10.2165/00019053-200725010-00002
#Equation 2; t = timefactor
time_factor_TA396 <- 12

Conditional_Death_df_temp$mix_haz <- -log(1-Conditional_Death_df_temp$mix_prob)/time_factor_TA396

gmp_haz_vec_TA396 = rep(Conditional_Death_df_temp$mix_haz,each = time_factor_TA396)
gmp_haz_df_TA396 <- data.frame(time = 1:length(gmp_haz_vec_TA396), hazard = gmp_haz_vec_TA396)

TA396_DT_OS_Initial_all_mods <- compare.surv.mods(TA396_DT_OS_Initial_model, max_predict = 7*12,
                                                   n.iter.jags = n.iter.jags, #Should be 5000
                                                   n.thin.jags = 1,
                                                   n.burnin.jags = n.burnin.jags,
                                                   chng.num = max_prob_chng,
                                                   gmp_haz_df =gmp_haz_df_TA396, 
                                                  col_km = "grey" )

TA396_DT_OS_gg_final <-compare_dco(TA396_DT_OS_Initial_all_mods, TA396_DT_OS_Initial_df,`TA396_D+T_OS_Update`)+ xlab("Time (Months)")
ggsave(filename =paste0(path,"pub_plots_tabs/TA396_DT_OS_Surv.png"),plot = TA396_DT_OS_gg_final, height = 5.59, width = 7)

#Compare the RMST time using the bootstrap of the KM curve
TA396_DT_OS_boot <- PiecewiseChangepoint:::compare_boot_sims(mod_parametric_orig = TA396_DT_OS_Initial_all_mods,
                                       follow_up_data =`TA396_D+T_OS_Update`)

write.xlsx(TA396_DT_OS_boot, paste0(path,"pub_plots_tabs/TA396_DT_OS_AUC.xlsx"))


# Revewer 2 comment: Add Many knots/change-points

km_update <- survfit(Surv(time, status)~1, data = `TA396_D+T_OS_Update`)
spline_many <- flexsurvspline(Surv(time, status)~1, data =TA396_DT_OS_Initial_df, k = 5)

gibbs_sampler<- PiecewiseChangepoint:::gibbs.sampler(df =TA396_DT_OS_Initial_df, num.breaks = 10, n.iter =5000,
                                     burn_in = 100, num.chains = 2, 
                                     n.thin = 1, alpha.hyper = 1, beta.hyper = 1)

colMeans(gibbs_sampler$changepoint)
max(TA396_DT_OS_Initial_df$time)
exp(spline_many$knots)
Surv_gibbs <- PiecewiseChangepoint:::get_Surv(gibbs_sampler, time = c(0:70))

plot(spline_many, t = c(0:70), xlim = c(0,70))
lines(y = km_update$surv,km_update$time, col = "blue" )
lines(x  = c(0:70),rowMeans(Surv_gibbs), col = "green")


rm(TA396_DT_OS_Initial_all_mods,TA396_DT_OS_boot, max_prob_chng) # Delete unneeded objects to save memory 

## 4.2 PFS ----


### Dabrafenib + Trametinib (DT) ----
#https://www.nejm.org/doi/pdf/10.1056/NEJMoa1904059?articleTools=true


TA396_DT_PFS_Initial_df <- rbind(`TA396_D+T_PFS_Initial_D`, `TA396_D+T_PFS_Initial_V`)


TA396_DT_PFS_Initial_model  <- collapsing.model(TA396_DT_PFS_Initial_df,
                                       n.iter = n.iter,
                                       burn_in = burn_in,
                                       n.chains = n.chains,
                                       timescale = "months",
                                       max.num.breaks = 15,
                                       seed.val = seed.val)

TA396_DT_PFS_Initial_posterior_chng <- plot.pos_changepoint(TA396_DT_PFS_Initial_model, breaks = seq(from = 0, to = 25, by = 5))+theme_light()

TA396_DT_PFS_Initial_plot <- plot(TA396_DT_PFS_Initial_model, max_predict  =7*12)

max_prob_chng <-   as.numeric(names(which.max(TA396_DT_PFS_Initial_model$prob.changepoint)))

TA396_DT_PFS_all_mods <- compare.surv.mods(TA396_DT_PFS_Initial_model, max_predict = 7*12,
                                           n.iter.jags = n.iter.jags,
                                           n.thin.jags = 1,
                                           n.burnin.jags = n.burnin.jags, 
                                           chng.num = max_prob_chng,
                                           gmp_haz_df =gmp_haz_df_TA396,
                                           col_km = "grey",
                                           final_chng_only=T)

TA396_DT_PFS_boot <- PiecewiseChangepoint:::compare_boot_sims(mod_parametric_orig = TA396_DT_PFS_all_mods,
                                      follow_up_data = `TA396_D+T_PFS_Update`)

write.xlsx(TA396_DT_PFS_boot, paste0(path,"pub_plots_tabs/TA396_DT_PFS_AUC.xlsx"))

TA396_DT_PFS_gg_final <-compare_dco(TA396_DT_PFS_all_mods, TA396_DT_PFS_Initial_df,`TA396_D+T_PFS_Update`)

TA396_DT_PFS_gg_final <- TA396_DT_PFS_gg_final+ xlab("Time (Months)")
ggsave(filename =paste0(path,"pub_plots_tabs/TA396_DT_PFS_Surv.png"),plot = TA396_DT_PFS_gg_final , height = 5.59, width = 7)




# Revewer 2 comment: Add Many knots/change-points
km_update <- survfit(Surv(time, status)~1, data = `TA396_D+T_PFS_Update`)
spline_many <- flexsurvspline(Surv(time, status)~1, data =TA396_DT_PFS_Initial_df, k = 5)

gibbs_sampler<- PiecewiseChangepoint:::gibbs.sampler(df =TA396_DT_PFS_Initial_df, num.breaks = 2, n.iter =5000,
                                                     burn_in = 100, num.chains = 2, 
                                                     n.thin = 1, alpha.hyper = 1, beta.hyper = 1)

colMeans(gibbs_sampler$changepoint)
max(TA396_DT_PFS_Initial_df$time)
exp(spline_many$knots)
Surv_gibbs <- PiecewiseChangepoint:::get_Surv(gibbs_sampler, time = c(0:70))

plot(spline_many, t = c(0:70), xlim = c(0,70))
lines(y = km_update$surv,km_update$time, col = "blue" )
lines(x  = c(0:70),rowMeans(Surv_gibbs), col = "green")

rm(TA396_DT_PFS_all_mods,TA396_DT_PFS_boot, max_prob_chng) # Delete unneeded objects to save memory 



### Comparator ----

TA396_COMP_PFS_Initial <- rbind(TA396_V_PFS_Initial, TA396_D_PFS_Initial)


TA396_COMP_PFS_model <- collapsing.model(TA396_COMP_PFS_Initial,
                                         n.iter = n.iter,
                                         burn_in = n.iter/2,
                                         n.chains = 2,
                                         timescale = "months",
                                         seed.val = seed.val,
                                         lambda.prior = 2)

# 5 TA428  ----

## 5.1 OS ----

TA428_PEM_OS_Initial_model <- collapsing.model(TA428_PEM_OS_Initial,
                                              n.iter = n.iter, # Change to 20750
                                              burn_in = burn_in,
                                              n.chains = n.chains,
                                              timescale = "months",
                                              seed.val = seed.val)

# Plot (and save) the posterior distribution of the changepoint
TA428_PEM_OS_Initial_posterior_chng <- plot.pos_changepoint(TA428_PEM_OS_Initial_model,
                                                           breaks = seq(from = 0, to = 25, by = 5))+theme_light()

max_prob_chng <-   as.numeric(names(which.max(TA428_PEM_OS_Initial_model$prob.changepoint)))

# Compare survival models

age_baseline_TA428 <- round(63.000, digits = 0)
prop_male_TA428 <- 0.616

Conditional_Death_df_temp <- Conditional_Death_df
Conditional_Death_df_temp[, "mix_prob"] <- Conditional_Death_df_temp[,2]*prop_male_TA428 + Conditional_Death_df_temp[,3]*(1-prop_male_TA428)
Conditional_Death_df_temp <- Conditional_Death_df_temp %>% filter(age >= age_baseline_TA428 & age <= time_horizon)

time_factor_TA428 <- 12

Conditional_Death_df_temp$mix_haz <- -log(1-Conditional_Death_df_temp$mix_prob)/time_factor_TA428

gmp_haz_vec_TA428 = rep(Conditional_Death_df_temp$mix_haz,each = time_factor_TA428)
gmp_haz_df_TA428 <- data.frame(time = 1:length(gmp_haz_vec_TA428), hazard = gmp_haz_vec_TA428)
#stroke = 2 change stroke
TA428_PEM_OS_Initial_all_mods <- compare.surv.mods(TA428_PEM_OS_Initial_model, max_predict = 7*12,
                                                  n.iter.jags = n.iter.jags, #Should be 5000
                                                  n.thin.jags = 1,
                                                  n.burnin.jags = n.burnin.jags,
                                                  chng.num = max_prob_chng,
                                                  gmp_haz_df =gmp_haz_df_TA428,
                                                  col_km = "grey",
                                                  km_risk =NULL)

# Compare with later data cut-off (TA396_D+T_OS_Update)
TA428_PEM_OS_gg_final <-compare_dco(TA428_PEM_OS_Initial_all_mods, TA428_PEM_OS_Initial, TA428_PEM_OS_Update,
                                    km_risk =NULL)

TA428_PEM_OS_gg_final <- TA428_PEM_OS_gg_final+ xlab("Time (Months)")
ggsave(filename =paste0(path,"pub_plots_tabs/TA428_PEM_OS_Surv.png"),plot = TA428_PEM_OS_gg_final , height = 5.59, width = 7)

#Compare the RMST time using the bootstrap of the KM curve
TA428_PEM_OS_boot <- PiecewiseChangepoint:::compare_boot_sims(mod_parametric_orig = TA428_PEM_OS_Initial_all_mods,
                                      follow_up_data = TA428_PEM_OS_Update)

write.xlsx(TA428_PEM_OS_boot, paste0(path,"pub_plots_tabs/TA428_PEM_OS_AUC.xlsx"))


rm(TA428_PEM_OS_Initial_all_mods,TA428_PEM_OS_boot, max_prob_chng) # Delete unneeded objects to save memory 

# 6 TA447 ----

## 6.1 TA447 ----

TA447_PEM_OS_Initial_model <- collapsing.model(TA447_PEM_OS_Initial,
                                            n.iter = n.iter,
                                            burn_in = burn_in,
                                            n.chains = n.chains,
                                            timescale = "months",
                                            seed.val = seed.val)
max_prob_chng <-   as.numeric(names(which.max(TA447_PEM_OS_Initial_model$prob.changepoint)))

#Does impact the survival probability give some timepoints

plt_val <- plot(TA447_PEM_OS_Initial_model, max_predict = 72, chng.num = 0)
plt_val2 <- plot(TA447_PEM_OS_Initial_model, max_predict = 72, chng.num = 1)

#Reviewer response query

times_validate <- c(24,36,48, 60)
TA447_PEM_OS_Surv_nochng <- get_Surv(TA447_PEM_OS_Initial_model, chng.num = 0,time = times_validate)
TA447_PEM_OS_Surv_1chng <- get_Surv(TA447_PEM_OS_Initial_model, chng.num = 1,time = times_validate)
apply(TA447_PEM_OS_Surv_nochng,1, quantile, probs = c(0.025,0.5,0.975))
apply(TA447_PEM_OS_Surv_1chng,1, quantile, probs = c(0.025,0.5,0.975))

age_baseline_TA447 <- round(65, digits = 0)
prop_male_TA447 <- 0.646

Conditional_Death_df_temp <- Conditional_Death_df
Conditional_Death_df_temp[, "mix_prob"] <- Conditional_Death_df_temp[,2]*prop_male_TA447 + Conditional_Death_df_temp[,3]*(1-prop_male_TA447)
Conditional_Death_df_temp <- Conditional_Death_df_temp %>% filter(age >= age_baseline_TA447 & age <= time_horizon)

time_factor_TA447 <- 12

Conditional_Death_df_temp$mix_haz <- -log(1-Conditional_Death_df_temp$mix_prob)/time_factor_TA447

gmp_haz_vec_TA447 = rep(Conditional_Death_df_temp$mix_haz,each = time_factor_TA447)
gmp_haz_df_TA447 <- data.frame(time = 1:length(gmp_haz_vec_TA447), hazard = gmp_haz_vec_TA447)


TA447_PEM_OS_Initial_all_mods <- compare.surv.mods(TA447_PEM_OS_Initial_model, max_predict = 7*12,
                                                n.iter.jags = n.iter.jags,
                                                n.thin.jags = 1,
                                                n.burnin.jags = n.burnin.jags,
                                                chng.num = max_prob_chng,
                                                gmp_haz_df =gmp_haz_df_TA447,
                                                col_km = "grey")

TA447_PEM_OS_gg_final <-compare_dco(TA447_PEM_OS_Initial_all_mods, TA447_PEM_OS_Initial, TA447_PEM_OS_Update)
TA447_PEM_OS_gg_final <- TA447_PEM_OS_gg_final+ xlab("Time (Months)")
ggsave(filename =paste0(path,"pub_plots_tabs/TA447_PEM_OS_Surv.png"),plot = TA447_PEM_OS_gg_final, height = 5.59, width = 7)


TA447_PEM_OS_boot <- PiecewiseChangepoint:::compare_boot_sims(mod_parametric_orig = TA447_PEM_OS_Initial_all_mods,
                                                              follow_up_data =TA447_PEM_OS_Update)

write.xlsx(TA447_PEM_OS_boot, paste0(path,"pub_plots_tabs/TA447_PEM_OS_AUC.xlsx"))

# 7  Extract Results for Table 1----

mod.names <- c("TA396_DT_OS_Initial_model",
               "TA396_DT_PFS_Initial_model",
               "TA396_COMP_PFS_model",
               "TA269_V_PFS_model",
               "TA269_D_PFS_model",
               "TA268_GP100_OS_model",
               "TA268_Ipi_OS_model",
               "TA268_Ipi_PFS_model",
               "TA268_GP100_PFS_model",
               "TA347_Doce_OS_model",
               "TA347_N_OS_model",
               "TA347_Doce_PFS_model",
               "TA347_N_PFS_model",
               "TA428_PEM_OS_Initial_model", 
               "TA447_PEM_OS_Initial_model")


TA_vec <- c("TA396","TA396", "TA396",
             "TA269","TA269",
            "TA268","TA268",
            "TA268","TA268",
            "TA347", "TA347","TA347", "TA347",
            "TA428",
            "TA447")

outcome <- c("OS", rep("PFS",4), rep("OS",2), rep("PFS",2),  rep("OS",2),rep("PFS",2), "OS", "OS")
arm <- c("DT", "DT", "Comp","BRAF_V",
         "BRAF_D", "GP100", "IPI", "IPI", "GP100",
         "D","N","D","N", "Pembro", "Pembro" )

#Custom function - extracts mean changepoint value

#undebug(summary.changepoint_mod)
summary.changepoint_mod <- function(object, chng.num = NULL){
  
  if (is.null(chng.num)) {
    chng.prob <- as.numeric(names(which.max(object$prob.changepoint)))
  } else {
    chng.prob <- chng.num
  }
  
  
  k <- object$k.stacked
  lambda <- object$lambda
  changepoint <- object$changepoint
  num.changepoints <- unlist(apply(k, 1, function(x) {
    length(na.omit(x))
  }))
  if(chng.prob == 0){
    changepoint_curr <- NULL
  }else{
    
    changepoint_curr <- data.frame(changepoint[which(num.changepoints == chng.prob), 1:chng.prob])
    
  }
  lambda_curr <- data.frame(lambda[which(num.changepoints == chng.prob), 1:(chng.prob + 1)])
  
  if(chng.prob == 0){
    
    return_vec <- c(NA,NA,NA,
                    round(tail(colMeans(lambda_curr), n = 1), 3))
    
  }else{
    quantile_vec <- apply(changepoint_curr, 2, quantile, probs = c(0.025,0.975))
    
    return_vec <- c(round(tail(colMeans(changepoint_curr), n = 1),2),
                    round(quantile_vec[,ncol(quantile_vec)],2),
                    round(tail(colMeans(lambda_curr), n = 1), 3))
    
    
  }
  return(return_vec)
  
}
output_mat <- matrix(ncol = 4, nrow = length(arm))
for(i in 1:nrow(output_mat)){
  output_mat[i, ] <- summary.changepoint_mod(get(mod.names[i]))
  
}

colnames(output_mat) <- c("final_changepoint_months","final_changepoint_months_0.025",
                          "final_changepoint_months_0.975", "final_lambda")

output <- data.frame(TA_vec, outcome, arm , output_mat) %>%
  mutate(final_changepoint_months =  ifelse(TA_vec == "TA268", final_changepoint_months*12,
                                            final_changepoint_months),
         final_changepoint_months_0.025 =  ifelse(TA_vec == "TA268",
                                                  (final_changepoint_months_0.025*12),
                                                  final_changepoint_months_0.025),
         final_changepoint_months_0.975 =  ifelse(TA_vec == "TA268",
                                                  (final_changepoint_months_0.975*12),
                                                  final_changepoint_months_0.975),
         final_lambda = ifelse(TA_vec == "TA268", final_lambda/12,
                               final_lambda)) %>%
  mutate(final_changepoint_months =  ifelse(TA_vec == "TA347"& outcome == "PFS",
                                            (final_changepoint_months*12)/365.25,
                                            final_changepoint_months),
         final_changepoint_months_0.025 =  ifelse(TA_vec == "TA347"& outcome == "PFS",
                                            (final_changepoint_months_0.025*12)/365.25,
                                            final_changepoint_months_0.025),
         final_changepoint_months_0.975 =  ifelse(TA_vec == "TA347"& outcome == "PFS",
                                            (final_changepoint_months_0.975*12)/365.25,
                                            final_changepoint_months_0.975),
         final_lambda = ifelse(TA_vec == "TA347", (final_lambda*365.25)/12,
                               final_lambda))

#Output Table 1 of results in Publication.
output_final <- output %>%
  arrange(match(TA_vec, c("TA268","TA269","TA347","TA396", "TA428")) ,outcome, arm) %>%
  mutate(final_changepoint_months = round(final_changepoint_months, digits = 1),
         final_changepoint_months_0.025  = round(final_changepoint_months_0.025 , digits = 1),
         final_changepoint_months_0.975   = round(final_changepoint_months_0.975  , digits = 1))


output_final %>%filter(outcome == "PFS")
output_final %>%filter(outcome == "OS")

write.xlsx(output, paste0(path, "pub_plots_tabs/Piecewise_TA_results.xlsx"))


# 8 Plot Cum Hazard and Hazard Functions for TA428 ----
list.of.packages <- c("Epi","bshazard","muhaz","ggpubr") 

#install.packages(list.of.packages)

library("Epi")
library("bshazard")
library("muhaz")
library("ggpubr")


TA428_PEM_OS_Update_model <- collapsing.model(TA428_PEM_OS_Update,n.iter =20750,
                                        burn_in = 750,
                                        n.chains = 2,
                                        timescale = "months")

data_vec <- c("TA428_PEM_OS_Initial","TA428_PEM_OS_Update" )
model_vec <- c("TA428_PEM_OS_Initial_model","TA428_PEM_OS_Update_model")
ggtitle_vec <- c("Hazard Functions for KEYNOTE-010 trial data available at TA428", "Hazard Functions for 5 year update of KEYNOTE-010 trial data")
cum_haz_vec <- c("Cum. Hazard Function for KEYNOTE-010 trial data available at TA428", "Cum. Hazard Functions for 5 year update of KEYNOTE-010 trial data")

for(i in 1:2){
  
  
  data_used <- get(paste0(data_vec[i]))
  model_used <- get(paste0(model_vec[i]))
  data_used$enter <- 0
  max_time <- max(data_used$time)
  model_best_fit <- as.numeric(names(which.max(model_used$prob.changepoint)))
  
  
  bshazard.fit <- bshazard(formula = Surv(enter, time,
                                          status) ~ 1,
                           data = data_used,
                           nbin = length(unique(data_used[which(data_used$status == 1),"time"])))
  
  
  cumhaz_plot <- survminer::ggsurvplot(survfit(formula = Surv(enter, time,
                                                              status) ~ 1,
                                               data = data_used),fun = "cumhaz")
  cumhaz_plot_df <- cumhaz_plot$plot$data
  
  if(i == 1){
    time_index <- 1
  }else{
    time_index <- 5
  }
  
  cumhaz_plot <- ggplot(cumhaz_plot_df, aes(x = time, y = surv))+
    geom_step(colour="red")+
    ylab("Cumulative hazard")+
    xlab("Time (Months)")+
    scale_x_continuous(expand = c(0, 0), limits = c(0, round(max_time,0)), breaks=seq(0, round(max_time,0), by = time_index))+
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=strata),
                alpha=0.3, colour=NA)+ theme_bw()
  assign(paste0("cumhaz_plot",i),cumhaz_plot +ggtitle(cum_haz_vec[i]) )
  

  bshazard.df  <- data.frame(time = bshazard.fit[["time"]],
                             hazard = bshazard.fit[["hazard"]],
                             upper.ci = bshazard.fit[["upper.ci"]],
                             lower.ci = bshazard.fit[["lower.ci"]])
  
  
  
  harzard.plt.mod <- function(time.vec, cens.vec, Outcome = "Survival",
                              lrg.haz.int = 3,  bw.method = "local"){
    
    #Plot the hazards
    max.time <- max(time.vec)
    
    result.hazard.pe.lrg <- pehaz(time.vec, cens.vec, width= lrg.haz.int, max.time=round(max_time,0))
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
  
  
  hazard.PFS.mod <-harzard.plt.mod(time.vec = data_used$time ,
                                   cens.vec =  data_used$status ,
                                   Outcome = "Death", lrg.haz.int = 3)
  
  df.haz <- data.frame(time = hazard.PFS.mod$Cuts[-1]-1,
                       hazard = c(hazard.PFS.mod$Hazard))
  
  
  
  haz_plot_collapsing <- plot(model_used, type = "hazard", max_predict = round(max_time,0),
                              chng.num = model_best_fit)
  haz_df <- haz_plot_collapsing$data %>% mutate(timepoints = timepoints,
                                                hazards = hazards)
  haz_df_summary <- haz_df %>% group_by(timepoints) %>% summarize(hazards.mean = mean(hazards),
                                                                  hazards.975 = quantile(hazards, 0.975),
                                                                  hazards.025 = quantile(hazards, 0.025))
  
  
  
  
  assign(paste0("plot",i),ggplot(bshazard.df[which(bshazard.df$time <= round(max_time,0)),], aes(x = time, y = hazard))+
           geom_line(aes(colour = "blue"), size = 1.3,linetype = "longdash")+
           geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill = "blue"), alpha = 0.2)+
           # geom_step(data = df.haz[which(df.haz$time <= round(max_time,0)),],  aes(time, hazard, colour = "black"), size = 1.1,linetype = "dotdash")+
           xlab("Time (Months)")+
           ylab("Hazard")+
           ggtitle(ggtitle_vec[i])+
           #geom_step(data = haz_df ,aes(timepoints, hazards,group = id), linetype = "dashed", alpha = 0.03, colour = "red")+
           geom_line(aes(timepoints, hazards.mean, colour = "purple"),data = haz_df_summary, size = 1.3)+
           ggplot2::geom_ribbon(aes(ymin = hazards.025, ymax = hazards.975, x = timepoints,fill = "purple"),data =haz_df_summary,
                                inherit.aes = F, alpha = 0.2)+
           #scale_y_continuous(breaks = seq(from = 0, to = 0.2, by =0.1))+
           scale_x_continuous(expand = c(0, 0), limits = c(0, round(max_time,0)), breaks= seq(0, round(max_time,0), by = time_index)) +
           scale_y_continuous(expand = c(0, 0), limits = c(0, NA), breaks = seq(from = 0, to = 0.2, by =0.01))+
           # coord_cartesian(ylim = c(0,.2),xlim = c(0,max_time))+
           scale_colour_manual(name = 'Hazard Functions',
                               values =c('black'='black','blue'='blue','purple'= 'purple'), labels = c('Interval hazard','bsharzard', 'PEM'))+
           scale_fill_identity(name = '95% Intervals', guide = 'legend',labels = c('bsharzard', "PEM")) +
           theme_light())
  
  
  
}

#Figure 2 in Publication
ggarrange( cumhaz_plot1,plot1, ncol = 1, labels = c("A", "B"), common.legend = T, legend = "bottom",
           heights = c(4, 4), nrow = 2, align = "v")
ggsave(paste0(path,"pub_plots_tabs/TA428_PEM_OS_Initial_hazards.png"), width = 15, height = 10)
