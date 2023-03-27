devtools::install_github("Anon19820/PiecewiseChangepoint", lib = "~/R-packages/")

library("survminer")
library("survival")
library("xlsx")
library("dplyr")
library("flexsurv")
library("PiecewiseChangepoint")
pathway <- "~/PhD/KM_Piecewise_Major_Review_final/Readme Update/"

library("PiecewiseChangepoint",lib.loc = "~/R-packages/")
## basic example code

set.seed(123)
n_obs =300
n_events_req=300
max_time =  24 # months

rate = c(0.75,0.25)/12 # we want to report on months
t_change =12 # change-point at 12 months

df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                   num.breaks = length(t_change),rate = rate ,
                   t_change = t_change, max_time = max_time)

fit <- survfit(Surv(time, status) ~ 1, data = df)
# Drawing curves
ggsurvplot(fit, palette = "#2E9FDF")




ggsurvplot(fit, palette = "#2E9FDF", fun = "cumhaz")

Collapsing_Model <- collapsing.model(df,
                                     n.iter = 2070,
                                     burn_in = 75,
                                     n.chains = 2,
                                     timescale = "months")

print(Collapsing_Model)                                    

plot(Collapsing_Model, max_predict = 60, chng.num = 1)+xlab("Time (Months)")
plot(Collapsing_Model, type = "hazard")+xlab("Time (Months)")+ylab("Hazards")

plot(Collapsing_Model, max_predict = 5, chng.num = 2)

age_baseline_example <- 55
Conditional_Death_df <- read.xlsx(paste0(pathway, "Examples/Conditional_Death_UK.xlsx"), 1) %>% 
                          filter(age >=age_baseline_example)
head(Conditional_Death_df)
prop_male <- 0.5
time_horizon <- 100 
time_factor <- 12
Conditional_Death_df_temp <- Conditional_Death_df
Conditional_Death_df_temp[, "mix_prob"] <- Conditional_Death_df_temp[,2]*prop_male + Conditional_Death_df_temp[,3]*(1-prop_male)
Conditional_Death_df_temp <- Conditional_Death_df_temp %>% filter(age >= age_baseline_example & age <= time_horizon)

Conditional_Death_df_temp$mix_haz <- -log(1-Conditional_Death_df_temp$mix_prob)/time_factor

gmp_haz_vec_example = rep(Conditional_Death_df_temp$mix_haz,each = time_factor)
#We now have the hazard at each timepoint
gmp_haz_df_example <- data.frame(time = 1:length(gmp_haz_vec_example), hazard = gmp_haz_vec_example)


#This can take a number of minutes 
undebug(compare.surv.mods)
undebug(PiecewiseChangepoint:::fit_surv_models)
library("rstan")
rps.stan <- "// Royston-Parmar splines model\n\nfunctions {\n  real rps_lpdf(vector t, vector d, vector gamma, matrix B, matrix DB, vector linpred) {\n    // t = vector of observed times\n    // d = event indicator (=1 if event happened and 0 if censored)\n    // gamma = M+2 vector of coefficients for the flexible part\n    // B = matrix of basis\n    // DB = matrix of derivatives for the basis\n    // linpred = fixed effect part\n    vector[num_elements(t)] eta;\n    vector[num_elements(t)] eta_prime;\n    vector[num_elements(t)] log_lik;\n    real lprob;\n    \n    eta = B*gamma + linpred;\n    eta_prime = DB*gamma;\n    log_lik = d .* (-log(t) + log(eta_prime) + eta) - exp(eta);\n    lprob = sum(log_lik);\n    return lprob;\n  }\n  \n  real Sind( vector gamma, row_vector B, real linpred) {\n    // t = vector of observed times\n    // gamma = M+2 vector of coefficients for the flexible part\n    // B = row_vector of basis\n    // linpred = fixed effect part\n    real eta;\n    real Sind_rtn;\n    \n    eta = B*gamma + linpred;\n    Sind_rtn = exp(-exp(eta));\n    return Sind_rtn;\n  }\n  \n  \n\n}\n\ndata {\n  int<lower=1> n;                   // number of observations\n  int<lower=0> M;                   // number of internal knots for the splines model\n  int<lower=1> H;                   // number of covariates in the (time-independent) linear predictor\n  vector<lower=0>[n] t;             // observed times (including censored values)\n  vector<lower=0,upper=1>[n] d;     // censoring indicator: 1 if fully observed, 0 if censored\n  matrix[n,H] X;                    // matrix of covariates for the (time-independent) linear predictor\n  matrix[n,M+2] B;                  // matrix with basis\n  matrix[n,M+2] DB;                 // matrix with derivatives of the basis\n  vector[H] mu_beta;                // mean of the covariates coefficients\n  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients\n  vector[M+2] mu_gamma;             // mean of the splines coefficients\n  vector<lower=0>[M+2] sigma_gamma; // sd of the splines coefficients\n  \n}\n\n\nparameters {\n  vector[M+2] gamma;\n  vector[H] beta;\n}\n\n\ntransformed parameters{\n  vector[n] linpred;\n  vector[n] mu;\n\n  linpred = X*beta;\n  for (i in 1:n) {\n    mu[i] = linpred[i];\n  }\n\n}\n\nmodel {\n  // Priors\n  gamma ~ normal(mu_gamma,sigma_gamma);\n  beta ~ normal(mu_beta,sigma_beta);\n  \n  // Data model\n  t ~ rps(d,gamma,B,DB,X*beta);\n  \n}"
rps.stan_mod <- rstan::stan_model(model_code = rps.stan)

mod_comp <- compare.surv.mods(Collapsing_Model, max_predict = 100,
                                                   n.iter.jags = 5000, #Run JAGS/Stan for 5000 samples
                                                   n.thin.jags = 1,
                                                   n.burnin.jags = 500,
                                                   chng.num = 1, #Using results from 1 change-point PEM
                                                   gmp_haz_df =gmp_haz_df_example) #GPM dataset 


#Returns a dataframe with the model fit results
mod_comp$mod.comp

mod_comp$plot_Surv_all




Bullshit hack to fix comparesurv mods
df_km <- data.frame(Surv = c(0), t_pred = 0)

plot_surv + geom_line(df_km , 
                      mapping = aes_string(y = "Surv", x = "t_pred", colour = "'KM curve'"), inherit.aes = F)


plot_surv <- plot_surv + geom_line(df_km , 
                                   mapping = aes_string(y = "Surv", x = "t_pred", colour = "'Piecewise Expo'"), inherit.aes = F)
