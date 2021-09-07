
remove.packages("PiecewiseChangepoint")
devtools::load_all()

set.seed(123)
n_obs =300
n_events_req=300
max_time =  2

rate = c(0.75,0.25)
t_change =1


#hesim::rpwexp
df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                   num.breaks = length(t_change),rate = rate ,
                   t_change = t_change, max_time = max_time)


Collapsing_Model <- collapsing.model(df,
                                     n.iter = 20750,
                                     burn_in = 750,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1)


fit <- survfit(Surv(time, status) ~ 1, data = df)
# Drawing curves
ggsurvplot(fit, palette = "#2E9FDF")
ggsurvplot(fit, palette = "#2E9FDF", fun = "cumhaz")


saveRDS(Collapsing_Model, file = "Examples/Collapsing_Model.rds")
# Restore the object
#readRDS(file = "Collapsing_Model.rds")
chain.mixing(Collapsing.Model)
print(Collapsing_Model, chng.num = 1)
plot(Collapsing_Model, chng.num = 1, add.post = T)
plot(Collapsing_Model, chng.num = 1,type = "hazard")


mod_comp <-compare.surv.mods(Collapsing_Model)
saveRDS(mod_comp$mod.comp, file = "Examples/mod_comp.rds")

df_true <- df
df_true$time <- df_true$time_event
df_true$status <- 1
df_true <- df_true %>% mutate(time = ifelse(time >10, 10, time))
df_true <- df_true %>% mutate(status = ifelse(time >=10, 0, status))
plot_surv_true <- add_km(mod_comp$plot_Surv_all, df_true, colour = "black")


ggsave(file = "Examples/plt_Survival.png")

saveRDS(mod_comp$plot_Surv_all, file = "Examples/plot_Surv_all.rds")


mod_comp$plot_Surv_all
#Test to ensure the surv probability is correct

test<- get_Surv(Collapsing_Model, max_predict = 2)
plot(x = as.numeric(rownames(test)), y = test[,sample(ncol(test),1)])


### Test versus Weibull Example
set.seed(123)
times <- rweibull(n = 300, shape =1.2, scale  = 1.3)
df.weibull <- data.frame(time = times)

df.weibull <- df.weibull %>% mutate(time = ifelse(time >=1, 1, time)) %>%
                             mutate(status = ifelse(time ==1,  0,1))

fit <- survfit(Surv(time, status) ~ 1, data = df.weibull)
# Drawing curves
ggsurvplot(fit, palette = "#2E9FDF")

Collapsing_Model_weibull <- collapsing.model(df.weibull,
                                     n.iter = 20750,
                                     burn_in = 750,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1)


#plot(Collapsing_Model_weibull, add.post = T)

mod_comp_webiull <-compare.surv.mods(Collapsing_Model_weibull)
saveRDS(mod_comp_webiull$mod.comp, file = "Examples/mod_comp_weibull.rds")
saveRDS(mod_comp_webiull$plot_Surv_all, file = "Examples/plot_Surv_all_weibull.rds")


mod_comp_webiull$mod.comp
fit <- survfit(Surv(time, status) ~ 1, data = df)

mod_comp_webiull$plot_Surv_all

# Drawing curves
ggsurvplot(fit, palette = "#2E9FDF")
ggsurvplot(fit, palette = "#2E9FDF", fun = "cumhaz")


# Scale to Years

set.seed(123)
n_obs =300
n_events_req=300
max_time =  2

rate = c(0.75,0.25)
t_change =1


#hesim::rpwexp
df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                   num.breaks = length(t_change),rate = rate ,
                   t_change = t_change, max_time = max_time)

df$time <- df$time*365
df$time_event <- df$time_event*365
set.seed(212)
Collapsing_Model_prior <- collapsing.model(df,
                                     n.iter = 20750,
                                     burn_in = 750,
                                     n.chains = 2,
                                     alpha.hyper = 1,
                                     beta.hyper1 = 1,
                                     beta.hyper2 = 1/365)


plot(Collapsing_Model_prior, max_predict = 10*365)

