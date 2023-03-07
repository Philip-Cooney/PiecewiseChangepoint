setwd("C:\\Users\\phili\\OneDrive\\PhD\\R Codes\\Changepoint Analysis\\")

#setwd("C:/Users/cooneph1/Desktop/Simulation Study/")
source("Functions Reversible Jump chains final.R")
source("Functions Gibbs Markdown.R")




#Not Run
n_obs <- n_events_req <- 500
rate <- c(0.75, 0.5,0.25)
num.breaks <- length(rate)-1
t_change <- c(0.5,1)
max_time <- 2
df1 <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                                       num.breaks = num.breaks,rate = rate ,
                                       t_change = t_change, max_time = max_time)

res_gibbs1 <- gibbs.changepoint_chains_optim(df1, num.breaks = 1, n.iter = 10000, burn_in = 100, num.chains = 1)
res_gibbs2 <- gibbs.changepoint_chains_optim(df1, num.breaks = 2, n.iter = 10000, burn_in = 100, num.chains = 2)
res_gibbs3 <- gibbs.changepoint_chains_optim(df1, num.breaks = 3, n.iter = 10000, burn_in = 100, num.chains = 1)
as.numeric(res_gibbs3[["mean.marg.lik.log"]]/res_gibbs1[["mean.marg.lik.log"]])

plot(survival::survfit(Surv(time,status)~1,data = df1))

df1_new <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
                    num.breaks = num.breaks,rate = rate ,
                    t_change = t_change, max_time = Inf)

df_plt<- dplyr::bind_rows(df1, df1_new,.id = "id")

#plot(survival::survfit(Surv(time,status)~1,data = df1))
plot(survival::survfit(Surv(time,status)~id,data = df_plt))
time<- c(seq(from = 0, to  = max_predict, by = 0.1))

survival::survdiff(Surv(time,status)~id,data = df_plt)

test <- RJMC.piece(10750, df1,n.chains = 1,  burn_in = 750 , obs_true =df1$time_event,
                   max_predict = 20)


##Kidney

kidney <-survival::kidney
plot(survfit(Surv(time,status)~1,data = kidney))
kidney2 <- kidney[,c("time","status")]
kidney2$time <- kidney2$time/365

test <- RJMC.piece(10750, kidney2,n.chains = 1,  burn_in = 750,
                   max_predict = 2, max_plt = 2)

#test[["output"]]["2"]
test[["Log.Lik.models"]]
test[["prob.changepoint"]]
test[["output"]][["0"]]
test[["output"]][["1"]]
num.changepoints <-unlist(apply(test[["k.stacked"]],1, function(x){length(na.omit(x))}))
colMeans(test[["changepoint"]][which(num.changepoints ==1),],na.rm= T)
colMeans(test[["lambda"]][which(num.changepoints ==1),],na.rm= T)


#Stanford

data("stanford2")
stanford2$time <- stanford2$time/365
surv.stanford <- survfit(Surv(time,status)~1,data = stanford2)
plot(surv.stanford)

cum_haz_stan <- ggsurvplot(surv.stanford, risk.table = TRUE, fun = "cumhaz",break.time.by =1)
ggsave(file = "cum_haz_stan.png", print(cum_haz_stan))
ggsave("cum_haz_stan.png")

test <- RJMC.piece(10750, stanford2,n.chains = 1,  burn_in = 750,
                   max_predict = 15, max_plt = 15)
test[["output"]][["2"]]

stanford_test <-stanford2
for(i in 1:nrow(stanford_test)){
  if(stanford_test$time[i] >= 2){
    stanford_test$time[i] <- 2
    stanford_test$status[i] <-0
  }
}


test <- RJMC.piece(10750, stanford_test,n.chains = 1,  burn_in = 750,
                   max_predict = 15, max_plt = 15, 
                   obs_true = stanford2[,c("time","status")])
test[["prob.changepoint"]]
num.changepoints <-unlist(apply(test[["k.stacked"]],1, function(x){length(na.omit(x))}))
colMeans(test[["changepoint"]][which(num.changepoints ==2),],na.rm= T)
colMeans(test[["lambda"]][which(num.changepoints ==2),],na.rm= T)
test[["output"]][["2"]]
test[["Log.Lik.models"]]


#Staph Aureus
data("burn")
burn$time <- burn$T3
burn$status <- burn$D3
burn$time <-burn$time/30.25
plot(survfit(Surv(time,status)~1,data = burn))


test <- RJMC.piece(3000, burn,n.chains = 1,  burn_in = 750,
                   max_predict = 10, max_plt = 10)

test[["prob.changepoint"]]

test[["output"]][["3"]]


#Leukemia dataset from Brazil

df1 <-read.xlsx("Achar Leukemia Data.xlsx",1)
df1$time <- df1$time_adj

df1$time <-df1$time/365.25
plot(survfit(Surv(time,status)~1,data = df1))

plot(y = log(survfit(Surv(time,status)~1,data = df1)$surv), survfit(Surv(time,status)~1,data = df1)$time, type = "l")

test <- RJMC.piece(30750, df1,n.chains = 1,  burn_in = 750,
                   max_predict = 10, max_plt = 10, obs_true = df1, alpha.hyper = )

model_expo <- flexsurvreg(Surv(time,status)~1, data = df1, dist = "exp")
model_weibull <- update(model_expo, dist = "weibull")


test[["prob.changepoint"]]
#Check why the RJMC.piece does not find the MLE.


time_diffs <- df_recast(df1)
times <- exposure_death_alt(time_diffs,NA)
marg.lik.eval.log(times,1,1)


res_gibbs <- gibbs.changepoint_chains_optim(df1, num.breaks = 1, n.iter = 10000, burn_in = 100, num.chains = 1)
res_gibbs[["plot_Surv"]]

res_gibbs2 <- gibbs.changepoint_chains_optim(df1, num.breaks = 2, n.iter = 10000, burn_in = 100, num.chains = 1)

as.numeric(res_gibbs2[["mean.marg.lik.log"]]/res_gibbs[["mean.marg.lik.log"]])
res_gibbs2[["plot_Surv"]]
plt.save <- res_gibbs2[["plot_Surv"]]
ggsave("plt.save.png")
exp(1.8616)
#par(mfrow = c(1,1))

test <- gibbs.changepoint_chains_optim(df1, num.breaks = 2, n.iter = 10000, burn_in = 0, num.chains = 2)

summary(test[["stacked_df"]])

t(apply(test[["stacked_df"]],2, FUN = quantile, probs = c(0.5,0.025,0.975)))

test[["plot_Surv"]]
raftery.diag(test[[1]]$chain_1, r = 0.005) 
raftery.diag(test[[1]]$chain_1, r = 0.01) 
gelman.diag(test[[1]])

geweke.diag(test[[1]]$chain_1)

df <- df1
num.breaks = 2
n.iter = 100
burn_in = 10
num.chains = 2


res_gibbs3[["plot_Surv"]]
library(gridExtra)
grid.arrange(res_gibbs[["plot_Surv"]], res_gibbs2[["plot_Surv"]],res_gibbs3[["plot_Surv"]], nrow = 2)
res_gibbs3 <- gibbs.changepoint_chains_optim(df1, num.breaks = 3, n.iter = 10000, burn_in = 100, num.chains = 1)
as.numeric(res_gibbs3[["mean.marg.lik.log"]]/res_gibbs2[["mean.marg.lik.log"]])
res_gibbs3[["plot_Surv"]]

num.changepoints <-unlist(apply(test[["k.stacked"]],1, function(x){length(na.omit(x))}))
colMeans(test[["changepoint"]][which(num.changepoints ==1),],na.rm= T)
colMeans(test[["lambda"]][which(num.changepoints ==1),],na.rm= T)

1.442965 *365.25
test[["output"]][["1"]]
704/365
test[["output"]][["2"]]
test[["output"]][["3"]]
df1$enter <- 0

fit <- phreg(Surv(enter, time, status) ~ 1, data = df1, dist = "pch", cuts = c(1.47))

fit <- phreg(Surv(enter, time, status) ~ 1, data = df1, dist = "pch", cuts = c(1.92))

df1 <-read.xlsx("Hacettepe Study.xlsx",1)
plot(survfit(Surv(time,status)~1,data = df1))

df1$time <-df1$time/12

plot(survfit(Surv(time,status)~1,data = df1))

test <- RJMC.piece(10750, df1,n.chains = 1,  burn_in = 750,
                   max_predict = 15, max_plt = 15,lambda.prior = 1)


#Need to be able to account for multiple events at the same timepoint
test[["prob.changepoint"]]
test[["output"]][["2"]]
df1$enter <- 0
expos_df<- exposure_death(df1,c( 2.935056, 3.015192, 6.786435, 7.004332, 8.897414, 9.090412))

expos_df[,1]/expos_df[,2]

num.changepoints <-unlist(apply(test[["k.stacked"]],1, function(x){length(na.omit(x))}))
colMeans(test[["changepoint"]][which(num.changepoints ==6),],na.rm= T)
colMeans(test[["lambda"]][which(num.changepoints ==6),],na.rm= T)

#file:///C:/Users/phili/OneDrive/PhD/Papers%20for%20PhD/Change%20Point%20analysis/A%20Constant%20Hazard%20Function%20Model%20with%20a%20change%20point%20-%20Includes%20data[Achcar%201998].pdf

#Initialize the packages if running from a Novartis Machine.
#Latex
#http://pages.stat.wisc.edu/~jgillett/371/RStudio/RMarkdown.pdf
#https://www.calvin.edu/~rpruim/courses/s341/S17/from-class/MathinRmd.html
#http://www.math.mcgill.ca/yyang/regression/RMarkdown/example.html
##https://oeis.org/wiki/List_of_LaTeX_mathematical_symbols



#https://acsjournals.onlinelibrary.wiley.com/doi/epdf/10.1002/cncr.32162
url.path <- "http://merlot.stat.uconn.edu/~mhchen/survbook/dataset/e1690.missing.dat"
E1690.dat <- read.delim(url(url.path), header = T, sep="", skip=12, as.is=TRUE)
#Drop PFS events with time equal zero
E1690.dat <- E1690.dat[-which(E1690.dat$FAILTIME ==0),]

#Convert to the correct notation for survival objects
E1690.dat[which(E1690.dat$FAILCENS == 1),"FAILCENS"] <-0
E1690.dat[which(E1690.dat$FAILCENS == 2),"FAILCENS"] <-1
E1690.dat[which(E1690.dat$SURVCENS == 1),"SURVCENS"] <-0
E1690.dat[which(E1690.dat$SURVCENS == 2),"SURVCENS"] <-1

#Create survival objects
fit.OS <- survfit(Surv(SURVTIME, SURVCENS)~TRT, 
                  data = E1690.dat)#[which(E1690.dat$STUDY == "1684"),])
fit.PFS <- survfit(Surv(FAILTIME, FAILCENS)~TRT, 
                   data = E1690.dat)#[which(E1690.dat$STUDY == "1684"),])

#Plot of Kaplan Meiers

OS.obs1 <- as.numeric(fit.OS$strata[1])
PFS.obs1 <- as.numeric(fit.PFS$strata[1])


plot(fit.OS$time[1:OS.obs1],fit.OS$surv[1:OS.obs1],
     col = "red", lty = 1, type = "l", ylim = c(0,1), main = "OS and PFS Survival",
     xlab = "Time in years", ylab = "Survival")
lines(fit.OS$time[OS.obs1+1:length(fit.OS$time)],fit.OS$surv[OS.obs1+1:length(fit.OS$time)],
      col = "blue", lty = 1)
lines(fit.PFS$time[1:PFS.obs1],fit.PFS$surv[1:PFS.obs1],
      col = "darkred", lty = 2, type = "l")
lines(fit.PFS$time[PFS.obs1+1:length(fit.PFS$time)],fit.PFS$surv[PFS.obs1+1:length(fit.PFS$time)],
      col = "darkblue", lty = 2)
legend("topright", legend = c("OS Trt 1 (OBS)", "OS Trt 2 (INF)", "PFS Trt 1 (OBS)", "PFS Trt 2 (INF)"),
       col = c("red", "blue", "darkred", "darkblue"), lty = c(1,1,2,2), cex = 0.8)


trt.arm <- E1690.dat[which(E1690.dat$TRT == 2),]
#trt.arm <- trt.arm %>% dplyr::rename(time = FAILTIME, status = FAILCENS)
trt.arm <- trt.arm %>% dplyr::rename(time = SURVTIME, status = SURVCENS)

fit.km <-survfit(Surv(time, status)~TRT, 
                 data = trt.arm)
ggsurvplot(fit.km,risk.table = T)





test <- RJMC.piece(10000, trt.arm,n.chains = 1,  burn_in = 750 , 
                   max_predict = 20, max_plt = 20, lambda = 1 )

test[["prob.changepoint"]]
#test[["output"]][["5"]]
test[["output"]][["3"]]
test[["output"]][["2"]]
num.changepoints <-unlist(apply(test[["k.stacked"]],1, function(x){length(na.omit(x))}))
colMeans(test[["changepoint"]][which(num.changepoints ==3),],na.rm= T)
colMeans(test[["lambda"]][which(num.changepoints ==3),],na.rm= T)




trt.arm <- E1690.dat[which(E1690.dat$TRT == 2),]
#trt.arm <- trt.arm %>% dplyr::rename(time = FAILTIME, status = FAILCENS)
trt.arm <- trt.arm %>% dplyr::rename(time = SURVTIME, status = SURVCENS)

fit.km <-survfit(Surv(time, status)~TRT, 
                 data = trt.arm)
ggsurvplot(fit.km,risk.table = T)


test <- RJMC.piece(10000, trt.arm,n.chains = 1,  burn_in = 750 , 
                   max_predict = 20, max_plt = 20 )
test[["prob.changepoint"]]

test[["output"]][["3"]]
test[["output"]][["2"]]
num.changepoints <-unlist(apply(test[["k.stacked"]],1, function(x){length(na.omit(x))}))
colMeans(test[["changepoint"]][which(num.changepoints ==2),],na.rm= T)
colMeans(test[["lambda"]][which(num.changepoints ==2),],na.rm= T)

mean(as.numeric(trt.arm$AGE), na.rm = T)



