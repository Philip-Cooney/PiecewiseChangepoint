path <- "C:/Users/phili/OneDrive/PhD/R codes/Changepoint Analysis/"

source(paste0(path, "Functions Gibbs Markdown changed December.R"))
source(paste0(path, "Functions Reversible Jump chains final hyper.R"))


#Other Datasets

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("curatedBreastData")
library(survminer)


#https://bioconnector.github.io/workshops/r-survival.html

# Load the bioconductor installer. 
# Try http:// if https:// doesn't work.

#BiocManager::install("RTCGA", lib = "C:/Users/phili/Documents/R/win-library/4.0")
# Install the main RTCGA package
library("RTCGA")
#BiocManager::install("RTCGA.clinical")
# Create the clinical data


#https://bioconnector.github.io/workshops/r-survival.html#tcga
#Survival data.
library(RTCGA.clinical)
library(survival)

clinfinal <- survivalTCGA(GBM.clinical,
                          extract.cols="admin.disease_code")

clinfinal$time <- clinfinal$times#/365
clinfinal[which(clinfinal$time ==0),"time"] <- 0.001 

clinfinal$status <- clinfinal$patient.vital_status

surv.fit_glio <- survfit(Surv(time, status)~1, data=clinfinal)

ggsurvplot(surv.fit_glio)


clin.RJMCMC <- collapsing.model(n.iter = 20750, df = clinfinal,n.chains = 1,  burn_in = 750,
                          max_predict = 4000,lambda.prior = 1, alpha.hyper =  1,beta.hyper1 = 3650, beta.hyper2 = 10,interval = 50)

#saveRDS(clin.RJMCMC, paste0(path_biometrics,"Glioblastoma Example\\Glio Collapsing.rds"))
clin.RJMCMC <- readRDS( paste0(path_biometrics,"Glioblastoma Example\\Glio Collapsing.rds"))

changepoint_num <- 2

chang_num_vec <- apply(clin.RJMCMC[["k.stacked"]],1,function(x){length(na.omit(x))})

#summary(clin.RJMCMC[["changepoint"]][which(chang_num_vec == changepoint_num),])
#summary(clin.RJMCMC[["lambda"]][which(chang_num_vec == changepoint_num),])

apply(clin.RJMCMC[["changepoint"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)
apply(clin.RJMCMC[["changepoint"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)
apply(clin.RJMCMC[["changepoint"]][which(chang_num_vec == changepoint_num),],2,quantile,na.rm =T, probs = c(0.025,0.5,0.975))

apply(clin.RJMCMC[["lambda"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)
apply(clin.RJMCMC[["lambda"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)
apply(clin.RJMCMC[["lambda"]][which(chang_num_vec == changepoint_num),],2,quantile,na.rm =T, probs = c(0.025,0.5,0.975))

apply(clin.RJMCMC[["lambda.df_mu"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)
apply(clin.RJMCMC[["lambda.df_mu"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)
apply(clin.RJMCMC[["lambda.df_mu"]][which(chang_num_vec == changepoint_num),],2,quantile,na.rm =T, probs = c(0.025,0.5,0.975))

apply(sqrt(clin.RJMCMC[["lambda.df_var"]][which(chang_num_vec == changepoint_num),]),2,mean,na.rm =T)
apply(clin.RJMCMC[["lambda"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)


#Gibbs

res_final2<- readRDS(paste0(path_biometrics,"Glioblastoma Example\\Glio Gibbs2.rds"))

summary(res_final2[["stacked_df"]])

apply(res_final2[["stacked_df"]],2,mean,na.rm =T)
apply(res_final2[["stacked_df"]],2,sd,na.rm =T)
apply(res_final2[["stacked_df"]],2,quantile,na.rm =T, probs = c(0.025,0.5,0.975))

Lambda_column_mu <- c(apply(res_final2[["stacked_df"]],2,mean,na.rm =T)[3:5], #Gibbs Sampler
                      na.omit(apply(clin.RJMCMC[["lambda"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)),#Uncollapsed Collapsing model
                      na.omit(apply(clin.RJMCMC[["lambda.df_mu"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)) #Collapsing Post-hoc Expect
)
index_1 <- seq(from = 1, to =length(Lambda_column_mu), by = 3)
index_2 <- seq(from = 2, to =length(Lambda_column_mu), by = 3)
index_3 <- seq(from = 3, to =length(Lambda_column_mu), by = 3)

Lambda_column_mu1 <- Lambda_column_mu[index_1]
Lambda_column_mu2 <- Lambda_column_mu[index_2]
Lambda_column_mu3 <- Lambda_column_mu[index_3]

Lambda_column_sd <- c(apply(res_final2[["stacked_df"]],2,sd,na.rm =T)[3:5], #Gibbs Sampler
                      na.omit(apply(clin.RJMCMC[["lambda"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)),#Uncollapsed Collapsing model
                      na.omit(apply(sqrt(clin.RJMCMC[["lambda.df_var"]][which(chang_num_vec == changepoint_num),]),2,mean,na.rm =T)) #Collapsing Post-hoc Expect
)

Lambda_column_sd1 <- Lambda_column_sd[index_1]
Lambda_column_sd2 <- Lambda_column_sd[index_2]
Lambda_column_sd3 <- Lambda_column_sd[index_3]


Changepoint_column <- c(apply(res_final2[["stacked_df"]],2,mean,na.rm =T)[1:2],
                        na.omit(apply(clin.RJMCMC[["changepoint"]][which(chang_num_vec == changepoint_num),],2,mean,na.rm =T)),
                        NA,NA)


Changepoint_sd <- c(apply(res_final2[["stacked_df"]],2,sd,na.rm =T)[1:2],
                    na.omit(apply(clin.RJMCMC[["changepoint"]][which(chang_num_vec == changepoint_num),],2,sd,na.rm =T)),
                    NA,NA)


index_1 <- seq(from = 1, to =length(Changepoint_column), by = 2)
index_2 <- seq(from = 2, to =length(Changepoint_column), by = 2)

Changepoint_column_mu1 <- Changepoint_column[index_1]
Changepoint_column_mu2 <- Changepoint_column[index_2]

Changepoint_column_sd1 <- Changepoint_sd[index_1]
Changepoint_column_sd2 <- Changepoint_sd[index_2]

results.vec_glio <- cbind(Lambda_column_mu1, Lambda_column_sd1,
                          Lambda_column_mu2, Lambda_column_sd2,
                          Lambda_column_mu3,Lambda_column_sd3,
                          Changepoint_column_mu1,Changepoint_column_sd1,
                          Changepoint_column_mu2, Changepoint_column_sd2)

results.vec_glio <- apply(results.vec_glio,2, round,digits = 5)


rownames(results.vec_glio) <- c("Gibbs", "Collapsing - Uncollapsed", "Collapsing - Post Hoc")
colnames(results.vec_glio) <- c("Lambda 1", "sd", "Lambda 2", "sd","Lambda 3", "sd", "Changepoint 1", "sd","Changepoint 2", "sd")




res_final1 <- gibbs.changepoint_chains_optim_rcpp(clinfinal, 1, 10100, 100,1,alpha.hyper = 1, beta.hyper = 365, max_predict = 2000, interval = 50)
res_final2 <- gibbs.changepoint_chains_optim_rcpp(clinfinal, 2, 10100, 100,1,alpha.hyper = 1, beta.hyper = 365, max_predict = 2000, interval = 50)
res_final3 <- gibbs.changepoint_chains_optim_rcpp(clinfinal, 3, 10100, 100,1,alpha.hyper = 1, beta.hyper = 365,max_predict = 2000, interval = 50)
res_final4 <- gibbs.changepoint_chains_optim_rcpp(clinfinal, 4, 10100, 100,1,alpha.hyper = 1, beta.hyper = 365,max_predict = 2000, interval = 50)
res_final5 <- gibbs.changepoint_chains_optim_rcpp(clinfinal, 5, 10100, 100,1,alpha.hyper = 1, beta.hyper = 365,max_predict = 2000, interval = 50)
res_final6 <- gibbs.changepoint_chains_optim_rcpp(clinfinal, 6, 10100, 100,1,alpha.hyper = 1, beta.hyper = 365,max_predict = 2000, interval = 50)

saveRDS(res_final1, paste0(path_biometrics,"Glioblastoma Example\\Glio Gibbs1.rds"))
saveRDS(res_final2, paste0(path_biometrics,"Glioblastoma Example\\Glio Gibbs2.rds"))
saveRDS(res_final3, paste0(path_biometrics,"Glioblastoma Example\\Glio Gibbs3.rds"))
saveRDS(res_final4, paste0(path_biometrics,"Glioblastoma Example\\Glio Gibbs4.rds"))
saveRDS(res_final5, paste0(path_biometrics,"Glioblastoma Example\\Glio Gibbs5.rds"))
saveRDS(res_final6, paste0(path_biometrics,"Glioblastoma Example\\Glio Gibbs6.rds"))

marg_res_days <- c(as.brob(exp(margin_lik_calc_log(df_recast(clinfinal),NA, alpha = 1,beta =365))),
                   res_final1[["mean.marg.lik.log"]],
                   res_final2[["mean.marg.lik.log"]],
                   res_final3[["mean.marg.lik.log"]],
                   res_final4[["mean.marg.lik.log"]],
                   res_final5[["mean.marg.lik.log"]],
                   res_final6[["mean.marg.lik.log"]])

gibbs.model.selc(marg_res_days)


#Same parameters for large sample sizes between Collapsing and 

round(t(apply(res_final2[["stacked_df"]],2,quantile, probs = c(0.025,0.5,0.975))),5)
num.changepoints <- apply(clin.RJMCMC[["changepoint"]],1,function(x){length(na.omit(x))})
rbind(t(apply(clin.RJMCMC[["changepoint"]][which(num.changepoints==2),1:2],2,quantile, probs = c(0.025,0.5,0.975))),
round(t(apply(clin.RJMCMC[["lambda"]][which(num.changepoints==2),1:3],2,quantile, probs = c(0.025,0.5,0.975))),5))

rbind(t(apply(clin.RJMCMC[["changepoint"]][which(num.changepoints==2),1:2],2,quantile, probs = c(0.025,0.5,0.975))),
      round(t(apply(clin.RJMCMC[["lambda"]][which(num.changepoints==2),1:3]*365.25,2,quantile, probs = c(0.025,0.5,0.975))),5))



mean(clin.RJMCMC[["changepoint"]][which(num.changepoints==1),1])

round(t(apply(res_final1[["stacked_df"]],2,mean)),5)









source(paste0(path, "Functions Gibbs Markdown.R"))
source(paste0(path, "Functions Reversible Jump chains final.R"))



# Proof that it is invariatant to the lambda prior
clin.RJMCMC2 <- RJMC.piece(10000, clinfinal,n.chains = 1,  burn_in = 750,
                           max_predict = 15, max_plt = 15, lambda.prior = 2)

clin.RJMCMC3 <- RJMC.piece(10000, clinfinal,n.chains = 1,  burn_in = 750,
                           max_predict = 15, max_plt = 15, lambda.prior = 3)

clin.RJMCMC4 <- RJMC.piece(10000, clinfinal,n.chains = 1,  burn_in = 750,
                           max_predict = 15, max_plt = 15, lambda.prior = 4)


clin.RJMCMC4[["prob.changepoint"]]




clinfinal_years <- clinfinal[,c("time","status")]

clinfinal_years$time <- clinfinal$time/365.25


glio.gibbs2 <- gibbs.changepoint_chains_optim_rcpp(clinfinal_years, 2, 10100, 100,alpha.hyper =1, beta.hyper = 1, max_predict = 10)




#num.changepoints <-unlist(apply(clin.RJMCMC[["k.stacked"]],1, function(x){length(na.omit(x))}))
#colMeans(clin.RJMCMC[["changepoint"]][which(num.changepoints ==2),],na.rm= T)
#colMeans(clin.RJMCMC[["lambda"]][which(num.changepoints ==2),],na.rm= T)
#apply(clin.RJMCMC[["lambda"]][which(num.changepoints ==2),],2, quantile, probs = c(0.025,0.975), na.rm = T)
#apply(clin.RJMCMC[["changepoint"]][which(num.changepoints ==2),],2, quantile, probs = c(0.025,0.975), na.rm = T)

#clin.RJMCMC[["output"]][["2"]]

#install.packages("bshazard")
library(bshazard)
clinfinal_years$enter <- 0
bshazard.fit <- bshazard(formula = Surv(enter, time,
                                        status) ~ 1,data = clinfinal_years,nbin = 
                           length(unique(clinfinal_years[which(clinfinal_years$status == 1),"time"])))


plot(bshazard.fit)

bshazard.df  <- data.frame(time = bshazard.fit[["time"]],
                           hazard = bshazard.fit[["hazard"]],
                           upper.ci = bshazard.fit[["upper.ci"]],
                           lower.ci = bshazard.fit[["lower.ci"]])



harzard.plt.mod <- function(time.vec, cens.vec, Outcome = "Survival",  
                            lrg.haz.int = 1,  bw.method = "local"){
  
  #Plot the hazards
  max.time <- max(time.vec)
  result.hazard.pe.lrg <- pehaz(time.vec, cens.vec, width= lrg.haz.int, max.time=max.time)
  
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


hazard.PFS.mod <-harzard.plt.mod(time.vec = clinfinal_years$time , cens.vec =  clinfinal_years$status , 
                                 Outcome = "Death", lrg.haz.int = 1)

df.haz <- data.frame(time = hazard.PFS.mod$Cuts[-1]-1, 
                     hazard = c(hazard.PFS.mod$Hazard))


samp.plot<- glio.gibbs2[["samp.plot"]]
id_retain <- sample(max(samp.plot$id),500)
samp.plot<- samp.plot[which(samp.plot$id %in% id_retain),]

df.summary <- glio.gibbs2[["df.summary"]]

df.summary[nrow(df.summary),1] <- 9
path_latex <- "C:/Users/phili/OneDrive/PhD/Latex Folder/"
ggplot(bshazard.df[which(bshazard.df$time <= 10),], aes(x = time, y = hazard))+
  geom_line(colour = "blue", size = 1.3,linetype = "longdash")+
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci), alpha = 0.2)+
  geom_step(data = df.haz[which(df.haz$time <= 9),],  aes(time, hazard), size = 1.1,linetype = "dotdash")+
  scale_x_continuous(breaks=0:9)+
  xlab("Time")+
  ylab("Hazard")+
  ggtitle("Hazard Functions for Gioblastoma data")+
  geom_step(data = samp.plot,aes(timepoints, hazards,group = id), linetype = "dashed", alpha = 0.03, colour = "red")+
  geom_step( aes(timepoints, hazards,group = id),data = df.summary, colour = "purple", size = 1.3)+
  scale_y_continuous(breaks = seq(from = 0, to = 1.7, by =0.1))+
  coord_cartesian(ylim = c(0,1.7),xlim = c(0,8.5))
#Non-parametric Hazards
ggsave(paste0(path_latex,"Hazards Gioblastoma.png"), width = 7, height = 5)


results.gibbs <- glio.gibbs2 #gibbs.changepoint_chains_optim(clinfinal,num.breaks =  2,10100,100,1, max_predict = 9 )

haz.df <- glio.gibbs2[["stacked_df"]][,3:5]
apply(results.gibbs[["stacked_df"]],2, quantile,probs = c(0.025,.5,0.975))


haz.df2 <- data.frame(hazard = c(haz.df[,1],haz.df[,2],haz.df[,3]),
                      interval = rep(1:3, each =nrow(haz.df)))

haz.df2$interval <- as.factor(haz.df2$interval) 

hazard.post  <- ggplot(haz.df2, aes(x=interval, y=hazard, colour = interval)) + 
  geom_boxplot(outlier.shape = NA)+
  theme(legend.position = "none")+
  xlab("Interval")+
  ylab("Hazard")+
  ggtitle("Hazards")+
  theme(plot.title = element_text(size = 11))+
  scale_y_continuous(breaks = seq(from = 0, to = 1.2, by =0.1))+
  coord_cartesian(ylim = c(0,1.2))


gg.post <- ggarrange(results.gibbs[["plot_change"]]+
                       ggtitle("Changepoints")+
                       theme(legend.position = "top",plot.title = element_text(size = 11)),
                     hazard.post)
annotate_figure(gg.post,top = text_grob("Posterior Distribution of Parameters",size = 14))
ggsave(paste0(path_latex,"Hazards Gioblastoma Change.png"), width = 10, height = 5)


#Interesting Datasets

#Don't need to load these datasets just examples
clin <- survivalTCGA(
  ACC.clinical,
  BLCA.clinical,
  BRCA.clinical,
  CESC.clinical,
  CHOL.clinical,
  COAD.clinical,
  COADREAD.clinical,
  DLBC.clinical,
  ESCA.clinical,
  GBM.clinical,
  GBMLGG.clinical,
  HNSC.clinical,
  KICH.clinical,
  KIPAN.clinical,
  KIRC.clinical,
  KIRP.clinical,
  LAML.clinical,
  LGG.clinical,
  LIHC.clinical,
  LUAD.clinical,
  LUSC.clinical,
  MESO.clinical,
  OV.clinical,
  PAAD.clinical,
  PCPG.clinical,
  PRAD.clinical,
  READ.clinical,
  SARC.clinical,
  SKCM.clinical,
  STAD.clinical,
  STES.clinical,
  STES.clinical,
  TGCT.clinical,
  THCA.clinical,
  THYM.clinical,
  UCEC.clinical,
  UCS.clinical,
  UVM.clinical, extract.cols="admin.disease_code")

#https://portal.gdc.cancer.gov/projects/TCGA-BRCA
#://wiki.cancerimagingarchive.net/display/Public/TCGA-BRCA


clin2 <- survivalTCGA(KIRC.clinical,
                      GBM.clinical,
                      OV.clinical,extract.cols="admin.disease_code")


sfit <- survfit(Surv(times, patient.vital_status)~admin.disease_code, data=clin2)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE)
ggsurvplot(sfit, conf.int=TRUE, fun = "cumhaz")


