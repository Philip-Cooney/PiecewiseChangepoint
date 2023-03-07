# #Required packages
#
# list.of.packages <- c("muhaz", "hesim","ggplot2", "Rcpp", "RcppArmadillo", "coda","flexsurv","eha","Brobdingnag", "survminer", "bindrcpp", "progress","DirichletReg","sfsmisc")
#
# #Check to see if these are installed and install if not
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)
#
# #load the packages
# lapply(list.of.packages, require, character.only = TRUE)
#
#
# piecewise_loglik <- function(df, changepoint, method = "ML", lambda = NULL) {
#
#   surv.object <- with(df, Surv(enter, time, status))
#   split <- SurvSplit(surv.object, changepoint)
#   n.ivl <- length(changepoint) + 1
#
#   T <- split$Y$exit - split$Y$enter
#   d <- split$Y$event
#   ivl <- split$ivl
#
#   if(method == "ML"){
#
#     alpha <- matrix(0, nrow = 1, ncol = n.ivl)
#     Dtot <- sum(d)
#     res <- -Dtot
#     for (i in seq_len(n.ivl)) {
#       indx <- ivl == i
#       D <- sum(d[indx])
#       if (D > 0) {
#         sumT <- sum(T[indx])
#         alpha[i] <- D/sumT
#         res <- res + D * (log(D) - log(sumT))
#       }
#     }
#     #return(list(loglik = res, hazards = alpha))
#     return(res)
#   }else{
#
#     death.loglik <-0
#     cum_haz.loglik <- 0
#
#     for (i in seq_len(n.ivl)) {
#       indx <- ivl == i
#       death.loglik <-  sum(d[indx])*log(lambda[i]) + death.loglik
#       cum_haz.loglik <- -sum(lambda[i]*T[indx]) +cum_haz.loglik
#
#     }
#     return(death.loglik +cum_haz.loglik)
#   }
# }
#
# #Predictive Likelihood
#
# k_indexes <- c(2,4)
#
# pred_lik <- function(time_diffs, k_indexes, hyper.param = 1){
#   n <- nrow(time_diffs)
#   split_vector <- split(1:n,cut(1:n,c(0,k_indexes, Inf)))
#   expo_death_df <- exposure_death_alt(time_diffs, k_indexes)
#   output_prob <- rep(NA, n)
#
#   for(i in 1:nrow(expo_death_df)){
#     alpha <- expo_death_df[i,1] +hyper.param
#     beta  <-expo_death_df[i,2] +hyper.param
#     vec <- split_vector[[i]]
#     for(j in 1:length(vec)){
#       #output_prob[vec[j]] <- (alpha*(beta^alpha))/((beta+ time_diffs[vec[j]])^(alpha+1))
#       output_prob[vec[j]] <- exp(log(alpha)+alpha*log(beta)-(alpha+1)*log(beta+ time_diffs[vec[j]]))
#     }
#
#   }
#   return(output_prob)
# }
#
#
#
# #Time and death intervals
# exposure_death <- function(df, changepoint) {
#
#   surv.object <- with(df, Surv(enter, time, status))
#   split <- SurvSplit(surv.object, changepoint)
#   n.ivl <- length(changepoint) + 1
#
#   T <- split$Y$exit - split$Y$enter
#   d <- split$Y$event
#   ivl <- split$ivl
#
#   result <- matrix(NA, ncol = 2, nrow = n.ivl)
#   for (i in seq_len(n.ivl)) {
#     indx <- ivl == i
#     result[i,1] <- sum(d[indx])
#     result[i,2] <- sum(T[indx])
#   }
#   return(result)
# }
#
# hazard.plot <- function(time.vec, cens.vec, lrg.haz.int = 1,  bw.method = "local"){
#
#   #Plot the hazards
#   max.time <- max(time.vec)
#
#
#   result.hazard.pe.lrg <- pehaz(time.vec, cens.vec, width= lrg.haz.int, max.time=max.time)
#
#   result.smooth <- muhaz(time.vec, cens.vec, b.cor="left",
#                          max.time=max.time, bw.method = bw.method)
#
#   #Plot the Survival function
#
#   plot(result.hazard.pe.lrg,  col="black")
#   lines(result.smooth, col = "red", lty= 2)
#
#   #return(result.smooth)
# }
#
#
#
#
# df_hazard_plot <- function(df, time.vec, num.breaks){
#
#   changepoint_index <-grep("changepoint", colnames(df))
#   lambda_index <-grep("lambda", colnames(df))
#
#   change_vec <- t(df[,changepoint_index])
#   change_vec2 <- data.frame(rbind(0,change_vec, max(time.vec)))
#
#   lambda_vec <- t(df[,lambda_index])
#   lambda_vec2 <- data.frame(rbind(lambda_vec,lambda_vec[num.breaks+1,]))
#
#   df.plot <- data.frame(timepoints = unlist(change_vec2),
#                         hazards = unlist(data.frame(lambda_vec2)),
#                         id = rep(1:nrow(df), each = (num.breaks+2)))
#   return(df.plot)
# }
#
#
# round_df <- function(x, digits) {
#   # round all numeric variables
#   # x: data frame
#   # digits: number of digits to round
#   #https://stackoverflow.com/questions/29875914/rounding-values-in-a-dataframe-in-r/29876220
#   numeric_columns <- sapply(x, mode) == 'numeric'
#   x[numeric_columns] <-  round(x[numeric_columns], digits)
#   x
#
#
# }
#
#
# df_recast <- function(df){
#
#   #Takes the data and reformats it into a time between each event format
#
#   df_event  <-df[which(df$status == 1),c("status","time")]
#   df_event <-  df_event_origin <- df_event[order(df_event$time),]
#   n_event <- nrow(df_event)
#
#   time_diff_event <- time_diff_orgin <- c(df_event$time[1]*n_event,diff(df_event$time)*((n_event-1):1))
#
#   if(any(time_diff_event==0)){
#
#
#     time_diff_distinct_event <- time_diff_event[-which(time_diff_event==0)]
#
#     time_diff_event <- time_diff_distinct_event
#   }
#
#   n_distinct_events <- length(time_diff_event)
#
#   if(length(which(df$status == 0))>0){
#     df_event <- unique(df_event)
#     ind_time_diff <- c(df_event$time[1],diff(df_event$time))
#     df_cens <- df[which(df$status == 0),]
#     n_cens <- length(which(df$status == 0))
#     time_diff_cens <- matrix(NA, nrow = n_cens,ncol = length(time_diff_event))
#
#
#     for(i in 1:n_cens){
#       cens_obs <- df_cens$time[i]
#       if(cens_obs <=  df_event$time[1]){
#         time_diff_cens[i,1] <- cens_obs
#         next
#       }
#       for(j in 1:n_distinct_events){
#         diff <-  cens_obs - df_event$time[j]
#         if(diff <= 0){
#           time_diff_cens[i,j] <- cens_obs - df_event$time[j-1]
#           break
#         }else if(j==n_distinct_events && cens_obs > df_event$time[j]){
#           time_diff_cens[i,j] <- ind_time_diff[j] + cens_obs - df_event$time[j]
#         }else{
#           time_diff_cens[i,j] <- ind_time_diff[j]
#         }
#
#       }
#     }
#     time_diff_cens_sum <- colSums(time_diff_cens, na.rm = T)
#     time_diff_event <- time_diff_event +time_diff_cens_sum
#   }
#   return(cbind(time_diff_event,table(df_event_origin$time)))
# }
#
#
#
#
# piecewise_loglik.indiv <- function(df, changepoint, lambda = NULL) {
#
#   surv.object <- with(df, Surv(enter, time, status))
#   split <- SurvSplit(surv.object, changepoint)
#   n.ivl <- length(changepoint) + 1
#
#   T <- split$Y$exit - split$Y$enter
#   d <- split$Y$event
#   ivl <- split$ivl
#
#
#   n.obs <- length(unique(split$idx))
#   log.lik.df <- array(NA, dim = c(n.obs, 1))
#
#   for(id in 1:n.obs){
#     death.loglik <- 0
#     cum_haz.loglik <- 0
#     for (i in seq_len(n.ivl)) {
#       indx <- ivl == i
#       id_index <- (split$idx == id) & indx
#
#       death.loglik <-  sum(d[id_index])*log(lambda[i]) + death.loglik
#       cum_haz.loglik <- -sum(lambda[i]*T[id_index]) +cum_haz.loglik
#
#     }
#
#     log.lik.df[id,1] <- death.loglik +cum_haz.loglik
#
#   }
#   return(log.lik.df)
#
# }
#
#
# #Full Suite of functions for the Gibbs Sampler
#
# PML.calc <- function(origin.df, output_df, num.breaks, plot.conver = F, boot.samps = 10000){
#   n <- nrow(origin.df)
#
#   loglik.indiv.df <- matrix(NA, nrow = n, ncol = nrow(output_df))
#
#   for(i in 1:nrow(output_df)){
#
#     loglik.indiv.df[,i] <- piecewise_loglik.indiv(df = origin.df, changepoint = output_df[i,1:num.breaks],
#                                                   output_df[i,(num.breaks+1):(2*(num.breaks)+1)])
#   }
#
#   #PML calculation as per Jackson
#   lik.indiv.df <- 1/exp(loglik.indiv.df)
#   PML.vec <- nrow(output_df)/apply(lik.indiv.df, 1, sum)
#
#   if(plot.conver == T){
#     lik.indiv.mean <- t(apply(lik.indiv.df, 1, cumsum)/(1:nrow(output_df)))
#     largest.val.index <- which(lik.indiv.df == max(lik.indiv.df), arr.ind = TRUE)[1]
#     #plot the convergence of the PML for the last event (largest)
#     plot(y = lik.indiv.mean[largest.val.index,], x = 1:nrow(output_df))
#
#   }
#
#
#   print(paste0("Log Marginal Likelihood is :", log(prod(PML.vec))))
#
#   #return(PML.vec)
#
#   #Assuming the values are stable we take the mean
#
#   log.PML.vec <- log(PML.vec)
#   log.PML <- rep(NA,boot.samps)
#   alpha <- c(rep(1,n))
#
#   for(i in 1:boot.samps){
#     weights <- DirichletReg::rdirichlet(1, alpha)
#     log.PML[i] <- n*sum(weights*log.PML.vec)
#   }
#   return(log.PML)
# }
#





gibbs.sampler <- function(df,
                          num.breaks, # Number of Changepoints
                          n.iter,
                          burn_in = 100,
                          num.chains = 2,
                          n.thin = 1,
                          alpha.hyper = 1, #Hyperparameter for Marginal Likelihood
                          beta.hyper = 1, #Hyperparameter for Marginal Likelihood
                          MLE = FALSE #If TRUE, considers on Uniform Prior
                          ){

  df <- df[order(df$time),]
  df_event <- df[which(df$status == 1),]
  df_event_unique <- unique(df_event[,c("time","status")])
  # data processing and sampler setup
  time_diffs <- df_recast(df)
  n <- nrow(time_diffs)
  n.vec <- 1:n # Need to update in the presence of same timepoints
  m <- n.iter #length of the chain
  n.chains <- num.chains
  num.breaks <- num.breaks



  output_df <- array(NA, dim = c(m,2*(num.breaks) + 1,n.chains))
  #objects and names for output
  mcmc.output <- list()
  #Create names for the outputs
  changepoint_names <- rep(NA, num.breaks)
  lambda_names <- rep(NA, num.breaks+1)
  for(i in 1:num.breaks){
    changepoint_names[i] <- paste0("changepoint_", i)
  }
  for(i in 1:(num.breaks+1)){
    lambda_names[i] <- paste0("lambda_", i)
  }
  names_vector <- c(changepoint_names,lambda_names)

  k <- array(NA, dim = c(m, num.breaks, n.chains))
  lambda  <- array(NA, dim = c(m, num.breaks+1,n.chains))


  for(c in 1:n.chains){
    #print(paste0("Chain Number ", c))
    #pb <- progress_bar$new(total = m)
    #Create vectors to store the changepoints and lambdas
    k_temp <- array(NA, dim = c(m, num.breaks))

    #lambda_temp <- array(NA, dim = c(m, num.breaks+1))


    int_unorder <- sample(2:(n-1), num.breaks, replace = FALSE)
    int_locs <- int_unorder[order(int_unorder)]

    while(pos_prob(int_locs, n, MLE) == 0){

      int_unorder <- sample(2:(n-1), num.breaks, replace = FALSE)
      int_locs <- int_unorder[order(int_unorder)]
    }

    k_temp[1,] <- int_locs



    #Initial alpha and beta Hyperparameters
    rnd.draws <- rgamma(m*(num.breaks + 1), shape = alpha.hyper, rate = beta.hyper)
    rnd.draws[which(rnd.draws < 5.701102e-247)] <- 5.701102e-247
    beta_array <- alpha_array <- array(rnd.draws,dim = c(m, (num.breaks+1)))

    #run the Gibbs sampler
    chain_res <- Gibbs_core(k_temp,num.breaks,time_diffs,alpha_array,beta_array,MLE, alpha.hyper, beta.hyper)
    k[,,c] <- chain_res[["k"]]
    lambda[,,c] <- chain_res[["lambda"]]


    #output_df[,,c] <-  cbind(k,lambda)


  }

  #Discard the iterations
  samp_retain <- seq(from = n.thin, to = m, by = n.thin)

    k <- k[samp_retain[which(samp_retain > burn_in)], , , drop = F]
    lambda <- lambda[samp_retain[which(samp_retain > burn_in)], , , drop = F]



  # Stack the outputs from the chains
  if (n.chains == 1) {
    k.stacked <- k[,,1]
    lambda.stacked <- lambda[,,1]

  }
  if (n.chains > 1) {
      k.stacked <- as.matrix(k[, , 1])
      lambda.stacked <- as.matrix(lambda[, , 1])
    for (i in 2:n.chains) {
      k.stacked <- rbind(k.stacked, as.matrix(k[, , i]))
      lambda.stacked <- rbind(lambda.stacked, as.matrix(lambda[, , i]))

    }
  }

  changepoint <- matrix(nrow = nrow(k.stacked), ncol = ncol(k.stacked))


  for (i in 1:nrow(k.stacked)) {

          changepoint[i, 1:num.breaks] <- df_event_unique[k.stacked[i, 1:num.breaks], "time"]

      }



  # for(c in 1:n.chains){
  #
  #   if(n.chains == 1){
  #     output_df_slice <-  output_df
  #   }else{
  #     output_df_slice <-  output_df[,,c]
  #   }
  #
  #   colnames(output_df_slice) <- names_vector
  #   chain_id <- paste0("chain_",c)
  #
  #   if(num.breaks == 1){
  #     output_df_slice[,1] <- sapply(output_df_slice[,1], function(x)df_event_unique[x,"time"] )
  #   }else{
  #     output_df_slice[,1:num.breaks] <- apply(output_df_slice[,1:num.breaks],c(1,2), function(x)df_event_unique[x,"time"] )
  #   }
  #
  #   mcmc.output[[chain_id]] <- as.mcmc(output_df_slice)
  #
  # }
  #
  # mcmc.output <- as.mcmc.list(mcmc.output)
  # stacked_df <- do.call(rbind, mcmc.output)


  # if(num.breaks == 1){
  #   change.point_df <- data.frame(changetime = stacked_df[,1],
  #                                 changepoint = 1)
  # }else{
  #   change.point_df <- data.frame(changetime = c(unlist(stacked_df[,1:num.breaks])),
  #                                 changepoint = rep(1:num.breaks,each = nrow(stacked_df)))
  # }
  #
  # change.point_df$changepoint <- factor(change.point_df$changepoint)
  #
  # change.point_plot <- change.point_df %>% group_by(changepoint, changetime) %>%
  #   dplyr::summarize(n = dplyr::n()) %>% mutate(perc = (n*100/nrow(stacked_df)))

  #Do ths outside!

  # plot_change  <- ggplot(change.point_plot[which(change.point_plot$perc >0.5),], aes(x = changetime, y = perc, color=changepoint))+
  #   geom_pointrange(aes(ymin=0, ymax=perc), size = 0.02)+
  #   scale_x_continuous(name="Time", breaks = round(seq(round(min(change.point_plot$changetime),1),
  #                                                      max(change.point_plot[which(change.point_plot$perc >0.5),
  #                                                                            "changetime"]), by = interval),1) )+
  #   scale_y_continuous(name="Probability of Changepoint (%)")+
  #   ggtitle("Posterior distribution of changepoints")


  gibbs.obj <- list(k.stacked = k.stacked,
                   changepoint = changepoint,
                   lambda = lambda.stacked,
                   k = k,
                   lambda.array = lambda,
                   df = df)

  class(gibbs.obj) <- c("gibbs_changepoint")

    return(gibbs.obj)


}





