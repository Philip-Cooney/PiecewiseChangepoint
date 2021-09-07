

object <- Collapsing.Model

plot.changepoint <- function(object, type = "survival",chng.num = "all", add.km = T, max_predict = 10, interval = 0.2,...){
 if(type == "survival"){
   St <- get_Surv(object, chng.num = chng.num, max_predict = max_predict, interval = interval)
   return(plot.Survival(St, add.km = add.km))
 }

  if(type == "hazard"){
    return(plot.hazard(object, chng.num = chng.num))
  }


}

plot.Survival <- function(St,max.num.post = 500, add.km = T, env = parent.frame()){

  nSims <- ncol(St)
  time <- as.numeric(rownames(St))

  mean.Surv_df <- data.frame(survival = apply(St,1,FUN = mean),
                             time = time)
  mod_quantile_df <- data.frame(cbind(time,t(apply(St, 1,FUN = quantile, c(0.025,0.975)))))


  colnames(mod_quantile_df) <- c("time", "lower", "upper")


  Surv.plot <- data.frame(Survival = c(unlist(St)),
                          time = rep(time,nSims),
                          id = rep(1:nSims, each = length(time)))

  if(max.num.post < nSims){
    post_id <-  sample(1:nSims, size = max.num.post)
    Surv.plot <-  dplyr::filter(Surv.plot,id %in% post_id)
  }


  plot_Surv <- ggplot(data = Surv.plot, mapping = aes(x = time, y = Survival, group = id))+
    geom_line(size = 0.1, alpha = 0.05, colour = "red")+
    geom_line(data = mean.Surv_df, aes(x = time, y = survival),size= 1, inherit.aes = F, colour = "purple")+
    geom_line(data = mod_quantile_df, aes(x = time, y = lower),linetype = "dashed",size= 1, inherit.aes = F, colour = "grey")+
    geom_line(data = mod_quantile_df, aes(x = time, y = upper),linetype = "dashed",size= 1, inherit.aes = F, colour = "grey")+
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),expand = c(0, 0))+
    scale_x_continuous(expand = c(0,0))+
    annotate(geom="segment", x=seq(0,max(time),1), xend = seq(0,max(time),1),
             y=0, yend= 0.01)+ theme_classic()

  if(env$chng.num != "all"){

    k <- env$object$k.stacked
    num.changepoints <-unlist(apply(k,1, function(x){length(na.omit(x))}))
    k_curr <- data.frame(k[which(num.changepoints == env$chng.num),1:env$chng.num])
    df <- env$object$df
    df_event  <-df[which(df$status == 1),c("status","time")]

    time.break <- df_event[apply(k_curr, 2,FUN = mean),"time"]
    survival.close <- sapply(time.break,FUN = function(x){which.min(abs(mean.Surv_df$time - x))})

    break.points.Surv <- data.frame(time = time.break,
                                    Survival =  mean.Surv_df$survival[survival.close])

    plot_Surv <-  plot_Surv+
      geom_point(data = break.points.Surv, aes(x = time, Survival),
                 shape=23, fill = "green",inherit.aes = F)

  }

  if(add.km){
    plot_Surv  <- add_km(plot_Surv,df)
  }

  return(plot_Surv)
}

plot.hazard <- function(object, chng.num = "all",max.num.post = 500){

  k <- object$k.stacked
  changepoint <- object$changepoint
  lambda <- object$lambda
  df <- object$df

  time.seq <- c(seq(from = 0, to  = max(df$time), by = max(df$time)/100))

  num.changepoints <-unlist(apply(k,1, function(x){length(na.omit(x))}))

  if(chng.num != "all"){
    lambda <- as.matrix(lambda[which(num.changepoints == chng.num),1:(chng.num+1)])
    changepoint <- as.matrix(changepoint[which(num.changepoints == chng.num),1:chng.num])
    num.changepoints <- num.changepoints[which(num.changepoints == chng.num)]
  }

  lambda_res_final <- NULL

  for(i in seq_along(unique(num.changepoints))){

    index <- unique(num.changepoints)[order(unique(num.changepoints))][i]

    # if(length(which(num.changepoints == index))<2){
    #   next
    # }

    if(index == 0){

      lambda_curr <- lambda[which(num.changepoints == index),1:(index+1)]
      lambda_res_final <- matrix(rep(lambda_curr, each = length(time.seq)),
                                 nrow = length(lambda_curr),
                                 ncol =  length(time.seq),
                                 byrow = T)

      df.changepoint <- data.frame(timepoints = rep(c(0, max(df$time)),
                                                    times = length(lambda_curr)),
                                    hazards = rep(lambda_curr, each = 2),
                                    id = rep(1:length(lambda_curr), each = 2))


    }else{

      changepoint_curr <- as.matrix(changepoint[which(num.changepoints == index),1:index])
      lambda_curr <- lambda[which(num.changepoints == index),1:(index+1)]

      lambda_res_curr <- matrix(nrow = nrow(changepoint_curr),ncol = length(time.seq) )
      changepoint_curr_samp <- cbind(changepoint_curr, Inf)

        for(j in 1:length(time.seq)){
          index.lambda<- apply(changepoint_curr_samp, 1, function(x){which.min(time.seq[j]>x)})
          lambda_res_curr[,j] <- lambda_curr[cbind(1:nrow(changepoint_curr_samp),index.lambda)]
        }


        #lambda.list[[i]] <- lambda_res_curr
        lambda_res_final <- rbind(lambda_res_final,lambda_res_curr)

    }

  }

  df.hazard<- data.frame(timepoints = rep(time.seq,by = nrow(lambda_res_final)),
                                               hazards = c(t(lambda_res_final)),
                                               id = rep(1:nrow(lambda_res_final), each = length(time.seq)))


  lambda_res_final_df <- data.frame(t(apply(lambda_res_final,2,  FUN = quantile, probs = c(0.05, 0.5, 0.95))), time = time.seq)

  if(max.num.post < nrow(changepoint)){
    post_id <-  sample(1:nrow(changepoint), size = max.num.post)
    df.plot.final<- dplyr::filter(df.hazard,id %in% post_id)

  }

  plot_haz <- ggplot(df.plot.final, aes(timepoints, hazards))+
    geom_step(aes(group = id), linetype = "dashed", alpha = 0.05, colour = "red")+
    geom_line(data = lambda_res_final_df, aes(time, X50.), size = 1.5)+
    geom_line(data = lambda_res_final_df, aes(time, X5.),linetype = "dotted", size = 1.25)+
    geom_line(data = lambda_res_final_df, aes(time, X95.), linetype = "dotted",size = 1.25)+
    scale_x_continuous(breaks = seq(0, max(df$time), by = interval), expand = c(0,0))+
    scale_y_continuous(expand = c(0,0), limits = c(0, NA))+
    theme_classic()

    return(plot_haz)


  }

ind.expo <- function(time,status, lambda){
  if(status ==0){
    return(pexp(time,rate = lambda, lower.tail = F, log.p = TRUE ))
  }else{
    return(dexp(time,rate = lambda, log = TRUE ))
  }
}

piecewise_loglik.indiv <- function(df, changepoint, lambda = NULL){

  df$enter <- 0
  surv.object <- with(df, Surv(enter, time, status))
  split <- eha::SurvSplit(surv.object, changepoint)
  n.ivl <- length(changepoint) + 1

  T <- split$Y$exit - split$Y$enter
  d <- split$Y$event
  ivl <- split$ivl


  n.obs <- length(unique(split$idx))
  log.lik.df <- array(NA, dim = c(n.obs, 1))

  for(id in 1:n.obs){
    death.loglik <- 0
    cum_haz.loglik <- 0
    for (i in seq_len(n.ivl)) {
      indx <- ivl == i
      id_index <- (split$idx == id) & indx

      death.loglik <-  sum(d[id_index])*log(lambda[i]) + death.loglik
      cum_haz.loglik <- -sum(lambda[i]*T[id_index]) +cum_haz.loglik

    }

    log.lik.df[id,1] <- death.loglik +cum_haz.loglik

  }
  return(log.lik.df)

}

get.loglik <- function(object, chng.num = "all"){

  k <- object$k.stacked
  changepoint <- object$changepoint
  lambda <- object$lambda
  df <- object$df


  num.changepoints <-unlist(apply(k,1, function(x){length(na.omit(x))}))

  if(chng.num != "all"){
    lambda <- as.matrix(lambda[which(num.changepoints == chng.num),1:(chng.num+1)])
    changepoint <- as.matrix(changepoint[which(num.changepoints == chng.num),1:chng.num])
    num.changepoints <- num.changepoints[which(num.changepoints == chng.num)]
  }

  lambda_res_final <- NULL

  indiv.log.lik.final <- NULL

  for(i in seq_along(unique(num.changepoints))){

    index <- unique(num.changepoints)[order(unique(num.changepoints))][i]

    # if(length(which(num.changepoints == index))<2){
    #   next
    # }

    if(index == 0){

      lambda_curr <- lambda[which(num.changepoints == index),1:(index+1)]

      indiv.log.lik.final <- matrix(NA, nrow = length(lambda_curr), ncol = nrow(df))

      for(q in 1:nrow(df)){
        indiv.log.lik.final[,q] <- sapply(lambda_curr, FUN = ind.expo, time = df$time[q],status = df$status[q] )
      }


    }else{

     lambda_curr <- lambda[which(num.changepoints == index),1:(index+1)]
     changepoint_curr <- changepoint[which(num.changepoints == index),1:index]
     indiv.log.lik <- matrix(NA, nrow = nrow(lambda_curr), ncol = nrow(df))

      for(x in 1:nrow(lambda_curr)){
        indiv.log.lik[x,] <- piecewise_loglik.indiv(df, as.numeric(data.frame(changepoint_curr)[x,]), lambda_curr[x,])

      }

      indiv.log.lik.final <- rbind(indiv.log.lik.final,indiv.log.lik)


    }

  }

  indiv.log.lik.final
}

add_km <- function(plt, df){

  result.km <- survfit(Surv(time,status)~1, data = df)
  km.data <- data.frame(cbind(result.km[[c("time")]],result.km[[c("surv")]], result.km[[c("upper")]],result.km[[c("lower")]]))
  colnames(km.data) <- c("time", "survival", "upper", "lower")

  plt+
    geom_step(data = km.data, aes(x = time, y = survival), inherit.aes = F )+
    geom_step(data = km.data, aes(x = time, y = upper),linetype = "dashed", inherit.aes = F )+
    geom_step(data = km.data, aes(x = time, y = lower), linetype = "dashed", inherit.aes = F )


}

get_Surv <- function(object,chng.num = "all", max_predict, interval ){

  time <- c(seq(from = 0, to  = max_predict, by = interval))
  k <- object$k.stacked
  lambda <- object$lambda
  changepoint <- object$changepoint
  num.changepoints <-unlist(apply(k,1, function(x){length(na.omit(x))}))

  if(chng.num != "all"){
    if(length(which(num.changepoints == chng.num)) < 2){
      print("Too few simulations for this change-point model")
      break
    }else{
      changepoints_eval <- chng.num

    }
  }else{
    changepoints_eval <- unique(num.changepoints)[order(unique(num.changepoints))]
  }

  St <- NULL
  for(i in 1:length(changepoints_eval)){

    if(changepoints_eval[i] == 0){

      num_zero <- length(which(num.changepoints == changepoints_eval[i]))
      lambda_0 <- data.frame(lambda[which(num.changepoints == changepoints_eval[i]),1])
      St <- cbind(St,surv_nochange(time,  num_zero, lambda_0[,1]))

    }else{
      k_curr <- data.frame(k[which(num.changepoints == changepoints_eval[i]),1:changepoints_eval[i]])
      changepoint_curr <- data.frame(changepoint[which(num.changepoints == changepoints_eval[i]),1:changepoints_eval[i]])
      lambda_curr <- lambda[which(num.changepoints == changepoints_eval[i]),1:(changepoints_eval[i]+1)]
      changepoint_df <- cbind(0,changepoint_curr[,1:changepoints_eval[i]])
      time.interval_df <- t(apply(changepoint_df,1,diff))

      if(changepoints_eval[i] == 1){
        cum_haz_df <-  time.interval_df*lambda_curr[,1]
        surv_df <-  cbind(1,t(exp(-cum_haz_df)))

      }else{
        #head(index2, -1) last hazard not needed for cumhaz calc
        cum_haz_df <- t(apply(time.interval_df*lambda_curr[,1:changepoints_eval[i]], 1,cumsum))
        surv_df <- cbind(1,exp(-cum_haz_df))

      }
      St <- cbind(St,surv_change(time,nrow(k_curr), lambda_curr, data.matrix(changepoint_df),surv_df))
    }
  }
  rownames(St) <- time

  St

}


surv.gg <- plot(Collapsing.Model, type = "survival", chng.num = 2, add.km = F, max_predict = 2,interval = 0.02)
add_km(surv.gg, df)


Collapsing.Model$changepoint

object <- Collapsing.Model

print.changepoint <- function(object, mod = NULL,  digits = min(3L, getOption("digits"))){
  cat("Posterior Change-point Probabilities:\n")

  names(attr(object$prob.changepoint,"dimnames")) <- NULL

  print.default(format(object$prob.changepoint, digits = digits), print.gap = 2L,
                quote = FALSE)

  if(is.null(mod)){
    chng.prob <- as.numeric(names(which.max(object$prob.changepoint)))
  }else{
    chng.prob <- mod
  }


  k <- object$k.stacked
  lambda <- object$lambda
  changepoint <- object$changepoint
  num.changepoints <-unlist(apply(k,1, function(x){length(na.omit(x))}))
  changepoint_curr <- data.frame(changepoint[which(num.changepoints == chng.prob),1:chng.prob])
  lambda_curr <- lambda[which(num.changepoints == chng.prob),1:(chng.prob+1)]

  lambda_summary <- summary(lambda_curr)
  colnames(lambda_summary) <- paste0(rep("lambda_",chng.prob+1), 1:(chng.prob+1))

  changepoint_summary <- summary(changepoint_curr)
  colnames(changepoint_summary) <- paste0(rep("changepoint_",chng.prob), 1:(chng.prob))

  cat(paste0("\n"))
  cat(paste0("Summary of ",chng.prob," changepoint model:\n"))
  cat(paste0("\n"))
  print.default(format(cbind(changepoint_summary,lambda_summary), digits = digits), print.gap = 2L,
                quote = FALSE)


}


#stats:::print.lm

print(Collapsing.Model, mod = 2)
