df_recast <- function(df) {

  # Takes the data and reformats it into a time between each event format

  df_event <- df[which(df$status == 1), c("status", "time")]
  df_event <- df_event_origin <- df_event[order(df_event$time), ]
  n_event <- nrow(df_event)

  time_diff_event <- time_diff_orgin <- c(df_event$time[1] * n_event, diff(df_event$time) * ((n_event - 1):1))

  if (any(time_diff_event == 0)) {
    time_diff_distinct_event <- time_diff_event[-which(time_diff_event == 0)]

    time_diff_event <- time_diff_distinct_event
  }

  n_distinct_events <- length(time_diff_event)

  if (length(which(df$status == 0)) > 0) {
    df_event <- unique(df_event)
    ind_time_diff <- c(df_event$time[1], diff(df_event$time))
    df_cens <- df[which(df$status == 0), ]
    n_cens <- length(which(df$status == 0))
    time_diff_cens <- matrix(NA, nrow = n_cens, ncol = length(time_diff_event))


    for (i in 1:n_cens) {
      cens_obs <- df_cens$time[i]
      if (cens_obs <= df_event$time[1]) {
        time_diff_cens[i, 1] <- cens_obs
        next
      }
      for (j in 1:n_distinct_events) {
        diff <- cens_obs - df_event$time[j]
        if (diff <= 0) {
          time_diff_cens[i, j] <- cens_obs - df_event$time[j - 1]
          break
        } else if (j == n_distinct_events && cens_obs > df_event$time[j]) {
          time_diff_cens[i, j] <- ind_time_diff[j] + cens_obs - df_event$time[j]
        } else {
          time_diff_cens[i, j] <- ind_time_diff[j]
        }
      }
    }
    time_diff_cens_sum <- colSums(time_diff_cens, na.rm = T)
    time_diff_event <- time_diff_event + time_diff_cens_sum
  }
  return(cbind(time_diff_event, table(df_event_origin$time)))
}


#' Chain mixing plot for piecewise exponential model
#'
#' @param object an object of class "changepoint".
#'
#' @return plot with the number of change-points at each simulation, coloured by chain number
#' @export
#'
#' @examples \dontrun{chain.mixing(Collapsing_Model)}

chain.mixing <- function(object) {
  k <- object$k
  k.stacked <- object$k.stacked
  chang.vals <- unique(apply(k.stacked, 1, function(x) {
    length(na.omit(x))
  }))
  if (is.na(dim(k)[3])) {
    n.chains <- 1
  } else {
    n.chains <- dim(k)[3]
  }

  k.first <- k[, , 1]
  plot(jitter(apply(k.first, 1, function(x) {
    length(na.omit(x))
  }), factor = 0.3), cex = 0.1, ylab = "Number of Changepoints", yaxt = "n")
  axis(side = 2, at = chang.vals[order(chang.vals)])
  if (n.chains > 1) {
    for (i in 2:n.chains) {
      points(jitter(apply(k[, , i], 1, function(x) {
        length(na.omit(x))
      }), factor = 0.3), cex = 0.1, col = i)
    }
  }
  legend("topright",
    inset = .02, title = "Chains",
    as.character(1:n.chains), fill = 1:n.chains, horiz = TRUE, cex = 0.8
  )
}




#' Piecewise Exponential Simulations
#' Generate piecewise exponential observations.
#' @param n_obs number of observations
#' @param n_events_req (approximate) number of events contained within the sample
#' @param num.breaks number of change-points, NA if no changepoints
#' @param rate vector of hazard rate
#' @param t_change vector of change-point times
#' @param max_time maximum time after which all observations are censored.
#'
#' @return a dataframe with three columns, the time of the events (time_event), time which is the minimum of the cenorsing time and event time. status is an indicator which is 0 censored observation, or 1 if event.
#' @export
#' @importFrom hesim rpwexp
#' @importFrom dplyr mutate
#' @examples
#' set.seed(123)
#' n_obs =300
#' n_events_req=300
#' max_time =  2
#' rate = c(0.75,0.25)
#' t_change =1
#' df <- gen_piece_df(n_obs = n_obs,n_events_req = n_events_req,
#'                   num.breaks = length(t_change),rate = rate ,
#'                   t_change = t_change, max_time = max_time)

gen_piece_df <- function(n_obs, n_events_req, num.breaks, rate, t_change, max_time = Inf) {
  n_cens_req <- n_obs - n_events_req
  ratemat <- matrix(rep(rate, n_obs / 2),
    nrow = n_obs,
    ncol = num.breaks + 1, byrow = TRUE
  )

  if (n_cens_req > 0) {
    if (num.breaks == 0) {
      samp_cens <- rexp(n_cens_req * 2, rate)
      samp <- rexp(n_obs, rate)
    } else {
      samp_cens <- rpwexp(n_cens_req * 2, ratemat, c(0, t_change)) # Assume that on average half the observations will be censors
      samp <- rpwexp(n_obs, ratemat, c(0, t_change))
    }
    samp_cens <- sapply(samp_cens, FUN = min, max_time)
    samp_cens <- sample(c(samp_cens, rep(max_time, n_obs - n_cens_req * 2))) # Randomized vector
    # http://www.cookbook-r.com/Manipulating_data/Randomizing_order/

    df <- data.frame(time_event = samp, time_cens = samp_cens)
    df$time <- apply(cbind(samp, samp_cens), 1, min)
    df <- df %>% mutate(status = ifelse(samp <= samp_cens, 1, 0), enter = 0)
  } else {
    if (num.breaks == 0) {
      samp <- rexp(n_obs, rate)
    } else {
      samp <- rpwexp(n_obs, ratemat, c(0, t_change))
    }
    df <- data.frame(time_event = samp)
    df <- df %>% mutate(
      status = ifelse(time_event > max_time, 0, 1),
      time = ifelse(time_event > max_time, max_time, time_event)
    )
  }
  df <- df[order(df$time), ]
  return(df)
}





posterior.changepoint <- function(changepoint.res, num.breaks) {
  num.changepoints <- unlist(apply(changepoint.res, 1, function(x) {
    length(na.omit(x))
  }))
  changepoint.res2 <- changepoint.res[which(num.changepoints == num.breaks), 1:num.breaks]

  ### Change plot function
  if (num.breaks == 1) {
    change.point_df <- data.frame(
      changetime = changepoint.res2,
      changepoint = 1
    )
  } else {
    change.point_df <- data.frame(
      changetime = c(unlist(changepoint.res2)),
      changepoint = rep(1:num.breaks, each = nrow(changepoint.res2))
    )
  }

  change.point_df$changepoint <- factor(change.point_df$changepoint)

  change.point_plot <- change.point_df %>%
    group_by(changepoint, changetime) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    mutate(perc = (n * 100 / nrow(changepoint.res2)))

  plot_change <- ggplot(change.point_plot, aes(x = changetime, y = perc, color = changepoint)) +
    geom_pointrange(aes(ymin = 0, ymax = perc), size = 0.02) +
    scale_x_continuous(name = "Time", breaks = round(seq(round(min(change.point_plot$changetime), 1),
      max(change.point_plot[
        which(change.point_plot$perc > 0),
        "changetime"
      ]),
      by = 0.5
    ), 1)) +
    scale_y_continuous(name = "Probability of Changepoint (%)") +
    ggtitle("Posterior Distribution of changepoints") +
    theme_bw()

  return(plot_change)
}


Mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

index.loc <- function(index, k.slice) {
  res <- max(which(index > k.slice))
  if (res == -Inf) {
    res <- 1
  } else {
    res <- res + 1
  }
  return(res)
}



#' Plot functions for change-point models
#'
#' @param object an object of class "changepoint".
#' @param type the type of plot to be drawn, default is the survival function. Also can plot the hazard function with "hazard".
#' @param chng.num value indicating the changepoint model to plotted, default is "all", in which all posterior simulations will be used.
#' @param add.km indicator to add Kaplan Meier curve. Only applicable to survival function. Default is true.
#' @param max_predict maximum time to be plotted. Default is 10, however, depending on the timescale this should be changed.
#' @param add.post indicator whether to plot the posterior simulations (a random sample of 500) for the survival function. Default is true.
#'
#' @return a ggplot object
#' @export
#' @examples \dontrun{
#' plot(Collapsing_Model, add.post = F)
#' plot(Collapsing_Model, type = "hazard")}

plot.changepoint <- function(object, type = "survival", chng.num = "all", add.km = T, max_predict = 10,  add.post = T, ...) {
  if (type == "survival") {
    St <- get_Surv(object, chng.num = chng.num, max_predict = max_predict)
    return(plot.Survival(St, add.km = add.km, add.post = add.post))
  }

  if (type == "hazard") {
    return(plot.hazard(object, chng.num = chng.num))
  }
}

plot.Survival <- function(St, max.num.post = 500, add.km, add.post, env = parent.frame()) {
  nSims <- ncol(St)
  time <- as.numeric(rownames(St))

  mean.Surv_df <- data.frame(
    survival = apply(St, 1, FUN = mean),
    time = time
  )
  mod_quantile_df <- data.frame(cbind(time, t(apply(St, 1, FUN = quantile, c(0.025, 0.975)))))


  colnames(mod_quantile_df) <- c("time", "lower", "upper")


  Surv.plot <- data.frame(
    Survival = c(unlist(St)),
    time = rep(time, nSims),
    id = rep(1:nSims, each = length(time))
  )

  if (max.num.post < nSims) {
    post_id <- sample(1:nSims, size = max.num.post)
    Surv.plot <- dplyr::filter(Surv.plot, id %in% post_id)
  }


  plot_Surv <- ggplot(data = Surv.plot, mapping = aes(x = time, y = Survival, group = id)) +
    geom_line(data = mean.Surv_df, aes(x = time, y = survival), size = 1, inherit.aes = F, colour = "purple") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    annotate(
      geom = "segment", x = seq(0, max(time), max(time)/50), xend = seq(0, max(time), max(time)/50),
      y = 0, yend = 0.01
    ) +
    theme_classic()

  if (add.post == T) {
    plot_Surv <- plot_Surv +
      geom_line(data = mod_quantile_df, aes(x = time, y = lower), linetype = "dashed", size = 1, inherit.aes = F, colour = "grey") +
      geom_line(data = mod_quantile_df, aes(x = time, y = upper), linetype = "dashed", size = 1, inherit.aes = F, colour = "grey") +
      geom_line(size = 0.1, alpha = 0.025, colour = "red")
  }



  if (env$chng.num != "all" && env$chng.num != 0) {
    k <- env$object$k.stacked
    num.changepoints <- unlist(apply(k, 1, function(x) {
      length(na.omit(x))
    }))
    k_curr <- data.frame(k[which(num.changepoints == env$chng.num), 1:env$chng.num])
    df <- env$object$df
    df_event <- df[which(df$status == 1), c("status", "time")]

    time.break <- df_event[apply(k_curr, 2, FUN = mean), "time"]
    survival.close <- sapply(time.break, FUN = function(x) {
      which.min(abs(mean.Surv_df$time - x))
    })

    break.points.Surv <- data.frame(
      time = time.break,
      Survival = mean.Surv_df$survival[survival.close]
    )

    plot_Surv <- plot_Surv +
      geom_point(
        data = break.points.Surv, aes(x = time, Survival),
        shape = 23, fill = "green", inherit.aes = F
      )
  }

  if (add.km) {
    plot_Surv <- add_km(plot_Surv, env$object$df)
  }

  return(plot_Surv)
}

plot.hazard <- function(object, chng.num = "all", max.num.post = 500, ...) {
  k <- object$k.stacked
  changepoint <- object$changepoint
  lambda <- object$lambda
  df <- object$df
  interval <- max(df$time) / 100

  time.seq <- c(seq(from = 0, to = max(df$time), by = interval))

  num.changepoints <- unlist(apply(k, 1, function(x) {
    length(na.omit(x))
  }))

  if (chng.num != "all") {
    lambda <- as.matrix(lambda[which(num.changepoints == chng.num), 1:(chng.num + 1)])
    changepoint <- as.matrix(changepoint[which(num.changepoints == chng.num), 1:chng.num])
    num.changepoints <- num.changepoints[which(num.changepoints == chng.num)]
  }

  lambda_res_final <- NULL

  for (i in seq_along(unique(num.changepoints))) {
    index <- unique(num.changepoints)[order(unique(num.changepoints))][i]

    # if(length(which(num.changepoints == index))<2){
    #   next
    # }

    if (index == 0) {
      lambda_curr <- lambda[which(num.changepoints == index), 1:(index + 1)]
      lambda_res_final <- matrix(rep(lambda_curr, each = length(time.seq)),
        nrow = length(lambda_curr),
        ncol = length(time.seq),
        byrow = T
      )

      df.changepoint <- data.frame(
        timepoints = rep(c(0, max(df$time)),
          times = length(lambda_curr)
        ),
        hazards = rep(lambda_curr, each = 2),
        id = rep(1:length(lambda_curr), each = 2)
      )
    } else {
      changepoint_curr <- as.matrix(changepoint[which(num.changepoints == index), 1:index])
      lambda_curr <- lambda[which(num.changepoints == index), 1:(index + 1)]

      lambda_res_curr <- matrix(nrow = nrow(changepoint_curr), ncol = length(time.seq))
      changepoint_curr_samp <- cbind(changepoint_curr, Inf)

      for (j in 1:length(time.seq)) {
        index.lambda <- apply(changepoint_curr_samp, 1, function(x) {
          which.min(time.seq[j] > x)
        })
        lambda_res_curr[, j] <- lambda_curr[cbind(1:nrow(changepoint_curr_samp), index.lambda)]
      }


      # lambda.list[[i]] <- lambda_res_curr
      lambda_res_final <- rbind(lambda_res_final, lambda_res_curr)
    }
  }

  df.hazard <- data.frame(
    timepoints = rep(time.seq, by = nrow(lambda_res_final)),
    hazards = c(t(lambda_res_final)),
    id = rep(1:nrow(lambda_res_final), each = length(time.seq))
  )


  lambda_res_final_df <- data.frame(t(apply(lambda_res_final, 2, FUN = quantile, probs = c(0.05, 0.5, 0.95))), time = time.seq)

  if (max.num.post < nrow(changepoint)) {
    post_id <- sample(1:nrow(changepoint), size = max.num.post)
    df.plot.final <- dplyr::filter(df.hazard, id %in% post_id)
  }

  plot_haz <- ggplot(df.plot.final, aes(timepoints, hazards)) +
    geom_step(aes(group = id), linetype = "dashed", alpha = 0.05, colour = "red") +
    geom_line(data = lambda_res_final_df, aes(time, X50.), size = 1.5) +
    geom_line(data = lambda_res_final_df, aes(time, X5.), linetype = "dotted", size = 1.25) +
    geom_line(data = lambda_res_final_df, aes(time, X95.), linetype = "dotted", size = 1.25) +
    scale_x_continuous(breaks = seq(0, max(df$time), length.out = 11), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    theme_classic()

  return(plot_haz)
}

ind.expo <- function(time, status, lambda) {
  if (status == 0) {
    return(pexp(time, rate = lambda, lower.tail = F, log.p = TRUE))
  } else {
    return(dexp(time, rate = lambda, log = TRUE))
  }
}

piecewise_loglik.indiv <- function(df, changepoint, lambda = NULL) {
  df$enter <- 0
  surv.object <- with(df, Surv(enter, time, status))
  split <- eha::SurvSplit(surv.object, changepoint)
  n.ivl <- length(changepoint) + 1

  T <- split$Y$exit - split$Y$enter
  d <- split$Y$event
  ivl <- split$ivl


  n.obs <- length(unique(split$idx))
  log.lik.df <- array(NA, dim = c(n.obs, 1))

  for (id in 1:n.obs) {
    death.loglik <- 0
    cum_haz.loglik <- 0
    for (i in seq_len(n.ivl)) {
      indx <- ivl == i
      id_index <- (split$idx == id) & indx

      death.loglik <- sum(d[id_index]) * log(lambda[i]) + death.loglik
      cum_haz.loglik <- -sum(lambda[i] * T[id_index]) + cum_haz.loglik
    }

    log.lik.df[id, 1] <- death.loglik + cum_haz.loglik
  }
  return(log.lik.df)
}

#' Pointwise log-likelihood for the change-point model
#'
#' Provides the pointwise log-likelihood contribution for each datapoint and simulation which is required to compute the Pseudo-Marginal Likelihood (PML) and Widely Applicable Information Criterion (WAIC).
#'
#' @param object of class "changepoint".
#' @param chng.numvalue  indicating the changepoint model to plotted, default is "all", in which all posterior simulations will be used.
#'
#' @return a matrix of size nSims by n_obs (number of simulations and number of observations respectively).
#' @export
#'
#' @examples \dontrun{
#' Can take a considerable time to evaulate.
#' log.lik.df <- get.loglik(Collapsing_Model)
#' }
get.loglik <- function(object, chng.num = "all") {
  k <- object$k.stacked
  changepoint <- object$changepoint
  lambda <- object$lambda
  df <- object$df


  num.changepoints <- unlist(apply(k, 1, function(x) {
    length(na.omit(x))
  }))

  if (chng.num != "all") {
    lambda <- as.matrix(lambda[which(num.changepoints == chng.num), 1:(chng.num + 1)])
    changepoint <- as.matrix(changepoint[which(num.changepoints == chng.num), 1:chng.num])
    num.changepoints <- num.changepoints[which(num.changepoints == chng.num)]
  }

  lambda_res_final <- NULL

  indiv.log.lik.final <- NULL

  for (i in seq_along(unique(num.changepoints))) {
    index <- unique(num.changepoints)[order(unique(num.changepoints))][i]

    # if(length(which(num.changepoints == index))<2){
    #   next
    # }

    if (index == 0) {
      lambda_curr <- lambda[which(num.changepoints == index), 1:(index + 1)]

      indiv.log.lik.final <- matrix(NA, nrow = length(lambda_curr), ncol = nrow(df))

      for (q in 1:nrow(df)) {
        indiv.log.lik.final[, q] <- sapply(lambda_curr, FUN = ind.expo, time = df$time[q], status = df$status[q])
      }
    } else {
      lambda_curr <- lambda[which(num.changepoints == index), 1:(index + 1)]
      changepoint_curr <- changepoint[which(num.changepoints == index), 1:index]

      indiv.log.lik <- matrix(NA, nrow = nrow(lambda_curr), ncol = nrow(df))

      for (x in 1:nrow(lambda_curr)) {
        indiv.log.lik[x, ] <- piecewise_loglik.indiv(df, as.numeric(data.frame(changepoint_curr)[x, ]), lambda_curr[x, ])
      }

      indiv.log.lik.final <- rbind(indiv.log.lik.final, indiv.log.lik)
    }
  }

  indiv.log.lik.final
}

add_km <- function(plt, df, colour = "black") {
  result.km <- survfit(Surv(time, status) ~ 1, data = df)
  km.data <- data.frame(cbind(result.km[[c("time")]], result.km[[c("surv")]], result.km[[c("upper")]], result.km[[c("lower")]]))
  colnames(km.data) <- c("time", "survival", "upper", "lower")

  plt +
    geom_step(data = km.data, aes(x = time, y = survival),colour = colour, inherit.aes = F) +
    geom_step(data = km.data, aes(x = time, y = upper), colour = colour,linetype = "dashed", inherit.aes = F) +
    geom_step(data = km.data, aes(x = time, y = lower), colour = colour, linetype = "dashed", inherit.aes = F)
}

get_Surv <- function(object, chng.num = "all", max_predict) {
  interval <- max_predict/100
  time <- c(seq(from = 0, to = max_predict, by = interval))
  k <- object$k.stacked
  lambda <- object$lambda
  changepoint <- object$changepoint
  num.changepoints <- unlist(apply(k, 1, function(x) {
    length(na.omit(x))
  }))

  if (chng.num != "all") {
    if (length(which(num.changepoints == chng.num)) < 2) {
      print("Too few simulations for this change-point model")
      break
    } else {
      changepoints_eval <- chng.num
    }
  } else {
    changepoints_eval <- unique(num.changepoints)[order(unique(num.changepoints))]
  }

  St <- NULL
  for (i in 1:length(changepoints_eval)) {
    if (changepoints_eval[i] == 0) {
      num_zero <- length(which(num.changepoints == changepoints_eval[i]))
      lambda_0 <- data.frame(lambda[which(num.changepoints == changepoints_eval[i]), 1])
      St <- cbind(St, surv_nochange(time, num_zero, lambda_0[, 1]))
    } else {
      k_curr <- data.frame(k[which(num.changepoints == changepoints_eval[i]), 1:changepoints_eval[i]])
      changepoint_curr <- data.frame(changepoint[which(num.changepoints == changepoints_eval[i]), 1:changepoints_eval[i]])
      lambda_curr <- lambda[which(num.changepoints == changepoints_eval[i]), 1:(changepoints_eval[i] + 1)]
      changepoint_df <- cbind(0, changepoint_curr[, 1:changepoints_eval[i]])
      time.interval_df <- t(apply(changepoint_df, 1, diff))

      if (changepoints_eval[i] == 1) {
        cum_haz_df <- time.interval_df * lambda_curr[, 1]
        surv_df <- cbind(1, t(exp(-cum_haz_df)))
      } else {
        # head(index2, -1) last hazard not needed for cumhaz calc
        cum_haz_df <- t(apply(time.interval_df * lambda_curr[, 1:changepoints_eval[i]], 1, cumsum))
        surv_df <- cbind(1, exp(-cum_haz_df))
      }
      St <- cbind(St, surv_change(time, nrow(k_curr), lambda_curr, data.matrix(changepoint_df), surv_df))
    }
  }
  rownames(St) <- time

  St
}

#' Printing of Changepoint model objects
#'
#' @param object of class "changepoint".
#' @param chng.num indicating the change-point model to summarized in terms of change-point(s) and hazards, if ommited the model with the highest posterior probability is presented.
#' @param digits number of digits to be printed.
#'
#' @export
#'
#' @examples \notrun{
#' Print a summary of the one change-point model
#' print(Collapsing_Model, mod = 1)
#' Print a summary of the model with the highest posterior probability
#' print(Collapsing_Model)
#' }
print.changepoint <- function(object, chng.num = NULL, digits = min(3L, getOption("digits"))) {
  cat("Posterior Change-point Probabilities:\n")

  names(attr(object$prob.changepoint, "dimnames")) <- NULL

  print.default(format(object$prob.changepoint, digits = digits),
    print.gap = 2L,
    quote = FALSE
  )

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

  lambda_summary <- summary(lambda_curr)
  colnames(lambda_summary) <- paste0(rep("lambda_", chng.prob + 1), 1:(chng.prob + 1))
  if(chng.prob != 0){
  changepoint_summary <- summary(changepoint_curr)
  colnames(changepoint_summary) <- paste0(rep("changepoint_", chng.prob), 1:(chng.prob))
  }
  cat(paste0("\n"))
  cat(paste0("Summary of ", chng.prob, " change-point model:\n"))
  cat(paste0("\n"))

  if(chng.prob != 0){
    print.output <- cbind(changepoint_summary, lambda_summary)
  }else{
    print.output <- lambda_summary
  }

  print.default(format(print.output, digits = digits),
    print.gap = 2L,
    quote = FALSE)
}

#' Fitting Bayesian Parametric Models with JAGS
#'
#' In order to compare the change-point model with other common parametric survival models using Pseudo-Marginal Likelihood (PML) and Widely Applicable Information Criterion (WAIC) we need to fit them in a Bayesian framework.
#' Requires Just Another Gibbs Sampler (JAGS) along with the packages \code{\link[R2jags]{rjags}} and \code{\link[loo]{waic}} to run. The JAGS models that are produced by this function should be assessed for convergence. Additionally the chains may need to be run longer.
#' The following models are fit (Note that the parameterization used in JAGS is not equivalent to the \code{\link[flexsurv]{flexsurvreg}} parameterization for the Weibull, Log-Logistic and Generalized Gamma):
#' \itemize{
#'   \item \strong{Exponential}
#'   \item \strong{Weibull}
#'   \item \strong{Log-Normal}
#'   \item \strong{Log-Logistic}
#'   \item \strong{Gompertz}
#'   \item \strong{Generalized Gamma}
#' }
#' @param df standard dataframe for time-to-event data. Two columns required, time (to event or censoring) and status (indicating event or censoring).
#' @param max_predict maximum survival time to be predicted from the JAGS model. Default is 10, however, depending on the timescale this should be changed.
#' @return A list of with the following items:
#'  \itemize{
#'   \item \strong{model.fit}: A dataframe with the PML and WAIC for the six parametric models fitted by JAGS.
#'   \item \strong{jags.models}: A list containing the posterior simulations of the 6 JAGS models (fit using the \code{\link[R2jags]{rjags}} function).
#'   \item \strong{jags.surv}: A list of the survival probabilities for the prespecified times from the 6 JAGS models.
#' }
#' @export
#' @examples \dontrun{
#' mod.comp.jags <-fit_surv_models(df)
#' }
fit_surv_models <- function(df, max_predict = 10) {
  interval <- max_predict/100
  cat(" \n")
  if ("rjags" %in% rownames(installed.packages()) == FALSE) {
    print("Need JAGS to run this function")
    break
  }


  if ("R2jags" %in% rownames(installed.packages()) == FALSE) {
    print("Need R2jags package to run this function")
    break
  }

  if ("loo" %in% rownames(installed.packages()) == FALSE) {
    print("Need loo package to evaluate WAIC")
    break
  }

  require(rjags)
  require(R2jags)
  require(loo)

  cat(crayon::blue("Fitting parametric models with JAGS... can take several minutes \n"))


  # Exponential Model

  expo <- "model{
  for(i in 1:N){
    is.censored[i]~dinterval(t[i],t.cen[i])
    t[i] ~ dexp(lambda)
    Like[i] <- ifelse(is.censored[i], 1- pexp(t.cen[i],lambda), dexp(t[i], lambda))
    invLik[i] <- 1/Like[i]

  }
  for(i in 1:length(t_pred)){
    St_pred[i] <- 1- pexp(t_pred[i],lambda)
  }
  lambda ~ dgamma(1,1)

}"

  # Weibull Model


  weibull <- "model{
  for(i in 1:N){
is.censored[i]~dinterval(t[i],t.cen[i])
t[i] ~ dweib(v,lambda)
Like[i] <- ifelse(is.censored[i], 1- pweib(t.cen[i],v,lambda), dweib(t[i],v, lambda))
invLik[i] <- 1/Like[i]

  }
for(i in 1:length(t_pred)){
St_pred[i] <- 1- pweib(t_pred[i],v,lambda)
}
lambda ~ dgamma(1,1)
v ~ dgamma(1,1)
}"

  # Weibull Model


  gamma.jags <- "model{
  for(i in 1:N){
is.censored[i]~dinterval(t[i],t.cen[i])
t[i] ~ dgamma(shape,lambda)
Like[i] <- ifelse(is.censored[i], 1- pgamma(t.cen[i],shape,lambda), dgamma(t[i],shape, lambda))
invLik[i] <- 1/Like[i]
  }
 for(i in 1:length(t_pred)){
    St_pred[i] <- 1- pgamma(t_pred[i],shape,lambda)
  }

lambda ~ dgamma(1,1)
shape ~ dgamma(1,1)
}"

  # Log-Normal Model

  lnorm.jags <- "model{
  for(i in 1:N){
is.censored[i]~dinterval(t[i],t.cen[i])
t[i] ~ dlnorm(mu,tau)
Like[i] <- ifelse(is.censored[i], 1- plnorm(t.cen[i],mu,tau), dlnorm(t[i],mu, tau))
invLik[i] <- 1/Like[i]
  }
 for(i in 1:length(t_pred)){
    St_pred[i] <- 1- plnorm(t_pred[i],mu,tau)
  }
mu ~ dnorm(0,0.1)
sd ~ dunif(0,5)
tau <- pow(sd,-2)
}"

  # Loglogistic Model

  llogis.jags <- "
model{
for(i in 1:N){
is.censored[i]~dinterval(t.log[i],t.cen.log[i])
t.log[i] ~ dlogis(mu,tau)
Like[i] <- ifelse(is.censored[i], 1/(1 + pow(exp(t.cen.log[i])/beta, alpha)),
          (alpha/beta)*pow(exp(t.log[i])/beta, alpha-1)/pow(1 + pow(exp(t.log[i])/beta,alpha),2))
invLik[i] <- 1/Like[i]

}
 for(i in 1:length(t_pred)){
    St_pred[i] <- 1/(1 + pow(t_pred[i]/beta, alpha))
  }

mu ~ dnorm(0,0.1)
scale ~ dgamma(1,1)
tau <- pow(scale,-1) # Inverse of scale which is beta on the log-logistic dist
beta <- exp(mu)
alpha <- tau

}"

  # Gompertz Model
  gompertz.jags <- "
data{
for(i in 1:N){
zero[i] <- 0}
}

model{

C <- 10000
for(i in 1:N){

logHaz[i] <- (log(b)+ a*time[i])*status[i]
logSurv[i] <- (-b/a)*(exp(a*time[i])-1)

LL[i] <- logHaz[i]+ logSurv[i]
Like[i] <- exp(LL[i])
invLik[i] <-pow(Like[i],-1)


zero[i] ~ dpois(zero.mean[i])
zero.mean[i] <- -logHaz[i]-logSurv[i] + C
}

for(i in 1:length(t_pred)){
    St_pred[i] <- exp((-b/a)*(exp(a*t_pred[i])-1))
  }


a ~ dnorm(0,0.01)
b ~ dunif(0,5)
}"

  # Generalized Gamma Model

  gen.gamma.jags <- "model{
  for(i in 1:N){
is.censored[i]~dinterval(t[i],t.cen[i])
t[i] ~ dgen.gamma(r,lambda,b)
Like[i] <- ifelse(is.censored[i], 1- pgen.gamma(t.cen[i],r,lambda,b), dgen.gamma(t[i],r,lambda,b))
invLik[i] <- 1/Like[i]

  }

 for(i in 1:length(t_pred)){
    St_pred[i] <- 1- pgen.gamma(t_pred[i],r,lambda,b)
  }
r ~ dgamma(1,1)
lambda ~ dgamma(1,1)
b ~ dgamma(1,1)
}"


  # Data

  data_new <- list()
  df_jags <- df[, c("time", "status")]
  df_jags$t <- df$time


  tinits1 <- df_jags$t + max(df$time)
  is.na(tinits1) <- df_jags$status == 1
  tinits2 <- tinits1 + 5

  is.na(df_jags$t) <- df_jags$status == 0
  df_jags$is.censored <- 1 - df_jags$status
  df_jags$t.cen <- df_jags$time + df_jags$status


  modelinits <- list(
    list(t = tinits1),
    list(t = tinits2)
  )

  logmodelinits <- list(
    list(t.log = log(tinits1)),
    list(t.log = log(tinits2))
  )

  data_jags <- list(
    N = nrow(df_jags),
    t.cen = df_jags$t.cen,
    is.censored = df_jags$is.censored,
    t = df_jags$t
  )

  data_jags$t_pred <- seq(0, max_predict, by = interval)


  ## Different data format needed for loglogistic

  data_jags_llogis <- data_jags
  data_jags_llogis$t.log <- log(data_jags$t)
  data_jags_llogis$t.cen.log <- log(data_jags$t.cen)

  `%!in%` <- Negate(`%in%`)
  data_jags_llogis <- data_jags_llogis[names(data_jags_llogis) %!in% c("t", "t.cen")]


  # Zeros trick for the Gompertz distribution since it is not in JAGS so different data format required

  data_gomp <- list()
  data_gomp$time <- df$time
  data_gomp$status <- df$status
  data_gomp$N <- nrow(df)

  data_gomp$t_pred <- data_jags$t_pred
  # Fit the Models

  cat(crayon::blue("Exponential Model \n"))

  expo.mod <- R2jags::jags(
    model.file = textConnection(expo),
    data = data_jags,
    inits = modelinits,
    n.chains = 2,
    parameters.to.save = c("Like", "lambda", "invLik", "St_pred")
  )




  # Same as approach PML.expo above
  # PML.expo <- nrow(Like.sims.expo)/colSums(1/Like.sims.expo)

  cat(crayon::blue("Weibull Model \n"))

  weibull.mod <- R2jags::jags(
    model.file = textConnection(weibull),
    data = data_jags,
    inits = modelinits,
    n.chains = 2,
    parameters.to.save = c("lambda", "v", "Like", "invLik", "St_pred")
  )

  cat(crayon::blue("Gamma Model \n"))
  gamma.mod <- R2jags::jags(
    model.file = textConnection(gamma.jags),
    data = data_jags,
    inits = modelinits,
    n.chains = 2,
    parameters.to.save = c("lambda", "shape", "Like", "invLik", "St_pred")
  )

  cat(crayon::blue("LogNormal Model \n"))
  lnorm.mod <- R2jags::jags(
    model.file = textConnection(lnorm.jags),
    data = data_jags,
    inits = modelinits,
    n.chains = 2,
    parameters.to.save = c("mu", "sd", "Like", "invLik", "St_pred")
  )

  cat(crayon::blue("LogLogistic Model \n"))
  llogis.mod <- R2jags::jags(
    model.file = textConnection(llogis.jags),
    data = data_jags_llogis,
    inits = logmodelinits,
    n.chains = 2,
    parameters.to.save = c("alpha", "beta", "Like", "invLik", "St_pred")
  )

  cat(crayon::blue("Gompertz Model \n"))
  gompertz.mod <- R2jags::jags(
    model.file = textConnection(gompertz.jags),
    data = data_gomp,
    inits = list(list(a = -1, b = 1 / mean(df$time)), list(a = 1, b = 1 / mean(df$time))),
    n.chains = 2,
    parameters.to.save = c("a", "b", "Like", "invLik", "St_pred")
  )

  cat(crayon::blue("Generalized Gamma Model \n"))
  gen.gamma.mod <- R2jags::jags(
    model.file = textConnection(gen.gamma.jags),
    data = data_jags,
    inits = modelinits,
    n.chains = 2,
    parameters.to.save = c("r", "lambda", "b", "Like", "invLik", "St_pred")
  )


  PML.expo <- 1 / expo.mod$BUGSoutput[["summary"]][grep("invLik", rownames(expo.mod$BUGSoutput[["summary"]])), 1]
  PML.expo.trans <- sum(log(PML.expo)) * (-2)
  Like.sims.expo <- expo.mod$BUGSoutput[["sims.matrix"]][, grep("Like", colnames(expo.mod$BUGSoutput[["sims.matrix"]]))]
  WAIC.expo <- waic(log(Like.sims.expo))
  Surv.expo <- expo.mod$BUGSoutput[["summary"]][grep("St_pred", rownames(expo.mod$BUGSoutput[["summary"]])), 1]


  PML.weib <- 1 / weibull.mod$BUGSoutput[["summary"]][grep("invLik", rownames(weibull.mod$BUGSoutput[["summary"]])), 1]
  PML.weib.trans <- sum(log(PML.weib)) * (-2)
  Like.sims.weib <- weibull.mod$BUGSoutput[["sims.matrix"]][, grep("Like", colnames(weibull.mod$BUGSoutput[["sims.matrix"]]))]
  WAIC.weib <- waic(log(Like.sims.weib))
  Surv.weib <- weibull.mod$BUGSoutput[["summary"]][grep("St_pred", rownames(weibull.mod$BUGSoutput[["summary"]])), 1]


  PML.gamma <- 1 / gamma.mod$BUGSoutput[["summary"]][grep("invLik", rownames(gamma.mod$BUGSoutput[["summary"]])), 1]
  PML.gamma.trans <- sum(log(PML.gamma)) * (-2)
  Like.sims.gamma <- gamma.mod$BUGSoutput[["sims.matrix"]][, grep("Like", colnames(gamma.mod$BUGSoutput[["sims.matrix"]]))]
  WAIC.gamma <- waic(log(Like.sims.gamma))
  Surv.gamma <- gamma.mod$BUGSoutput[["summary"]][grep("St_pred", rownames(gamma.mod$BUGSoutput[["summary"]])), 1]



  PML.lnorm <- 1 / lnorm.mod$BUGSoutput[["summary"]][grep("invLik", rownames(lnorm.mod$BUGSoutput[["summary"]])), 1]
  PML.lnorm.trans <- sum(log(PML.lnorm)) * (-2)
  Like.sims.lnorm <- lnorm.mod$BUGSoutput[["sims.matrix"]][, grep("Like", colnames(lnorm.mod$BUGSoutput[["sims.matrix"]]))]
  WAIC.lnorm <- waic(log(Like.sims.lnorm))
  Surv.lnorm <- lnorm.mod$BUGSoutput[["summary"]][grep("St_pred", rownames(lnorm.mod$BUGSoutput[["summary"]])), 1]


  PML.llogis <- 1 / llogis.mod$BUGSoutput[["summary"]][grep("invLik", rownames(llogis.mod$BUGSoutput[["summary"]])), 1]
  PML.llogis.trans <- sum(log(PML.llogis)) * (-2)
  Like.sims.llogis <- llogis.mod$BUGSoutput[["sims.matrix"]][, grep("Like", colnames(llogis.mod$BUGSoutput[["sims.matrix"]]))]
  WAIC.llogis <- waic(log(Like.sims.llogis))
  Surv.llogis <- llogis.mod$BUGSoutput[["summary"]][grep("St_pred", rownames(llogis.mod$BUGSoutput[["summary"]])), 1]


  PML.gomp <- 1 / gompertz.mod$BUGSoutput[["summary"]][grep("invLik", rownames(gompertz.mod$BUGSoutput[["summary"]])), 1]
  PML.gomp.trans <- sum(log(PML.gomp)) * (-2)
  Like.sims.gomp <- gompertz.mod$BUGSoutput[["sims.matrix"]][, grep("Like", colnames(gompertz.mod$BUGSoutput[["sims.matrix"]]))]
  WAIC.gomp <- waic(log(Like.sims.gomp))
  Surv.gomp <- gompertz.mod$BUGSoutput[["summary"]][grep("St_pred", rownames(gompertz.mod$BUGSoutput[["summary"]])), 1]


  PML.gen.gamma <- 1 / gen.gamma.mod$BUGSoutput[["summary"]][grep("invLik", rownames(gen.gamma.mod$BUGSoutput[["summary"]])), 1]
  PML.gen.gamma.trans <- sum(log(PML.gen.gamma)) * (-2)
  Like.sims.gen.gamma <- gen.gamma.mod$BUGSoutput[["sims.matrix"]][, grep("Like", colnames(gen.gamma.mod$BUGSoutput[["sims.matrix"]]))]
  WAIC.gen.gamma <- waic(log(Like.sims.gen.gamma))
  Surv.gen.gamma <- gen.gamma.mod$BUGSoutput[["summary"]][grep("St_pred", rownames(gen.gamma.mod$BUGSoutput[["summary"]])), 1]




  model.fit <- data.frame(
    Model = c("Exponential", "Weibull", "Gamma", "Log-Normal", "Log-Logistic", "Gompertz", "Generalized Gamma"),
    minustwo_logPML = c(PML.expo.trans, PML.weib.trans, PML.gamma.trans, PML.lnorm.trans, PML.llogis.trans, PML.gomp.trans, PML.gen.gamma.trans),
    WAIC = c(WAIC.expo$estimates[3, 1], WAIC.weib$estimates[3, 1], WAIC.gamma$estimates[3, 1], WAIC.lnorm$estimates[3, 1], WAIC.llogis$estimates[3, 1], WAIC.gomp$estimates[3, 1], WAIC.gen.gamma$estimates[3, 1])
  )


  jags_output <- list(
    model.fit = model.fit,
    jags.models = list(
      expo.mod,
      weibull.mod,
      gamma.mod,
      lnorm.mod,
      llogis.mod,
      gompertz.mod,
      gen.gamma.mod
    ),
    jags.surv = list(
      Surv.expo = Surv.expo,
      Surv.weib = Surv.weib,
      Surv.gamma = Surv.gamma,
      Surv.lnorm = Surv.lnorm,
      Surv.llogis = Surv.llogis,
      Surv.lnorm = Surv.lnorm,
      Surv.gomp = Surv.gomp,
      Surv.gen.gamma = Surv.gen.gamma
    )
  )

  return(jags_output)


  # PML.gomp <- nrow(Like.sims.gomp)/colSums(1/Like.sims.gomp)
  # flexsurvreg(Surv(time,status)~1,data = df, dist = "llogis")
  # flexsurvreg(Surv(time,status)~1,data = df, dist = "gompertz")

  # inits for flexsurvreg
  # flexsurv:::flexsurv.dists
}


#' Comparing Piecewise Exponential model with other Parametric models
#'
#' Compares the piecewise exponential model with 6 other standard parameteric models (see \code{\link[=fit_surv_models]{fit_surv_models()}}). Note that exponential is a special case of the piecewise exponential, however, it is refit in JAGS to highlight the difference in statistical fit.
#' This functions computes the individual log-likelihood for the piecewise exponential model (see \code{\link[=get.loglik]{get.loglik()}}) and compares the Widely Applicable Information Criterion (using the \code{\link[loo:loo-package]{loo::loo-package}}) and Pseudo-Marginal Likelihood (PML) with the other standard parametric models.
#'
#' @param object of class "changepoint".
#' @param max_predict maximum survival time to be predicted for the survival models. Default is 10, however, depending on the timescale this should be changed.
#'
#' @return A list of with the following items:
#'  \itemize{
#'   \item \strong{model.comp}: A dataframe with the PML and WAIC for the piecewise exponential model and the six parametric models fitted by JAGS.
#'   \item \strong{jags.models}: A list containing the posterior simulations of the 6 JAGS models (fit using the R2jags:jags function).
#'   \item \strong{plot_Surv_all}: A ggplot with the posterior mean survival probabilities for the time specified by max_predict.
#' }
#' @importFrom loo waic
#' @export
#' @md
#' @examples \dontrun{
#' mod.comp <-compare.surv.mods(Collapsing_Model, df)}
compare.surv.mods <- function(object, max_predict = 10) {
  interval <- max_predict/100
  df <- object$df
  cat(crayon::blue("Evaluating Individual log-likelihood for changepoint model \n ... can take several minutes"))

  log.lik.piece <- get.loglik(object)

  indiv.lik.piece <- exp(log.lik.piece)
  PML.indiv.piece <- 1 / indiv.lik.piece
  PML.piece <- nrow(PML.indiv.piece) / colSums(PML.indiv.piece)

  minus2logPML.piece <- -2 * sum(log(PML.piece))
  WAIC.piece <- waic(log.lik.piece)$estimate[3, 1]

  jags_output <- fit_surv_models(df, max_predict = max_predict)

  piecewise.mod.fit <- data.frame(
    Model = "Piecewise Exponential",
    minustwo_logPML = minus2logPML.piece,
    WAIC = WAIC.piece
  )
  mod.comp <- rbind(jags_output$model.fit, piecewise.mod.fit)
  colnames(mod.comp) <- c("Model", "-2log(PML)", "WAIC")


  plot_surv <- plot(object, add.post = F, max_predict = max_predict, interval = interval)

  t_pred <- seq(0, to = max_predict, by = interval)

  df_surv_expo <- data.frame(jags_output$jags.surv$Surv.expo, t_pred)
  df_surv_weib <- data.frame(jags_output$jags.surv$Surv.weib, t_pred)
  df_surv_gamma <- data.frame(jags_output$jags.surv$Surv.gamma, t_pred)
  df_surv_llogis <- data.frame(jags_output$jags.surv$Surv.llogis, t_pred)
  df_surv_lnorm <- data.frame(jags_output$jags.surv$Surv.lnorm, t_pred)
  df_surv_gomp <- data.frame(jags_output$jags.surv$Surv.gomp, t_pred)
  df_surv_gen.gamma <- data.frame(jags_output$jags.surv$Surv.gen.gamma, t_pred)



  colors <- c(
    "KM curve" = "black",
    "Piecewise Expo" = "purple",
    "Exponential" = "red",
    "Weibull" = "cyan",
    "Gamma" = "brown",
    "Log-Logistic" = "blue",
    "Log-Normal" = "pink",
    "Gompertz" = "green",
    "Gen. Gamma" = "orange"
  )
  # https://stackoverflow.com/questions/13701347/force-the-origin-to-start-at-0



  plot_Surv_all <- plot_surv +
    geom_line(df_surv_expo, mapping = aes(y = jags_output.jags.surv.Surv.expo, x = t_pred, col = "Exponential"), inherit.aes = F) +
    geom_line(df_surv_weib, mapping = aes(y = jags_output.jags.surv.Surv.weib, x = t_pred, col = "Weibull"), inherit.aes = F) +
    geom_line(df_surv_gamma, mapping = aes(y = jags_output.jags.surv.Surv.gamma, x = t_pred, col = "Gamma"), inherit.aes = F) +
    geom_line(df_surv_llogis, mapping = aes(y = jags_output.jags.surv.Surv.llogis, x = t_pred, col = "Log-Logistic"), inherit.aes = F) +
    geom_line(df_surv_lnorm, mapping = aes(y = jags_output.jags.surv.Surv.lnorm, x = t_pred, col = "Log-Normal"), inherit.aes = F) +
    geom_line(df_surv_gomp, mapping = aes(y = jags_output.jags.surv.Surv.gomp, x = t_pred, col = "Gompertz"), inherit.aes = F) +
    geom_line(df_surv_gen.gamma, mapping = aes(y = jags_output.jags.surv.Surv.gen.gamma, x = t_pred, col = "Gen. Gamma"), inherit.aes = F) +
    scale_color_manual(values = colors) +
    labs(
      x = "Time",
      y = "Survival",
      color = "Survival Functions"
    ) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.8)) +
    geom_vline(xintercept = c(max(df$time)), linetype = "dotted")


  return(list(mod.comp = mod.comp, jag.models = jags_output$jags.models, plot_Surv_all = plot_Surv_all))
}
