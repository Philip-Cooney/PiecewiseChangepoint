
#Helper integration function... maybe add to the package
integrate.xy <- function(x,fx, a,b, use.spline = TRUE, xtol = 2e-8){
  if(is.list(x)) {
    fx <- x$y; x <- x$x
    if(length(x) == 0)
      stop("list 'x' has no valid $x component")
  }
  if((n <- length(x)) != length(fx))
    stop("'fx' must have same length as 'x'")
  
  if(is.unsorted(x)) { i <- sort.list(x); x <- x[i]; fx <- fx[i] }
  if(any(i <- duplicated(x))) {
    n <- length(x <- x[!i])
    ## we might have to check that the same fx[] are duplicated
    ## otherwise either give an error or take the mean() of those...
    fx <- fx[!i]
  }
  if(any(diff(x) == 0))
    stop("bug in 'duplicated()' killed me: have still multiple x[]!")
  
  if(missing(a)) a <- x[1]
  else if(any(a < x[1])) stop("'a' must NOT be smaller than min(x)")
  if(missing(b)) b <- x[n]
  else if(any(b > x[n])) stop("'b' must NOT be larger  than max(x)")
  if(length(a) != 1 && length(b) != 1 && length(a) != length(b))
    stop("'a' and 'b' must have length 1 or same length !")
  else {
    k <- max(length(a),length(b))
    if(any(b < a))    stop("'b' must be elementwise >= 'a'")
  }
  
  if(use.spline) {
    xy <- spline(x,fx, n = max(1024, 3*n))
    ##-- Work around spline(.) BUG:  (ex.:  range(spline(1:20,1:20,n=95)))
    if(xy$x[length(xy$x)] < x[n]) {
      if(TRUE) cat("working around spline(.) BUG --- hmm, really?\n\n")
      xy$x <- c(xy$x,  x[n])
      xy$y <- c(xy$y, fx[n])
    }
    ## END if work around ----
    x <- xy$x; fx <- xy$y
    n <- length(x)
  }
  
  ab <- unique(c(a,b))
  BB <- abs(outer(x,ab,"-")) < (xtol * max(b - a))
  if(any(j <- 0 == colSums(BB))) { # the j-th element(s) of ab are not in x[]
    y <- approx(x,fx, xout = ab[j])$y
    x <- c(ab[j],x)
    i <- sort.list(x)
    x <- x[i];  fx <- c(y,fx)[i];  n <- length(x)
  }
  
  ##--- now we could use 'Simpson's formula IFF the x[i] are equispaced... --
  ##--- Since this may well be wrong, just use 'trapezoid formula':
  
  dig0 <- floor(-log10(xtol)) #
  f.match <- function(x,table,dig) match(signif(x,dig), signif(table,dig))
  ## was (S+) f.match <- function(x,table) match(as.single(x), as.single(table))
  
  d <- dig0; while(anyNA(ai <- f.match(a,x, d))) d <- d - 1/8 ; ai <- rep_len(ai, k)
  d <- dig0; while(anyNA(bi <- f.match(b,x, d))) d <- d - 1/8 ; bi <- rep_len(bi, k)
  dfx <- fx[-c(1,n)] * diff(x,lag = 2)
  r <- numeric(k)
  for (i in 1:k) {
    a <- ai[i];  b <- bi[i]
    r[i] <- (x[a+1] - x[a])*fx[a] + (x[b] - x[b-1])*fx[b] +
      sum(dfx[seq(a, length = max(0,b-a-1))])
  }
  r/2
}


rpwexp <- function(n, lam, s){
  
  U = runif(n, 0, 1)
  X = rep(NA,n)
  haz_seg <- diff(c(0,s))*lam[-length(lam)]
  cum_haz <- cumsum(haz_seg)
  St_thres <- exp(-cum_haz)
  #https://rdrr.io/cran/CPsurv/src/R/sim.survdata.R
  for(i in 1:n){
    int <- which(U[i] < St_thres)
    if(length(int) == 0){
      X[i] = qexp(U[i], rate = lam[1], lower.tail = F)
    }else{
      X[i] = s[max(int)] + qexp(U[i]/St_thres[max(int)], rate = lam[max(int)+1], lower.tail = F)
    }
    
  }
  return(X)
  
}


compare_dco <-function (all_surv_mods, old_dco, new_dco, km_risk = 0.1){
  result.km_old <- survfit(Surv(time, status) ~ 1, data = old_dco)
  if (!is.null(km_risk)) {
    max_time <- result.km_old$time[max(which(result.km_old$n.risk/result.km_old$n >= km_risk))]
  }
  
  result.km <- survfit(Surv(time, status) ~ 1, data = new_dco)
  km.data <- data.frame(cbind(result.km[[c("time")]], result.km[[c("surv")]], 
                              result.km[[c("upper")]], result.km[[c("lower")]]))
  colnames(km.data) <- c("time", "survival", "upper", "lower")
  if (is.null(km_risk)) {
    km.data <- km.data %>% filter(time >= max(old_dco$time))
  }
  else {
    km.data <- km.data %>% filter(time >= max_time)
  }
  all_surv_mods$plot_Surv_all + geom_step(data = km.data, aes(x = time, y = survival), 
                                          colour = "black", inherit.aes = F) + geom_step(data = km.data, 
                                                                                         aes(x = time, y = upper), colour = "black", linetype = "dashed", 
                                                                                         inherit.aes = F) + geom_step(data = km.data, aes(x = time, 
                                                                                                                                          y = lower), 
                                                                                                                      colour = "black", linetype = "dashed", inherit.aes = F)
}




SurvSplit <- function (Y, cuts){
  #Taken from eha
  if (NCOL(Y) == 2)
    Y <- cbind(rep(0, NROW(Y)), Y)
  indat <- cbind(Y, 1:NROW(Y), rep(-1, NROW(Y)))
  colnames(indat) <- c("enter", "exit", "event", "idx", "ivl")
  n <- length(cuts)
  cuts <- sort(cuts)
  if ((cuts[1] <= 0) || (cuts[n] == Inf))
    stop("'cuts' must be positive and finite.")
  cuts <- c(0, cuts, Inf)
  n <- n + 1
  out <- list()
  indat <- as.data.frame(indat)
  for (i in 1:n) {
    out[[i]] <- age.window(indat, cuts[i:(i + 1)])
    out[[i]]$ivl <- i
  }
  Y <- do.call(rbind, out)
  colnames(Y) <- colnames(indat)
  list(Y = Y[, 1:3], ivl = Y[, 5], idx = Y[, 4])
}

age.window <- function (dat, window, surv = c("enter", "exit", "event")){
  #Taken from eha
  if (!is.data.frame(dat))
    stop("dat must be a data frame")
  if (length(surv) != 3)
    stop("surv must have length 3")
  fixed.names <- names(dat)
  surv.indices <- match(surv, fixed.names)
  if (length(which(is.na(surv.indices)))) {
    x <- which(is.na(surv.indices))
    stop(paste(surv[x], " is not a name in the data frame."))
  }
  enter <- dat[[surv.indices[1]]]
  exit <- dat[[surv.indices[2]]]
  event <- dat[[surv.indices[3]]]
  who <- (exit > window[1]) & (enter < window[2])
  if (sum(who) > 0.5) {
    enter <- enter[who]
    exit <- exit[who]
    event <- event[who]
    event[exit > window[2]] <- 0
    exit[exit > window[2]] <- window[2]
    enter[enter < window[1]] <- window[1]
    dat <- dat[who, ]
    dat[surv.indices[1]] <- enter
    dat[surv.indices[2]] <- exit
    dat[surv.indices[3]] <- event
  }
  else {
    dat <- NULL
  }
  dat
}


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
  #ratemat <- matrix(rep(rate, n_obs / 2),
  #  nrow = n_obs,
  #  ncol = num.breaks + 1, byrow = TRUE
  #)
  
  if (n_cens_req > 0) {
    if (num.breaks == 0) {
      samp_cens <- rexp(n_cens_req * 2, rate)
      samp <- rexp(n_obs, rate)
    } else {
      samp_cens <- rpwexp(n_cens_req * 2, rate, t_change) # Assume that on average half the observations will be censors
      samp <- rpwexp(n_obs, rate,  t_change)
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
      samp <- rpwexp(n_obs, rate,  t_change)
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
plot.changepoint <- function (object, type = "survival", chng.num = "all", add.km = T, 
                              max_predict = 10, add.post = T, alpha.pos = NULL, t_pred = NULL,
                              final_chng_only =F,col_km = "black",km_risk = NULL,
                              ...){
  if (type == "survival") {
    St <- get_Surv(object, chng.num = chng.num, max_predict = max_predict, 
                   time = t_pred)
    return(plot.Survival(St, add.km = add.km, add.post = add.post, 
                         alpha.pos = alpha.pos))
  }
  if (type == "hazard") {
    return(plot.hazard(object, chng.num = chng.num, alpha.pos = alpha.pos))
  }
}


plot.Survival<- function (St, max.num.post = 500, add.km, add.post, alpha.pos, 
                          env = parent.frame(), final_chng_only = F,col_km = "black", km_risk = NULL){
  nSims <- ncol(St)
  time <- as.numeric(rownames(St))
  mean.Surv_df <- data.frame(survival = apply(St, 1, FUN = mean), 
                             time = time)
  mod_quantile_df <- data.frame(cbind(time, t(apply(St, 1, 
                                                    FUN = quantile, c(0.025, 0.975)))))
  colnames(mod_quantile_df) <- c("time", "lower", "upper")
  Surv.plot <- data.frame(Survival = c(unlist(St)), time = rep(time, 
                                                               nSims), id = rep(1:nSims, each = length(time)))
  if (max.num.post < nSims) {
    post_id <- sample(1:nSims, size = max.num.post)
    Surv.plot <- dplyr::filter(Surv.plot, id %in% post_id)
  }
  plot_Surv <- ggplot(data = Surv.plot, mapping = aes(x = time, 
                                                      y = Survival, group = id)) + geom_line(data = mean.Surv_df, 
                                                                                             aes(x = time, y = survival), size = 1, inherit.aes = F, 
                                                                                             colour = "purple") + scale_y_continuous(breaks = seq(0,1, by = 0.1),
                                                                                                                                     expand = c(0, 0)) + 
    scale_x_continuous(expand = c(0,0)) + annotate(geom = "segment", x = seq(0, max(time),max(time)/50), xend = seq(0, max(time), max(time)/50),
                                                   y = 0, yend = 0.01) + theme_classic()
  if (add.post == T) {
    if (is.null(alpha.pos)) {
      alpha.pos <- 0.025
    }
    else {
      alpha.pos <- alpha.pos
    }
    plot_Surv <- plot_Surv + geom_line(data = mod_quantile_df, 
                                       aes(x = time, y = lower), linetype = "dashed", size = 1, 
                                       inherit.aes = F, colour = "grey") + 
      geom_line(data = mod_quantile_df, aes(x = time, y = upper), linetype = "dashed", size = 1, 
                inherit.aes = F, colour = "grey") + geom_line(size = 0.1,alpha = alpha.pos, colour = "red")
  }
  if (env$chng.num != "all" && env$chng.num != 0) {
    k <- env$object$k.stacked
    num.changepoints <- unlist(apply(k, 1, function(x) {
      length(na.omit(x))
    }))
    k_curr <- data.frame(k[which(num.changepoints == env$chng.num), 
                           1:env$chng.num])
    df <- env$object$df
    df_event <- unique(df[which(df$status == 1), c("status", 
                                                   "time")])
    time.break <- df_event[apply(k_curr, 2, FUN = mean), 
                           "time"]
    survival.close <- sapply(time.break, FUN = function(x) {
      which.min(abs(mean.Surv_df$time - x))
    })
    break.points.Surv <- data.frame(time = mean.Surv_df$time[survival.close], Survival = mean.Surv_df$survival[survival.close])
    
    if(env$final_chng_only){
      break.points.Surv <- break.points.Surv[nrow(break.points.Surv),,drop = F]
    }
    
    plot_Surv <- plot_Surv + geom_point(data = break.points.Surv, 
                                        aes(x = time, Survival), shape = 23, fill = "green", 
                                        color = "darkred", size = 5, inherit.aes = F, stroke  = 2.5)
  }
  if (add.km) {
    plot_Surv <- add_km(plot_Surv, env$object$df, colour = env$col_km, km_risk  = env$km_risk)
  }
  return(plot_Surv)
}

plot.hazard <- function(object, chng.num = "all", max.num.post = 500, alpha.pos = NULL, ...) {
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
  
  if(is.null(alpha.pos)){
    alpha.pos <- 0.025
  }else{
    alpha.pos <- alpha.pos
  }
  
  plot_haz <- ggplot(df.plot.final, aes(timepoints, hazards)) +
    geom_step(aes(group = id), linetype = "dashed", alpha = alpha.pos, colour = "red") +
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
  split <- SurvSplit(surv.object, changepoint)
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
#'
get.loglik <- function(df,lambda_df,changepoint_df ){
  df_event <- df %>% filter(status == 1)
  
  Surv_mat <- haz_mat <- matrix(nrow = nrow(lambda_df), ncol = nrow(df))
  
  for(i in 1:nrow(lambda_df)){
    haz_mat[i,] <-   GetHazPEH(df$time,na.omit(changepoint_df[i,]),na.omit(lambda_df[i,]))
    Surv_mat[i,] <-    GetSurvPEH(df$time,na.omit(changepoint_df[i,]),na.omit(lambda_df[i,]))
    
  }
  
  log_haz_mat <- log(haz_mat)
  log_Surv_mat <- log(Surv_mat)
  log_haz_mat[,which(df$status == 0)] <- 0
  
  return(log_haz_mat + log_Surv_mat)
  
}

get.loglik_redundant <- function(object, chng.num = "all") { # Redundant because it is very slow
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
      lambda_curr <- lambda[which(num.changepoints == index), 1:(index + 1), drop = F]
      changepoint_curr <- changepoint[which(num.changepoints == index), 1:index, drop = F]
      
      indiv.log.lik <- matrix(NA, nrow = nrow(lambda_curr), ncol = nrow(df))
      
      for (x in 1:nrow(lambda_curr)) {
        indiv.log.lik[x, ] <- piecewise_loglik.indiv(df, as.numeric(data.frame(changepoint_curr)[x, ]), lambda_curr[x, ])
      }
      
      indiv.log.lik.final <- rbind(indiv.log.lik.final, indiv.log.lik)
    }
  }
  
  indiv.log.lik.final
}




GetHazPEH = function(x, s, lam) {
  #x is times
  #s is changepoints
  #lam is lambda (although in other of his equations they are log lambda)
  y = x
  J = length(s)
  s <- c(0,s,Inf)
  for (m in 1:length(x)) {
    for (k in 1:(J + 1)) {
      if ((x[m] > s[k]) && (x[m] <= s[k + 1])) {
        y[m] = lam[k]
      }
    }
  }
  return(y)
}

GetSurvPEH = function(x, s, lam) {
  y = x
  J = length(s)
  s <- c(0,s,Inf)
  
  for (m in 1:length(x)) {
    for (k in 1:(J + 1)) {
      if ((x[m] > s[k]) && (x[m] <= s[k + 1])) {
        if (k > 1) {
          y[m] = exp(-lam[k] * (y[m] - s[k]) -
                       sum(lam[1:(k - 1)] * (s[2:k] - s[1:(k-1)])))
        }
        else {
          y[m] = exp(-lam[k] * (y[m] - s[k]))
        }
      }
    }
  }
  return(y)
}

get.loglik_ind <- function(df,lambda_df,changepoint_df ){
  df_event <- df %>% filter(status == 1)
  
  Surv_mat <- haz_mat <- matrix(nrow = nrow(lambda_df), ncol = nrow(df))
  
  for(i in 1:nrow(lambda_df)){
    haz_mat[i,] <-   GetHazPEH(df$time,na.omit(changepoint_df[i,]),na.omit(lambda_df[i,]))
    Surv_mat[i,] <-    GetSurvPEH(df$time,na.omit(changepoint_df[i,]),na.omit(lambda_df[i,]))
    
  }
  
  log_haz_mat <- log(haz_mat)
  log_Surv_mat <- log(Surv_mat)
  log_haz_mat[,which(df$status == 0)] <- 0
  
  return(log_haz_mat + log_Surv_mat)
  
}




add_km <- function (plt, df, colour = "black", km_risk = NULL){
  result.km <- survfit(Surv(time, status) ~ 1, data = df)
  km.data <- data.frame(cbind(result.km[[c("time")]], result.km[[c("surv")]], 
                              result.km[[c("upper")]], result.km[[c("lower")]]))
  colnames(km.data) <- c("time", "survival", "upper", "lower")
  if (!is.null(km_risk)) {
    max_time <- result.km$time[max(which(result.km$n.risk/result.km$n >= 
                                           km_risk))]
    km.data <- km.data %>% filter(time <= max_time)
  }
  
  plt <- plt + geom_step(data = km.data, aes(x = time, y = survival), 
                         colour = colour, inherit.aes = F) + geom_step(data = km.data, 
                                                                       aes(x = time, y = upper), colour = colour, linetype = "dashed", 
                                                                       inherit.aes = F) + 
    geom_step(data = km.data, aes(x = time,y = lower), colour = colour, linetype = "dashed", inherit.aes = F) 
  if (!is.null(km_risk)) {
    plt+geom_vline(xintercept = max_time, linetype = "dotted")
  }else{
    plt
  }
}


get_Surv <- function(object, chng.num = "all", max_predict = NULL, time = NULL) {
  
  if(is.null(time)){
    interval <- max_predict/100
    time <- c(seq(from = 0, to = max_predict, by = interval))
  }
  
  #k <- object$k.stacked
  lambda <- object$lambda
  changepoint <- object$changepoint
  num.changepoints <- unlist(apply(changepoint, 1, function(x) {
    length(na.omit(x))
  }))
  
  if (chng.num != "all") {
    if (length(which(num.changepoints == chng.num)) < 2) {
      stop("Too few simulations for this change-point model")
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
      #k_curr <- data.frame(k[which(num.changepoints == changepoints_eval[i]), 1:changepoints_eval[i]])
      changepoint_curr <- data.frame(changepoint[which(num.changepoints == changepoints_eval[i]), 1:changepoints_eval[i],
                                                 drop=FALSE])
      lambda_curr <- lambda[which(num.changepoints == changepoints_eval[i]), 1:(changepoints_eval[i] + 1), drop=FALSE]
      changepoint_df <- cbind(0, changepoint_curr[, 1:changepoints_eval[i], drop=FALSE])
      time.interval_df <- t(apply(changepoint_df, 1, diff))
      
      if (changepoints_eval[i] == 1) {
        cum_haz_df <- time.interval_df * lambda_curr[, 1]
        surv_df <- cbind(1, t(exp(-cum_haz_df)))
      } else {
        # head(index2, -1) last hazard not needed for cumhaz calc
        cum_haz_df <- t(apply(time.interval_df * lambda_curr[, 1:changepoints_eval[i],drop=FALSE], 1, cumsum))
        surv_df <- cbind(1, exp(-cum_haz_df))
      }
      
      St <- cbind(St, surv_change(time, nrow(changepoint_curr), lambda_curr, data.matrix(changepoint_df), surv_df))
      
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
#' Requires Just Another Gibbs Sampler (JAGS) along with the packages \code{\link[R2jags]{rjags}},\code{\link[rstan]{rstan}}  and \code{\link[loo]{waic}} to run. The JAGS models that are produced by this function should be assessed for convergence. Additionally the chains may need to be run longer.
#' The following models are fit (Note that the parameterization used in JAGS is not equivalent to the \code{\link[flexsurv]{flexsurvreg}} parameterization for the Weibull, Log-Logistic and Generalized Gamma):
#' \itemize{
#'   \item \strong{Exponential}
#'   \item \strong{Weibull}
#'   \item \strong{Log-Normal}
#'   \item \strong{Log-Logistic}
#'   \item \strong{Gompertz}
#'   \item \strong{Generalized Gamma}
#'   \item \strong{Royston-Parmar Cubic Spline (1 or 2 know)}
#' }
#' 
#' Number of knots for Royston-Parmar is made by assessing, and finding which model gives the lowest WAIC.
#' 
#' @param df standard dataframe for time-to-event data. Two columns required, time (to event or censoring) and status (indicating event or censoring).
#' @param max_predict maximum survival time to be predicted from the survival models. Default is 10, however, depending on the timescale this should be changed.
#' @return A list of with the following items:
#'  \itemize{
#'   \item \strong{model.fit}: A dataframe with the PML and WAIC for the seven parametric models fitted by JAGS/Stan.
#'   \item \strong{jags.models}: A list containing the posterior simulations of the 6 JAGS models (fit using the \code{\link[R2jags]{rjags}} function).
#'   \item \strong{jags.surv}: A list of the survival probabilities for the prespecified times from the 7 JAGS/Stan models.
#' }
#' @export
fit_surv_models <- function (df, max_predict = 10, n.iter.jags = 2000, n.thin.jags = NULL, 
                             n.burnin.jags = NULL, gof, inc_waic = T, t_pred =NULL){
  if (is.null(n.burnin.jags)) {
    n.burnin.jags = floor(n.iter.jags/2)
  }
  if (is.null(n.thin.jags)) {
    n.thin.jags <- max(1, floor((n.iter.jags - n.burnin.jags)/1000))
  }
  
  if (is.null(t_pred)) {
    t_pred <- seq(0, to = max_predict, length.out =100 )
  }
  
  
  
  
  cat(" \n")
  if ("rjags" %in% rownames(installed.packages()) == FALSE) {
    stop("Need JAGS to run this function")
  }
  if ("R2jags" %in% rownames(installed.packages()) == FALSE) {
    stop("Need R2jags package to run this function")
  }
  if ("loo" %in% rownames(installed.packages()) == FALSE) {
    stop("Need loo package to evaluate WAIC")
  }
  if ("rstan" %in% rownames(installed.packages()) == FALSE) {
    stop("Need rstan package to evaluate Royston-Parmar models")
  }
  require("rjags")
  require("R2jags")
  require("loo")
  require("rstan")
  inits_list <- function(mod, n.chains = 2) {
    list_return <- list()
    for (i in 1:n.chains) {
      list_inits <- list()
      list_inits$t <- tinits1 + runif(1)
      if (mod == "exp") {
        list_inits$lambda = 1/mean(df$time)
      }
      if (mod == "weibull") {
        lt <- log(df$time[df$time > 0])
        shape <- 1.64/var(lt)
        scale <- exp(mean(lt) + 0.572)
        list_inits$v <- shape
        list_inits$lambda <- scale^{
          -shape
        }
      }
      if (mod == "gompertz") {
        list_inits$a = 0.001
        list_inits$b = 1/mean(df$time)
        list_inits <- list_inits[names(list_inits) %!in% 
                                   c("t")]
      }
      if (mod == "lnorm") {
        lt <- log(df$time[df$time > 0])
        list_inits$mu <- mean(lt)
        list_inits$sd <- sd(lt)
      }
      if (mod == "llogis") {
        lt <- log(df$time[df$time > 0])
        list_inits$mu <- mean(lt)
        list_inits$scale <- 3 * var(lt)/(pi^2)
        list_inits$t.log <- log(tinits1 + runif(1))
        list_inits <- list_inits[names(list_inits) %!in% 
                                   c("t")]
      }
      if (mod == "gengamma") {
        list_inits$r <- 1
        list_inits$lambda <- 1/mean(df$time)
        list_inits$b <- 1
      }
      if (mod == "gamma") {
        list_inits$lambda = sum(df$time)
        list_inits$shape = sum(df$status)
      }
      list_return[[i]] <- list_inits
    }
    return(list_return)
  }
  cat("Fitting parametric models with JAGS... can take several minutes \n")
  expo <- "model{\n  for(i in 1:N){\n    is.censored[i]~dinterval(t[i],t.cen[i])\n    t[i] ~ dexp(lambda)\n    Like[i] <- ifelse(is.censored[i], 1- pexp(t.cen[i],lambda), dexp(t[i], lambda))\n    #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)\n  }\n  for(i in 1:length(t_pred)){\n    St_pred[i] <- 1- pexp(t_pred[i],lambda)\n  }\n  lambda ~ dgamma(0.001,0.001)\n  total_LLik <- sum(log(Like))\n}"
  weibull <- "model{\n  for(i in 1:N){\n    is.censored[i]~dinterval(t[i],t.cen[i])\n    t[i] ~ dweib(v,lambda)\n    Like[i] <- ifelse(is.censored[i], 1- pweib(t.cen[i],v,lambda), dweib(t[i],v, lambda))\n    #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)\n  }\nfor(i in 1:length(t_pred)){\nSt_pred[i] <- 1- pweib(t_pred[i],v,lambda)\n}\nlambda ~ dgamma(0.001,0.001)\nv ~ dgamma(0.001,0.001)\n  total_LLik <- sum(log(Like))\n}"
  gamma.jags <- "model{\n  for(i in 1:N){\n    is.censored[i]~dinterval(t[i],t.cen[i])\n    t[i] ~ dgamma(shape,lambda)\n    Like[i] <- ifelse(is.censored[i], 1- pgamma(t.cen[i],shape,lambda), dgamma(t[i],shape, lambda))\n   #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)\n  }\n for(i in 1:length(t_pred)){\n    St_pred[i] <- 1- pgamma(t_pred[i],shape,lambda)\n  }\nlambda ~ dgamma(0.01,0.01)\nshape ~dgamma(0.01,0.01)\n  total_LLik <- sum(log(Like))\n}"
  lnorm.jags <- "model{\n  for(i in 1:N){\n    is.censored[i]~dinterval(t[i],t.cen[i])\n    t[i] ~ dlnorm(mu,tau)\n    Like[i] <- ifelse(is.censored[i], 1- plnorm(t.cen[i],mu,tau), dlnorm(t[i],mu, tau))\n   #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)\n  }\n for(i in 1:length(t_pred)){\n    St_pred[i] <- 1- plnorm(t_pred[i],mu,tau)\n  }\nmu ~ dnorm(0,0.001)\nsd ~ dunif(0,10)\ntau <- pow(sd,-2)\n  total_LLik <- sum(log(Like))\n}"
  llogis.jags <- "\nmodel{\nfor(i in 1:N){\n    is.censored[i]~dinterval(t.log[i],t.cen.log[i])\n    t.log[i] ~ dlogis(mu,tau)\n    Like[i] <- ifelse(is.censored[i], 1/(1 + pow(exp(t.cen.log[i])/beta, alpha)),\n          (alpha/beta)*pow(exp(t.log[i])/beta, alpha-1)/pow(1 + pow(exp(t.log[i])/beta,alpha),2))\n   #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)\n  }\n for(i in 1:length(t_pred)){\n    St_pred[i] <- 1/(1 + pow(t_pred[i]/beta, alpha))\n  }\nmu ~ dnorm(0,0.001)\nscale ~ dgamma(0.001,0.001)\ntau <- pow(scale,-1) # Inverse of scale which is beta on the log-logistic dist\nbeta <- exp(mu)\nalpha <- tau\n  total_LLik <- sum(log(Like))\n}"
  gompertz.jags <- "\ndata{\nfor(i in 1:N){\nzero[i] <- 0}\n}\nmodel{\nC <- 10000\nfor(i in 1:N){\nlogHaz[i] <- (log(b)+ a*time[i])*status[i]\nlogSurv[i] <- (-b/a)*(exp(a*time[i])-1)\nLL[i] <- logHaz[i]+ logSurv[i]\nLike[i] <- exp(LL[i])\n#invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)\nzero[i] ~ dpois(zero.mean[i])\nzero.mean[i] <- -logHaz[i]-logSurv[i] + C\n}\nfor(i in 1:length(t_pred)){\n    St_pred[i] <- exp((-b/a)*(exp(a*t_pred[i])-1))\n  }\na ~ dnorm(0,0.001)\nb ~ dunif(0,10)\n  total_LLik <- sum(log(Like))\n}"
  gen.gamma.jags <- "model{\n    for(i in 1:N){\n    is.censored[i]~dinterval(t[i],t.cen[i])\n    t[i] ~ dgen.gamma(r,lambda,b)\n    Like[i] <- ifelse(is.censored[i], 1- pgen.gamma(t.cen[i],r,lambda,b), dgen.gamma(t[i],r,lambda,b))\n   #invLik[i] <- 1/Like[i] Unstable for some datasets (Will calculate outside JAGS)\n  }\n for(i in 1:length(t_pred)){\n    St_pred[i] <- 1- pgen.gamma(t_pred[i],r,lambda,b)\n  }\n    r ~ dgamma(0.001,0.001)\n    lambda ~ dgamma(0.001,0.001)\n    b ~ dgamma(0.001,0.001)\n     total_LLik <- sum(log(Like))\n}"
  rps.stan <- "// Royston-Parmar splines model\n\nfunctions {\n  real rps_lpdf(vector t, vector d, vector gamma, matrix B, matrix DB, vector linpred) {\n    // t = vector of observed times\n    // d = event indicator (=1 if event happened and 0 if censored)\n    // gamma = M+2 vector of coefficients for the flexible part\n    // B = matrix of basis\n    // DB = matrix of derivatives for the basis\n    // linpred = fixed effect part\n    vector[num_elements(t)] eta;\n    vector[num_elements(t)] eta_prime;\n    vector[num_elements(t)] log_lik;\n    real lprob;\n    \n    eta = B*gamma + linpred;\n    eta_prime = DB*gamma;\n    log_lik = d .* (-log(t) + log(eta_prime) + eta) - exp(eta);\n    lprob = sum(log_lik);\n    return lprob;\n  }\n  \n  real Sind( vector gamma, row_vector B, real linpred) {\n    // t = vector of observed times\n    // gamma = M+2 vector of coefficients for the flexible part\n    // B = row_vector of basis\n    // linpred = fixed effect part\n    real eta;\n    real Sind_rtn;\n    \n    eta = B*gamma + linpred;\n    Sind_rtn = exp(-exp(eta));\n    return Sind_rtn;\n  }\n  \n  \n\n}\n\ndata {\n  int<lower=1> n;                   // number of observations\n  int<lower=0> M;                   // number of internal knots for the splines model\n  int<lower=1> H;                   // number of covariates in the (time-independent) linear predictor\n  vector<lower=0>[n] t;             // observed times (including censored values)\n  vector<lower=0,upper=1>[n] d;     // censoring indicator: 1 if fully observed, 0 if censored\n  matrix[n,H] X;                    // matrix of covariates for the (time-independent) linear predictor\n  matrix[n,M+2] B;                  // matrix with basis\n  matrix[n,M+2] DB;                 // matrix with derivatives of the basis\n  vector[H] mu_beta;                // mean of the covariates coefficients\n  vector<lower=0> [H] sigma_beta;   // sd of the covariates coefficients\n  vector[M+2] mu_gamma;             // mean of the splines coefficients\n  vector<lower=0>[M+2] sigma_gamma; // sd of the splines coefficients\n  \n}\n\n\nparameters {\n  vector[M+2] gamma;\n  vector[H] beta;\n}\n\n\ntransformed parameters{\n  vector[n] linpred;\n  vector[n] mu;\n\n  linpred = X*beta;\n  for (i in 1:n) {\n    mu[i] = linpred[i];\n  }\n\n}\n\nmodel {\n  // Priors\n  gamma ~ normal(mu_gamma,sigma_gamma);\n  beta ~ normal(mu_beta,sigma_beta);\n  \n  // Data model\n  t ~ rps(d,gamma,B,DB,X*beta);\n  \n}"
  rps.stan_mod <- rstan::stan_model(model_code = rps.stan)
  data_new <- list()
  df_jags <- df[, c("time", "status")]
  df_jags$t <- df$time
  tinits1 <- df_jags$t + max(df$time)
  is.na(tinits1) <- df_jags$status == 1
  tinits2 <- tinits1 + 5
  is.na(df_jags$t) <- df_jags$status == 0
  df_jags$is.censored <- 1 - df_jags$status
  df_jags$t.cen <- df_jags$time + df_jags$status
  data_jags <- list(N = nrow(df_jags), t.cen = df_jags$t.cen, 
                    is.censored = df_jags$is.censored, t = df_jags$t)
  data_jags$t_pred <- t_pred
  data_jags_llogis <- data_jags
  data_jags_llogis$t.log <- log(data_jags$t)
  data_jags_llogis$t.cen.log <- log(data_jags$t.cen)
  `%!in%` <- Negate(`%in%`)
  data_jags_llogis <- data_jags_llogis[names(data_jags_llogis) %!in% 
                                         c("t", "t.cen")]
  data_gomp <- list()
  data_gomp$time <- df$time
  data_gomp$status <- df$status
  data_gomp$N <- nrow(df)
  data_gomp$t_pred <- data_jags$t_pred
  n.chains = 2
  cat("Exponential Model \n")
  expo.mod <- R2jags::jags(model.file = textConnection(expo), 
                           data = data_jags, inits = inits_list("exp", n.chains), 
                           n.chains = n.chains, parameters.to.save = c("Like", 
                                                                       "lambda", "St_pred", "total_LLik"), n.iter = n.iter.jags, 
                           n.thin = n.thin.jags, n.burnin = n.burnin.jags)
  cat("Weibull Model \n")
  weib.mod <- R2jags::jags(model.file = textConnection(weibull), 
                           data = data_jags, inits = inits_list("weibull", n.chains), 
                           n.chains = n.chains, parameters.to.save = c("lambda", 
                                                                       "v", "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags, 
                           n.thin = n.thin.jags, n.burnin = n.burnin.jags)
  cat("Gamma Model \n")
  gamma.mod <- R2jags::jags(model.file = textConnection(gamma.jags), 
                            data = data_jags, inits = inits_list("gamma", n.chains), 
                            n.chains = n.chains, parameters.to.save = c("lambda", 
                                                                        "shape", "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags, 
                            n.thin = n.thin.jags, n.burnin = n.burnin.jags)
  cat("LogNormal Model \n")
  lnorm.mod <- R2jags::jags(model.file = textConnection(lnorm.jags), 
                            data = data_jags, inits = inits_list("lnorm", n.chains), 
                            n.chains = n.chains, parameters.to.save = c("mu", "sd", 
                                                                        "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags, 
                            n.thin = n.thin.jags, n.burnin = n.burnin.jags)
  cat("LogLogistic Model \n")
  llogis.mod <- R2jags::jags(model.file = textConnection(llogis.jags), 
                             data = data_jags_llogis, inits = inits_list("llogis", 
                                                                         n.chains), n.chains = n.chains, parameters.to.save = c("alpha", 
                                                                                                                                "beta", "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags, 
                             n.thin = n.thin.jags, n.burnin = n.burnin.jags)
  cat("Gompertz Model \n")
  gomp.mod <- R2jags::jags(model.file = textConnection(gompertz.jags), 
                           data = data_gomp, inits = inits_list("gompertz", n.chains), 
                           n.chains = n.chains, parameters.to.save = c("a", "b", 
                                                                       "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags, 
                           n.thin = n.thin.jags, n.burnin = n.burnin.jags)
  cat("Generalized Gamma Model \n")
  gen.gamma.mod <- R2jags::jags(model.file = textConnection(gen.gamma.jags), 
                                data = data_jags, inits = inits_list("gengamma", n.chains), 
                                n.chains = n.chains, parameters.to.save = c("r", "lambda", 
                                                                            "b", "Like", "St_pred", "total_LLik"), n.iter = n.iter.jags, 
                                n.thin = n.thin.jags, n.burnin = n.burnin.jags)
  data <- df[, c("time", "status")]
  formula <- Surv(time, status) ~ 1
  formula_temp <- stats::update(formula, paste(all.vars(formula, 
                                                        data)[1], "~", all.vars(formula, data)[2], "+."))
  mf <- tibble::as_tibble(stats::model.frame(formula_temp, 
                                             data)) %>% dplyr::rename(time = 1, event = 2) %>% dplyr::rename_if(is.factor, 
                                                                                                                .funs = ~gsub("as.factor[( )]", "", .x)) %>% dplyr::rename_if(is.factor, 
                                                                                                                                                                              .funs = ~gsub("[( )]", "", .x)) %>% dplyr::bind_cols(tibble::as_tibble(stats::model.matrix(formula_temp, 
                                                                                                                                                                                                                                                                         data)) %>% dplyr::select(contains("Intercept"))) %>% 
    dplyr::select(time, event, contains("Intercept"), everything()) %>% 
    tibble::rownames_to_column("ID")
  AIC_rps_vec <- BIC_rps_vec <- pml_vec <- waic_vec <- rep(NA, 
                                                           2)
  knots_list <- list()
  cat("Royston-Parmar Spline Model \n")
  for (i in 1:2) {
    mle.ests_rps <- flexsurv::flexsurvspline(Surv(time, 
                                                  status) ~ 1, data = df, k = i)
    init_fun_rps <- function(...) {
      list(gamma = as.numeric(mvtnorm::rmvnorm(n = 1, 
                                               mean = mle.ests_rps$res[, 1], sigma = mle.ests_rps$cov)))
    }
    k <- i
    knots <- quantile(log((mf %>% filter(event == 1))$time), 
                      seq(0, 1, length = k + 2))
    knots_list[[i]] <- knots
    B <- flexsurv::basis(knots, log(mf$time))
    DB <- flexsurv::dbasis(knots, log(mf$time))
    mm <- stats::model.matrix(formula, data)[, -1]
    if (length(mm) < 1) {
      mm <- matrix(rep(0, nrow(mf)), nrow = nrow(mf), 
                   ncol = 2)
    }
    if (is.null(dim(mm))) {
      mm <- cbind(mm, rep(0, length(mm)))
    }
    data.stan <- list(t = mf$time, d = mf$event, n = nrow(mf), 
                      M = k, X = mm, H = ncol(mm), B = B, DB = DB, mu_gamma = rep(0, 
                                                                                  k + 2), sigma_gamma = rep(5, k + 2), knots = knots)
    data.stan$mu_beta = rep(0, data.stan$H)
    data.stan$sigma_beta = rep(20, data.stan$H)
    assign(paste0("rps.", i), rstan::sampling(rps.stan_mod, 
                                              data.stan, chains = n.chains, iter = n.iter.jags, 
                                              warmup = n.burnin.jags, thin = 1, init = init_fun_rps))
    temp_gamma <- rstan::extract(get(paste0("rps.", i)), 
                                 pars = "gamma")[["gamma"]]
    LL_rps <- apply(temp_gamma, 1, function(x) {
      dsurvspline(x = df$time, gamma = x, knots = knots, 
                  log = T) * df$status + psurvspline(q = df$time, 
                                                     gamma = x, knots = knots, lower.tail = FALSE, 
                                                     log.p = T) * (1 - df$status)
    })
    LL_rps <- t(LL_rps)
    waic_vec[i] <- waic(LL_rps)[["estimates"]][3, 1]
    pml_vec[i] <- -2 * sum(log(nrow(LL_rps)/colSums(1/exp(LL_rps))))
    lp_vec <- rstan::extract(get(paste0("rps.", i)), pars = "lp__")[["lp__"]]
    LL_max_rps <- mean(lp_vec) + var(lp_vec)
    if (i == 1) {
      parm_rps <- 3
    }
    else {
      parm_rps <- 4
    }
    BIC_rps_vec[i] <- -2 * LL_max_rps + parm_rps * log(sum(df$status))
    AIC_rps_vec[i] <- -2 * LL_max_rps + parm_rps * 2
  }
  if (waic_vec[1] <= waic_vec[2]) {
    BIC_rps <- BIC_rps_vec[1]
    AIC_rps <- AIC_rps_vec[1]
    waic_rps <- waic_vec[1]
    pml_rps <- pml_vec[1]
    rps.mod <- rps.1
    knot_used <- knots_list[[1]]
    knot_num <- 1
  }
  else {
    BIC_rps <- BIC_rps_vec[2]
    AIC_rps <- AIC_rps_vec[2]
    rps.mod <- rps.2
    waic_rps <- waic_vec[2]
    pml_rps <- pml_vec[2]
    knot_used <- knots_list[[2]]
    knot_num <- 2
  }
  gamma_rps <- rstan::extract(rps.mod, pars = "gamma")[["gamma"]]
  psa_rps <- apply(gamma_rps, 1, function(x) {
    psurvspline(q = t_pred, gamma = x, knots = knot_used, 
                lower.tail = FALSE)
  })
  Surv.rps <- data.frame(time = t_pred, St_rps = rowMeans(psa_rps))
  jags.models = list(expo.mod, weib.mod, gamma.mod, lnorm.mod, 
                     llogis.mod, gomp.mod, gen.gamma.mod, rps.mod)
  AIC_vec <- BIC_vec <- rep(NA, 8)
  AIC_vec[8] <- AIC_rps
  BIC_vec[8] <- BIC_rps
  mod.names <- c("expo", "weib", "gamma", "lnorm", "llogis", 
                 "gomp", "gen.gamma")
  num.param <- c(1, 2, 2, 2, 2, 2, 3)
  PML_calc <- function(jags.mod) {
    Like_vec <- jags.mod$BUGSoutput$sims.matrix[, grep("Like", 
                                                       colnames(jags.mod$BUGSoutput$sims.matrix))]
    return(as.numeric(nrow(1/Like_vec)/colSums(1/Like_vec)))
  }
  for (i in 1:length(num.param)) {
    index <- grep("total_LLik", rownames(jags.models[[i]][["BUGSoutput"]][["summary"]]))
    LL_max <- jags.models[[i]][["BUGSoutput"]][["summary"]][index, 
                                                            1] + (jags.models[[i]][["BUGSoutput"]][["summary"]][index, 
                                                                                                                2])^2
    AIC_vec[i] <- -2 * LL_max + 2 * num.param[i]
    BIC_vec[i] <- -2 * LL_max + num.param[i] * log(sum(df$status))
    mod.temp <- get(paste0(mod.names[i], ".mod"))
    PML.temp <- assign(paste0("PML.", mod.names[i]), PML_calc(mod.temp))
    assign(paste0("PML.", mod.names[i], ".trans"), sum(log(PML.temp)) * 
             (-2))
    assign(paste0("Like.sims.", mod.names[i]), mod.temp$BUGSoutput[["sims.matrix"]][, 
                                                                                    grep("Like", colnames(mod.temp$BUGSoutput[["sims.matrix"]]))])
    if (inc_waic == F) {
      assign(paste0("WAIC.", mod.names[i], ".trans"), 
             sum(log(PML.temp)) * (-2))
    }
    else {
      assign(paste0("WAIC.", mod.names[i]), waic(log(get(paste0("Like.sims.", 
                                                                mod.names[i])))))
    }
    assign(paste0("Surv.", mod.names[i]), mod.temp$BUGSoutput[["summary"]][grep("St_pred", 
                                                                                rownames(mod.temp$BUGSoutput[["summary"]])), 1])
  }
  model.fit <- data.frame(Model = c("Exponential", "Weibull", 
                                    "Gamma", "Log-Normal", "Log-Logistic", "Gompertz", "Generalized Gamma", 
                                    paste0("Royston-Parmar ", knot_num, " knot")), minustwo_logPML = c(PML.expo.trans, 
                                                                                                       PML.weib.trans, PML.gamma.trans, PML.lnorm.trans, PML.llogis.trans, 
                                                                                                       PML.gomp.trans, PML.gen.gamma.trans, pml_rps), WAIC = c(WAIC.expo$estimates[3, 
                                                                                                                                                                                   1], WAIC.weib$estimates[3, 1], WAIC.gamma$estimates[3, 
                                                                                                                                                                                                                                       1], WAIC.lnorm$estimates[3, 1], WAIC.llogis$estimates[3, 
                                                                                                                                                                                                                                                                                             1], WAIC.gomp$estimates[3, 1], WAIC.gen.gamma$estimates[3, 
                                                                                                                                                                                                                                                                                                                                                     1], waic_rps), AIC = AIC_vec, BIC = BIC_vec)
  jags_output <- list(model.fit = model.fit, jags.models = list(expo.mod, 
                                                                weib.mod, gamma.mod, lnorm.mod, llogis.mod, gomp.mod, 
                                                                gen.gamma.mod, rps.mod), jags.surv = list(Surv.expo = Surv.expo, 
                                                                                                          Surv.weib = Surv.weib, Surv.gamma = Surv.gamma, Surv.lnorm = Surv.lnorm, 
                                                                                                          Surv.llogis = Surv.llogis, Surv.lnorm = Surv.lnorm, 
                                                                                                          Surv.gomp = Surv.gomp, Surv.gen.gamma = Surv.gen.gamma, 
                                                                                                          Surv.rps = Surv.rps))
  return(jags_output)
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
compare.surv.mods <- function (object, max_predict = 10, chng.num = "all", plot.best = 3, 
                               n.iter.jags = 2000, n.thin.jags = NULL, n.burnin.jags = NULL, 
                               gof = "WAIC", inc_waic = TRUE, km_risk = 0.1, gmp_haz_df = NULL, 
                               gpm_post_data = TRUE,  col_km = "black"){ 
  df <- object$df
  if (!is.null(gmp_haz_df)) {
    gmp_haz_df[nrow(gmp_haz_df) + 1, ] <- 0
    if (max(gmp_haz_df$time) < max_predict) {
      stop("You are predicting survival beyond the time that you have provided general population mortaility")
    }
    gmp_haz_df <- gmp_haz_df %>% arrange(time) %>% filter(time <= 
                                                            max_predict)
    if (gpm_post_data) {
      gmp_haz_df[which(gmp_haz_df$time <= max(df$time)), 
                 "hazard"] <- 0
    }
    gmp_haz_df$Cum_Haz_gmp <- cumsum(gmp_haz_df$hazard)
    t_pred <- gmp_haz_df$time
  }
  else {
    t_pred <- seq(0, max_predict, length.out = 100)
  }
  cat("Evaluating Individual log-likelihood for changepoint model \n ... can take several minutes")
  log.lik.piece <- get.loglik(object$df, object$lambda, object$changepoint)
  num.chng <- apply(object$k.stacked, 1, function(x) {
    length(na.omit(x))
  })
  model_most_prob <- as.numeric(names(which.max(object$prob.changepoint)))
  num.param_piece <- mean(model_most_prob * 2 + 1)
  log.lik.piece_most_prob <- log.lik.piece[which(num.chng == 
                                                   model_most_prob), ]
  if (chng.num != "all") {
    log.lik.piece <- log.lik.piece[which(num.chng == chng.num), 
    ]
  }
  indiv.lik.piece <- exp(log.lik.piece)
  PML.indiv.piece <- 1/indiv.lik.piece
  PML.piece <- nrow(PML.indiv.piece)/colSums(PML.indiv.piece)
  minus2logPML.piece <- -2 * sum(log(PML.piece))
  WAIC.piece <- waic(log.lik.piece)$estimate[3, 1]
  LL_max_piece <- mean(rowSums(log.lik.piece_most_prob)) + 
    var(rowSums(log.lik.piece_most_prob))
  AIC_piece <- -2 * LL_max_piece + 2 * num.param_piece
  BIC_piece <- -2 * LL_max_piece + num.param_piece * log(sum(df$status))
  jags_output <- fit_surv_models(df, max_predict = max_predict, 
                                 n.iter.jags, n.thin.jags, n.burnin.jags, gof = gof, 
                                 inc_waic = inc_waic, t_pred = t_pred)
  piecewise.mod.fit <- data.frame(Model = "Piecewise Exponential", 
                                  minustwo_logPML = minus2logPML.piece, WAIC = WAIC.piece, 
                                  AIC = AIC_piece, BIC = BIC_piece)
  mod.comp <- rbind(jags_output$model.fit, piecewise.mod.fit)
  colnames(mod.comp) <- c("Model", "-2log(PML)", "WAIC", "AIC", 
                          "BIC")
  plot_surv <- plot.changepoint(object, add.post = F, chng.num = chng.num, 
                                max_predict = max_predict, t_pred = t_pred, km_risk = km_risk,
                                col_km = col_km)
  df_surv_expo <- data.frame(Surv = jags_output$jags.surv$Surv.expo, 
                             t_pred)
  df_surv_weib <- data.frame(Surv = jags_output$jags.surv$Surv.weib, 
                             t_pred)
  df_surv_gamma <- data.frame(Surv = jags_output$jags.surv$Surv.gamma, 
                              t_pred)
  df_surv_llogis <- data.frame(Surv = jags_output$jags.surv$Surv.llogis, 
                               t_pred)
  df_surv_lnorm <- data.frame(Surv = jags_output$jags.surv$Surv.lnorm, 
                              t_pred)
  df_surv_gomp <- data.frame(Surv = jags_output$jags.surv$Surv.gomp, 
                             t_pred)
  df_surv_gen.gamma <- data.frame(Surv = jags_output$jags.surv$Surv.gen.gamma, 
                                  t_pred)
  df_surv_rps <- data.frame(Surv = jags_output$jags.surv$Surv.rps)
  colnames(df_surv_rps) <- c("t_pred", "Surv")
  df_surv_vec <- c("df_surv_expo", "df_surv_weib", "df_surv_gamma", 
                   "df_surv_llogis", "df_surv_lnorm", "df_surv_gomp", "df_surv_gen.gamma", 
                   "df_surv_rps")
  if (!is.null(gmp_haz_df)) {
    plot_surv[["layers"]][[1]]$data <- plot_surv[["layers"]][[1]]$data %>% 
      mutate(Cum_Haz_surv = -log(survival)) %>% left_join(gmp_haz_df, 
                                                          by = "time") %>% mutate(survival = exp(-(Cum_Haz_surv + 
                                                                                                     Cum_Haz_gmp))) %>% select(survival, time)
    for (q in 1:length(df_surv_vec)) {
      df_temp <- get(df_surv_vec[q])
      df_temp <- df_temp %>% left_join(gmp_haz_df, by = c(t_pred = "time")) %>% 
        mutate(Cum_Haz_surv = -log(Surv), Surv = exp(-(Cum_Haz_surv + 
                                                         Cum_Haz_gmp))) %>% select(Surv, t_pred)
      assign(df_surv_vec[q], df_temp)
    }
  }
  mu_surv_list <- list()
  for (q in 1:length(df_surv_vec)) {
    mu_surv_list[[q]] <- get(df_surv_vec[q])
  }
  mu_surv_list[[length(df_surv_vec) + 1]] <- plot_surv[["layers"]][[1]]$data
  df_order <- c("df_surv_expo", "df_surv_weib", "df_surv_gamma", 
                "df_surv_lnorm", "df_surv_llogis", "df_surv_gomp", "df_surv_gen.gamma", 
                "df_surv_rps")
  df_selc <- df_order[order(jags_output$model.fit[, gof])[1:plot.best]]
  model_names <- gsub("df_surv_", "", df_order)[order(jags_output$model.fit[, 
                                                                            gof])[1:plot.best]]
  col_vec <- c("black", "purple", "yellow", "cyan", "brown", 
               "blue", "pink", "green", "orange", "red")
  col_name <- c("KM curve", "Piecewise Expo", "Exponential", 
                "Weibull", "Gamma", "Log-Logistic", "Log-Normal", "Gompertz", 
                "Gen. Gamma", "Royston-Parmar")
  colors <- col_vec
  names(colors) <- col_name
  col_selc <- c(1, 2, order(jags_output$model.fit[, gof])[1:plot.best] + 
                  2)
  for (i in 1:plot.best) {
    plot_surv <- plot_surv + geom_line(get(df_selc[i]), 
                                       mapping = aes_string(y = "Surv", x = "t_pred", colour = paste0("'", 
                                                                                                      col_name[col_selc[i + 2]], "'")), inherit.aes = F)
  }
  if (!is.null(km_risk)) {
    result.km <- survfit(Surv(time, status) ~ 1, data = df)
    max_time <- result.km$time[max(which(result.km$n.risk/result.km$n >= 
                                           km_risk))]
  }
  else {
    max_time <- max(df$time)
  }
  plot_Surv_all <- plot_surv + scale_color_manual(values = colors[col_selc]) + 
    labs(x = "Time", y = "Survival", color = "Survival Functions") + 
    theme_classic() + theme(legend.position = c(0.9, 0.8)) + 
    geom_vline(xintercept = max_time, linetype = "dotted")
  return(list(mod.comp = mod.comp, jag.models = jags_output$jags.models, 
              jags.surv = jags_output$jags.surv, plot_Surv_all = plot_Surv_all, 
              mu_surv_list = mu_surv_list))
}


plot.pos_changepoint <- function(obj, breaks = NULL, probs = c(0.025, 0.975)){
  
  num.breaks <- as.numeric(names(which.max(obj$prob.changepoint)))
  num.changepoints <- apply(obj$k.stacked, 1, function(x){length(na.omit(x))})
  index_selc <- which(num.changepoints == num.breaks)
  
  if(num.breaks == 1){
    change.point_df <- data.frame(changetime = obj$changepoint[index_selc,1],
                                  changepoint = 1)
  }else{
    change.point_df <- data.frame(changetime = c(unlist(obj$changepoint[index_selc,1:num.breaks])),
                                  changepoint = rep(1:num.breaks,each = length(index_selc)))
  }
  
  print(change.point_df %>% filter(changepoint ==num.breaks) %>% pull(changetime)%>% quantile(probs = probs) %>% round(digits = 2))
  
  
  change.point_df$changepoint <- factor(change.point_df$changepoint)
  
  change.point_plot <- change.point_df %>% group_by(changepoint, changetime) %>%
    dplyr::summarize(n = dplyr::n()) %>% mutate(perc = (n*100/length(index_selc)))
  
  if(is.null(breaks)){
    ggplot(change.point_plot %>% filter(perc > .5), aes(x = changetime, y = perc, color=changepoint))+
      geom_pointrange(aes(ymin=0, ymax=perc), size = 0.02)+
      scale_y_continuous(name="Probability of Change-point (%)")+
      ggtitle("Posterior Distribution of Change-points")+
      scale_x_continuous(name="Time", breaks = round(seq(round(min(change.point_plot$changetime),1),
                                                         max(change.point_plot[which(change.point_plot$perc >0),
                                                                               "changetime"]), by = 0.5),1) )
  }else{
    ggplot(change.point_plot %>% filter(perc > .5), aes(x = changetime, y = perc, color=changepoint))+
      geom_pointrange(aes(ymin=0, ymax=perc), size = 0.02)+
      scale_y_continuous(name="Probability of Change-point (%)")+
      ggtitle("Posterior Distribution of Change-points")+
      scale_x_continuous(name="Time", breaks = breaks )
    
  }
  
  
  
  
}


compare_boot_sims <- function (mod_parametric_orig, follow_up_data) {
  
  t <- mod_parametric_orig$mu_surv_list[[1]][,2]
  mod_names <- c("expo", "weibull", "gamma", "llogis", "lnorm", 
                 "gomp", "gen.gamma", "rps", "piecewise")
  
  
  for (i in 1:length(mod_names)) {
    assign(paste0("Surv.", mod_names[i]), mod_parametric_orig$mu_surv_list[[i]][,1])
    
  }
  
  n.boots <- 1000
  bs <- sjstats::bootstrap(follow_up_data, n.boots)
  AUC_diff <- AUC_diff2 <- matrix(nrow = n.boots, ncol = 10)
  AUC_true <- rep(NA, n.boots)
  surv_km <- survfit(Surv(time, status) ~ 1, data = follow_up_data)
  for (b in 1:n.boots) {
    index <- b
    df_surv_boot <- bs[[1]][[index]][["data"]][bs[[1]][[index]]$id, 
    ]
    surv_boot_km <- survfit(Surv(time, status) ~ 1, data = df_surv_boot)
    surv_boot_km.df <- data.frame(cbind(surv_boot_km[[c("time")]], 
                                        surv_boot_km[[c("surv")]], surv_boot_km[[c("upper")]], 
                                        surv_boot_km[[c("lower")]]))
    colnames(surv_boot_km.df) <- c("time", "survival", "upper", 
                                   "lower")
    surv_km_time <- rep(NA, length(t))
    for (i in 1:length(t)) {
      if (t[i] < surv_boot_km.df$time[1]) {
        surv_km_time[i] <- 1
      }
      else if (t[i] > surv_boot_km.df$time[nrow(surv_boot_km.df)]) {
        surv_km_time[i] <- NA
      }
      else {
        surv_km_time[i] <- surv_boot_km.df$survival[max(which(surv_boot_km.df$time <= 
                                                                t[i]))]
      }
    }
    t_eval <- !is.na(surv_km_time)
    AUC_true[b] <- integrate.xy(t[t_eval], surv_km_time[t_eval])
    AUC_diff[b, ] = abs(c(integrate.xy(t[t_eval], Surv.expo[t_eval]), 
                          integrate.xy(t[t_eval], Surv.weibull[t_eval]),
                          integrate.xy(t[t_eval], Surv.gamma[t_eval]),
                          integrate.xy(t[t_eval],Surv.lnorm[t_eval]), 
                          integrate.xy(t[t_eval],Surv.llogis[t_eval]),
                          integrate.xy(t[t_eval],Surv.gomp[t_eval]), 
                          integrate.xy(t[t_eval], Surv.gen.gamma[t_eval]),
                          integrate.xy(t[t_eval], Surv.rps[t_eval]),
                          integrate.xy(t[t_eval], Surv.piecewise[t_eval]), 
                          NA) - AUC_true[b])
    Surv.piecewise
    AUC_diff2[b, ] = c(integrate.xy(t[t_eval], abs(Surv.expo[t_eval] - surv_km_time[t_eval])),
                       integrate.xy(t[t_eval], abs(Surv.weibull[t_eval] - surv_km_time[t_eval])), 
                       integrate.xy(t[t_eval], abs(Surv.gamma[t_eval] - surv_km_time[t_eval])),
                       integrate.xy(t[t_eval], abs(Surv.lnorm[t_eval] - surv_km_time[t_eval])), 
                       integrate.xy(t[t_eval], abs(Surv.llogis[t_eval] - surv_km_time[t_eval])),
                       integrate.xy(t[t_eval], abs(Surv.gomp[t_eval] - surv_km_time[t_eval])),
                       integrate.xy(t[t_eval], abs(Surv.gen.gamma[t_eval] - surv_km_time[t_eval])),
                       integrate.xy(t[t_eval], abs(Surv.rps[t_eval] - surv_km_time[t_eval])), 
                       integrate.xy(t[t_eval], abs(Surv.piecewise[t_eval] - surv_km_time[t_eval])), NA)
  }
  index_selc <- which(t <= max(follow_up_data$time))
  surv_KM <- survival:::survmean(surv_km, rmean = max(surv_km$time))
  rmean_KM <- as.numeric(surv_KM[[1]][min(grep("rmean", names(surv_KM[[1]])))])
  AUC = c(integrate.xy(t[index_selc], Surv.expo[index_selc]), 
          integrate.xy(t[index_selc], Surv.weibull[index_selc]), 
          integrate.xy(t[index_selc], Surv.gamma[index_selc]), 
          integrate.xy(t[index_selc], Surv.lnorm[index_selc]), 
          integrate.xy(t[index_selc], Surv.llogis[index_selc]), 
          integrate.xy(t[index_selc], Surv.gomp[index_selc]), 
          integrate.xy(t[index_selc], Surv.gen.gamma[index_selc]), 
          integrate.xy(t[index_selc], Surv.rps[index_selc]), integrate.xy(t[index_selc], 
                                                                          Surv.piecewise[index_selc]), rmean_KM)
  model_names_full <- c("Exponential", "Weibull", "Gamma", 
                        "Log-Normal", "Log-Logistic", "Gompertz", "Generalized Gamma", 
                        "Royston-Parmar", "Piecewise")
  index_vals <- rep(NA, length(model_names_full))
  for (i in 1:length(model_names_full)) {
    index_vals[i] <- grep(paste0("^", model_names_full[i]), 
                          mod_parametric_orig$mod.comp$Model)
  }
  mod_compfinal <- rbind(mod_parametric_orig$mod.comp[index_vals, 
                                                      c("Model", "WAIC")], c(NA, NA, NA))
  mod_compfinal$Model[nrow(mod_compfinal)] <- "True Observations"
  mod_compfinal <- cbind(mod_compfinal, AUC, AUC_diff = colMeans(AUC_diff))
  mod_compfinal <- mod_compfinal[order(mod_compfinal$AUC_diff),] %>% mutate_if(is.numeric, round, digits = 2)
  mod_compfinal
}
