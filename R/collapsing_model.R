

#' Fitting Bayesian Piecewise Exponential Models
#'
#' Implementation of piecewise exponential models typically assume that the number of a location of change-points is known. In many cases both of these quantities are unknown. This function estimates the posterior probability of change-point numbers and locations.
#'
#' @param df A dataframe with survival times and status indicator.
#' @param n.iter number of total iterations per chain excluding burn in (default 10750).
#' @param burn_in length of burn in, i.e. number of iterations to discard at the beginning of each chain (default 750).
#' @param n.chains number of Markov chains (default 3).
#' @param max.num.breaks maximum number of allowed change-points (default 6).
#' @param n.thin thinning rate. Must be a positive integer. Set n.thin > 1 to save memory and computation time if n.iter is large (default 1).
#' @param lambda.prior rate parameter for the Poisson prior for number of change-points (default 1).
#' @param timescale Scale of time. Can be "years", "months" or "days". Defines the values of the alpha and beta parameters so that the analysis on different timescales give the same posterior model probabilities.  
#' @param seed.val set seed value to get the same resuls
#' @return A list of the class "changepoint" with the following items:
#'  \itemize{
#'   \item \strong{k.stacked}: Matrix of change-point indices. Number refers to the number of events.
#'   \item \strong{changepoint}: Matrix of change-point indices in terms of time (rather than number of events).
#'   \item \strong{lambda}: Matrix with post-hoc estimation of the hazards.
#'   \item \strong{prob.changepoint}: Posterior probabilities of changepoint models.
#'   \item \strong{k}: k.stacked but in array form with each slice indicating a chain.
#'   \item \strong{beta.array} Posterior samples of the rate of the gamma prior for the (\eqn{\lambda}). Note that this is the rate of the (lower level) prior and not the hyperprior.
#' }
#'
#'
#' @importFrom survival survfit Surv
#' @importFrom magrittr %>%
#' @import ggplot2
#' @export
#' @md
#' @examples \dontrun{Collapsing_Model <- collapsing.model(df,
#'                                                         n.iter = 5000,
#'                                                         burn_in = 750,
#'                                                         n.chains = 2,
#'                                                         timescale = "months",
#'                                                         beta.hyper1 = 1,
#'                                                         beta.hyper2 = 1)}

collapsing.model <- function(df,
                             n.iter = 5750,
                             burn_in = 750,
                             n.chains = 2,
                             max.num.breaks = 6,
                             n.thin = 1,
                             lambda.prior = 1,
                             timescale = "months",
                             seed.val = NULL
) {
  
  
  df <- df[order(df$time), ]
  m <- n.iter
  df_event <- df[which(df$status == 1), ]
  df_event_unique <- unique(df_event[, c("time", "status")])
  time_diffs <- df_recast(df)
  n <- nrow(time_diffs)
  # Depreciated MLE should always be false
  if(is.numeric(seed.val)){
    set.seed(seed.val)  
  }
  
  
  MLE = FALSE
  
  # @param alpha.hyper value of shape parameter for the gamma prior for the hazard (\eqn{\lambda}). Should always be set to 1. As noted there is a gamma prior on \eqn{\lambda} but also a gamma prior (hyperprior) on the rate parameter of this prior.
  # @param beta.hyper1 first hyperprior (shape) for the rate parameter of the gamma prior for \eqn{\lambda}. Should be set to 1.
  # @param beta.hyper2 first hyperprior (rate) for the rate parameter of the gamma prior for \eqn{\lambda}. To insure that inferences are the same when different timescales (i.e. months vs years) are chosen, this hyperprior must be changed. If the analysis is in years this value should be equal to 1, if in months it should be 1/12 and if in days should be 1/365.
  

  beta.hyper1 <-  alpha.hyper <- 1
  if(timescale =="years"){
    beta.hyper2 <- 1
  }
  if(timescale =="months"){
    beta.hyper2 <- 1/12
  }
  if(timescale =="days"){
    beta.hyper2 <- 1/365
  }

  

  changepoint <- k <- array(NA, dim = c(m, max.num.breaks, n.chains))
  beta_array <- array(NA, dim = c(m, 1, n.chains))
  
  for (c in 1:n.chains) {
    k_curr <- array(NA, dim = c(m, max.num.breaks))
    
    initial.breaks <- min(rpois(1, lambda = 1), 6)
    if (initial.breaks != 0) {
      int_unorder <- sample(2:(n - 1), initial.breaks, replace = FALSE)
      int_locs <- int_unorder[order(int_unorder)]
      
      while (pos_prob(int_locs, n, MLE) == 0) {
        int_unorder <- sample(2:(n - 1), initial.breaks, replace = FALSE)
        int_locs <- int_unorder[order(int_unorder)]
      }
      
      
      k_curr[1, 1:initial.breaks] <- int_locs
    }
    core_output <- RJMCM_core(
      k_curr, max.num.breaks, time_diffs, MLE, alpha.hyper, beta.hyper1, beta.hyper2,
      lambda.prior
    )
    
    k[, , c] <- core_output[["k"]]
    beta_array[, , c] <- core_output[["beta"]]
  }
  
  samp_retain <- seq(from = n.thin, to = m, by = n.thin)
  
  if (n.chains == 1) {
    k <- k[samp_retain[which(samp_retain > burn_in)], , 1]
    beta_array <- beta_array[samp_retain[which(samp_retain > burn_in)], , 1]
  } else {
    k <- k[samp_retain[which(samp_retain > burn_in)], , ]
    beta_array <- array(beta_array[samp_retain[which(samp_retain > burn_in)], , ],
                        dim = c(length(samp_retain[which(samp_retain > burn_in)]), 1, n.chains)
    )
  }
  
  
  # Stack the k chains
  if (n.chains == 1) {
    k.stacked <- k
    beta.stacked <- beta_array
  }
  if (n.chains > 1) {
    k.stacked <- as.matrix(k[, , 1])
    beta.stacked <- as.matrix(beta_array[, , 1])
    for (i in 2:n.chains) {
      k.stacked <- rbind(k.stacked, as.matrix(k[, , i]))
      beta.stacked <- rbind(beta.stacked, as.matrix(beta_array[, , i]))
    }
  }
  
  # Calculate the lambdas
  num.changepoints <- unlist(apply(k.stacked, 1, function(x) {
    length(na.omit(x))
  }))
  lambda.df_var <- lambda.df_mu <- lambda.df <- changepoint <- array(NA, dim = c(nrow(k.stacked), max.num.breaks + 1))
  
  for (i in 1:nrow(k.stacked)) {
    if (num.changepoints[i] == 0) {
      res_matrix_temp <- exposure_death_alt(time_diffs, NA)
    } else {
      res_matrix_temp <- exposure_death_alt(time_diffs, k.stacked[i, 1:num.changepoints[i]])
    }
    # Calculate the lambda at the posterior values of the hyperpriors (alpha is fixed, beta is not)
    
    lambda.df[i, 1:(num.changepoints[i] + 1)] <- apply(res_matrix_temp, 1, function(x) {
      rgamma(1, x[1] + alpha.hyper, x[2] + beta.stacked[i])
    }) # Hyperparameters alpha and beta = 1
    lambda.df_mu[i, 1:(num.changepoints[i] + 1)] <- apply(res_matrix_temp, 1, function(x) {
      (x[1] + alpha.hyper) / (x[2] + beta.stacked[i])
    })
    lambda.df_var[i, 1:(num.changepoints[i] + 1)] <- apply(res_matrix_temp, 1, function(x) {
      (x[1] + alpha.hyper) / (x[2] + beta.stacked[i])^2
    })
    
    if (num.changepoints[i] != 0) {
      changepoint[i, 1:(num.changepoints[i])] <- df_event_unique[k.stacked[i, 1:num.changepoints[i]], "time"]
    }
  }
  
  
  prob.changepoint <- table(num.changepoints) / sum(table(num.changepoints))
  
  piecewise.obj <- list(
    k.stacked = k.stacked,
    changepoint = changepoint,
    lambda = lambda.df,
    prob.changepoint = prob.changepoint,
    #lambda.df_mu = lambda.df_mu,
    #lambda.df_var = lambda.df_var,
    k = k,
    beta.array = beta_array,
    df = df
  )
  class(piecewise.obj) <- "changepoint"
  return(piecewise.obj)
}
