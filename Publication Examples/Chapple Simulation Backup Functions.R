GetSurvPEH = function(x, s, lam, J) {
  y = x
  for (m in 1:length(x)) {
    for (k in 1:(J + 1)) {
      if ((x[m] > s[k]) && (x[m] <= s[k + 1])) {
        if (k > 1) {
          y[m] = exp(-exp(lam[k]) * (y[m] - s[k]) -
                       sum(exp(lam[1:(k - 1)]) * (s[2:k] - s[1:(k -
                                                                  1)])))
        }
        else {
          y[m] = exp(-exp(lam[k]) * (y[m] - s[k]))
        }
      }
    }
  }
  return(y)
}



GetSurvPEH_collapse <- function (x, mod) {

  # sample.index <- sample(nrow(mod$changepoint), size = 1000)
  # mod$k.stacked <-  mod$k.stacked[sample.index,]
  # mod$changepoint <-  mod$changepoint[sample.index,]
  # mod$lambda <-  mod$lambda[sample.index,]


  SurvHold = matrix(nrow = nrow(mod$changepoint), ncol = length(x))
  j.val <- apply(mod$k.stacked,1, function(x){length(na.omit(x))})
  n_rows <- nrow(mod$changepoint)
  log_lambda <- log(mod$lambda)
  for (b in 1:n_rows) {
    s = na.omit(c(0,mod$changepoint[b, ],100))
    lam =  na.omit(log_lambda[b, ])
    J = j.val[b]
    SurvHold[b, ] = GetSurvPEH(x, s, lam, J)
  }
  return(SurvHold)
}




GetCensorTimePIECE=function(PCT, ##Desired Percentage of
                            n, ##Number of patients
                            Rate_ACC, ##Accrual Per unit time
                            Rate_piece,
                            chng,
                            B, #Bootstraps,
                            factor_pool = 10 #Assume accruals happen together
){

  if(PCT==1){

    return(10000)

  }else{

    Storage =rep(NA,B)

    for(b in 1:B){
      ##First Let's Simulate accural times


      ACC=cumsum(rexp(n/factor_pool,Rate_ACC))
      ACC <-sample(ACC, n, replace = T)

      ##Now get the trial Times
      TIMES = PiecewiseChangepoint:::rpwexp(n, lam=Rate_piece, s = chng)

      ##Trial Times at death
      Y=TIMES+ACC
      Storage[b]=quantile(Y,prob=PCT)

    }
    ##Ok Now we have the Storage
    ##Lets find the Cut Trial time cut point to get desired censoring.

    return(mean(Storage))
  }


}

