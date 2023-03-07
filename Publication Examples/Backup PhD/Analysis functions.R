#Describe functions
#Hazard.plt:

#Function which takes a vector of times & status indicators along with interval widths and returns
# a plot of these hazards along with a smoothed hazard.


harzard.plt <- function(time.vec, cens.vec, Outcome = "Survival", trt = "",
                        lrg.haz.int = 1, sml.haz.int =3,  bw.method = "local"){

  #Plot the hazards
  max.time <- max(time.vec)
  sml.haz.int.inv <- sml.haz.int/12

  result.hazard.pe.lrg <- pehaz(time.vec, cens.vec, width= lrg.haz.int, max.time=max.time)

  result.smooth <- muhaz(time.vec, cens.vec, b.cor="left",
                         max.time=max.time, bw.method = bw.method)

  #Plot the Survival function

  plot(result.hazard.pe.lrg,  col="black", main = paste0("Hazards for the ", Outcome,
                                                         " of ",trt,  " arm"))
  result.hazard.pe.sml <- pehaz(time.vec, cens.vec, width=sml.haz.int.inv, max.time=max.time)
  lines(result.hazard.pe.sml,  col="blue")

  lines(result.smooth, col = "red", lty= 2)
  legend("topright", legend = c(paste0(lrg.haz.int," year hazard step function"),
                                paste0(sml.haz.int," month hazard step function"),
                                "Smoothed hazard function"),
         col = c("black", "blue", "red"), lty = c(1,2,2), cex = 0.8)



  return(result.smooth)
}



smoothed.survival <- function(result.smooth, time.vec, cens.vec,Outcome = "Survival", trt = ""){

  haz <- result.smooth$haz.est
  times <- result.smooth$est.grid
  surv <- exp(-cumsum(haz[1:(length(haz)-1)]*diff(times)))
  result.km <- survfit(Surv(time.vec, cens.vec) ~ 1,
                       conf.type="none")
  max.time <- max(time.vec)

  plot(result.km, conf.int=F, mark="|", xlab="Time",
       xlim=c(0,max.time), ylab="Survival probability",
       main = paste0("Survival for the ", Outcome,
                     " of ",trt,  " arm"))
  lines(surv ~ times[1:(length(times) - 1)], col = "blue", lty = 2)
  legend("topright", legend = c("Kaplan Meier",
                                "Survival derived from smoothed hazards"),
         col = c("black", "blue"), lty = c(1,2), cex = 0.8)

}


#From viewing the data we could suggest that patients may be "cured",
#or exposed to much lower hazards. However, for the purposes of illustrations
#I consider the piecewise exponential models

#Consider Piecewise Exponential Models

#Frequentist Approach , exponential chi sqaure distribution
#https://cran.r-project.org/web/packages/pssm/pssm.pdf
#Joint model, frailty model
#Do plots with ASAUR --- Done
#See chapter on Klein

#breaks should be based on Expected number of events


predict.piecewise <- function(piecewise.model, max.predict ){

  break.times  <- as.numeric(piecewise.model$breaks)[-1]
  inter.length <- attributes(piecewise.model$breaks)$h
  cum.time <- cumsum(inter.length)
  cum.hzd      <- piecewise.model$Lambda[1,]
  surv.break <- exp(-piecewise.model$Lambda[1,])
  hzd <- piecewise.model$lambda[1,]
  time <- seq(from = 0, to  = max.predict, by = 0.1)

  St <- rep(NA, length(time))
  t.len <- length(time)
  for (i in 1:t.len){
    if(length(which(time[i] > cum.time))== 0) {
      St[i] <-  exp(-hzd[1]*time[i])

    }else if(length(which(time[i] > cum.time))== length(hzd)){
      time.diff <- time[i] -cum.time[length(cum.time)]
      St[i] <- exp(-cum.hzd[length(cum.time)])*exp(-hzd[length(hzd)]*time.diff)

    }else{
      index <- max(which(time[i] > cum.time)) +1
      time.diff <- time[i] -cum.time[index]
      St[i] <- exp(-cum.hzd[index])*exp(-hzd[index]*time.diff)

    }

  }

  fit.plot <- survfit(Surv(FAILTIME, FAILCENS)~TRT, data = trt.df)
  fit.plot <- survfit(piecewise.model$mf[,1]~piecewise.model$mf[,2])
  plot(fit.plot, xlim = c(0,max.predict),
       xlab = "Time",
       ylab = "Survival",
       main = "Piecewise Exponential Model and Prediction")
  lines(y = St, x = time, col = "red", lty =2)
  points(y = surv.break, x = break.times, pch = 23, cex = 1, bg = "green", col = "black")
  legend( x="topright",
          legend=c("Kaplan Meier","Piecewise exponential","Breakpoints"),
          col=c("black","red","green"), lwd=1, lty=c(1,2,NA),
          pch=c(NA,NA,18),merge=FALSE)

}


#Piecewise Exponential Model function


# time = vector of times
# status = status of patient, censored or event
# breakpoints = times at which we want a breakpoint

piecewise.exponential <- function(breakpoints = NULL, time , status){

  time <- time
  status <- status
  df <- data.frame(time= time,
                   status = status)
  df <- df[order(time), ]


    log.likelihood <- -sum(df[,"status"])



  ##Note some of the code is redundant and refers to when the Nelder-Mead algorithm was applied.
  #Errors can occur if the algorithm picks a value that is outside the time horizon
  #or if a smaller time is picked for a later interval (or if the same time is picked for both!)
  if(length(breakpoints)>0){

  breaks <- c(0,breakpoints, max(time)+.001)

  #Because we use left breakpoints the last observation is NA, could avoid this
  #by adding an "epsilon" (small value to ensure it is above max value)
  df$timepoints <- cut(df$time, breaks = breaks, right = FALSE)

  #Create a dataframe with 1's and 0's for each interval
  df <- cbind(df, model.matrix(~ timepoints + 0, data=df))
  #Drop redundant columns
  df <- df[,-which(colnames(df) == "timepoints")]

  #Make sure that middle intervals at least 1 event
    if(length(breakpoints)> 1 && any(apply(df[,3:(ncol(df)-1)],2, function(x)all(x== 0)))){

      log.likelihood = NA
      return(log.likelihood)


  }else{

  #Define some placeholder vectors
  ncol.df <- ncol(df)
  piecewise.haz <- rep(NA, length(breakpoints)+1)

  nrows <- nrow(df)

  #Define some counters
  j <- 0
  t.prev <- 0
  #For 3rd column (first time interval) to the final time interval
  for(i in 3:ncol.df){

    j <- j + 1
    # browser()
    temp.df <- df[min(which(df[,i] == 1)):nrow(df),]
    nrows <- nrow(temp.df)
    end.point <- max(which(temp.df[,i] == 1))

    #If we are in an interval before the last one all observations
    # after the breakpoints are censored
    if(end.point < nrows){
      temp.df[end.point:nrows, "status"] <- 0
      temp.df[end.point:nrows, "time"] <- breakpoints[j]
    }
    # Center the time to be from the start of that interval
    temp.df$time <-  temp.df$time - t.prev

    #As defined in the chapter the log likelihood of a exponential constant hazard is
    #(Number of Events)*log(Number of Events/Exposure time in the interval)

    log.likelihood.new <- sum(temp.df[,"status"])*log(sum(temp.df[,"status"])/sum(temp.df$time))


    #If no events in the final interval then the likelihood is zero
    if(is.nan(log.likelihood.new)){
      log.likelihood.new <- 0
    }


    log.likelihood <- log.likelihood.new + log.likelihood
    #The current time is now the previous time for the next interval
    t.prev <- breakpoints[j]

  }

  return(log.likelihood)

  }
  }else{
    #If no breakpoints provided by the use just fit a one parameter exponential model
    log.likelihood <- sum(df[,"status"])*log(sum(df[,"status"])/sum(df$time)) + log.likelihood

    return(log.likelihood)

  }

}



grid.search.piecewise <- function(min.break, max.break, grid.width, num.breaks, min.break.width, time, status){

  seq.new <- seq(min.break, max.break, by = grid.width)
  combn.df <- t(combn(seq.new, num.breaks))
  combn.df.diff <- t(apply(combn.df, 1, diff))
  drop.index <- which(apply(combn.df.diff, 1, function(x)(any(x < min.break.width))) == TRUE)

  if(length(drop.index) !=0){
    combn.df <- data.frame(combn.df[-drop.index,])
  }

  combn.df <- data.frame(combn.df)
  print(paste0("Number of evaluations = ", nrow(combn.df)))
  combn.df$Likelihood <- apply(combn.df, 1, piecewise.exponential, time, status)

  #drop.index2 <- which(combn.df$Likelihood == -1.000000e+20)

  #if(length(drop.index2) !=0){
  #combn.df <- combn.df[-drop.index2,]
  #}

  if(num.breaks == 2){

    #https://www.r-bloggers.com/creating-a-matrix-from-a-long-data-frame/

    combn.df$X1_mod <- match(combn.df$X1, sort(unique(combn.df$X1)))
    combn.df$X2_mod <- match(combn.df$X2, sort(unique(combn.df$X2)))

    n1 <- length(unique(combn.df$X1))
    n2 <- length(unique(combn.df$X2))

    matrix.triangle  <- with(combn.df, {
      M <- matrix(nrow=n1, ncol=n2)
      M[cbind(X1_mod, X2_mod)] <- Likelihood
      M
    })

    #kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
    #persp(x =combn.df$X1_mod,  y =combn.df$X2_mod, z = M)
    output <- list()

    output[[1]] <- plot_ly(x = unique(combn.df$X1), y = unique(combn.df$X2), z = matrix.triangle) %>% add_surface() %>%
      layout(scene = list(xaxis=list(title = "Time of 1st Breakpoint"),
                          yaxis=list(title = "Time of 2nd Breakpoint"),
                          zaxis = list(title = "Log-Likelihood")))

    combn.df <- subset(combn.df, select = -c(X1_mod,X2_mod))
    output[[2]] <-  combn.df

    return(output)

  }

  return(combn.df)

}



### New Piecewise model
grid.search.piecewise.pchreg <- function(min.break, max.break, grid.width, num.breaks, min.break.width, time, status){

  seq.new <- seq(min.break, max.break, by = grid.width)
  combn.df <- t(combn(seq.new, num.breaks))
  combn.df <- data.frame(cbind(0,combn.df,max(time)))
  combn.df.diff <- t(apply(combn.df, 1, diff))
  drop.index <- which(apply(combn.df.diff, 1, function(x)(any(x < min.break.width))) == TRUE)

  if(length(drop.index) !=0){
    combn.df <- data.frame(combn.df[-drop.index,])
  }

  print(paste0("Number of evaluations = ", nrow(combn.df)))

  Surv.formula <- as.formula("Surv(time, status) ~ 1")
  #https://www.r-bloggers.com/skip-errors-in-r-loops-by-not-writing-loops/
  possibly_some_function = possibly(function(x)pchreg(Surv.formula,breaks = x)$logLik,
                                    otherwise = NA)
  combn.df$Likelihood <- apply(combn.df, 1, possibly_some_function)

  #combn.df$Likelihood <- apply(combn.df, 1, function(x)pchreg(Surv.formula,breaks = x)$logLik)

  #drop.index2 <- which(combn.df$Likelihood == -1.000000e+20)

  #if(length(drop.index2) !=0){
  #combn.df <- combn.df[-drop.index2,]
  #}

  if(num.breaks == 2){

    #https://www.r-bloggers.com/creating-a-matrix-from-a-long-data-frame/

    combn.df$X2_mod <- match(combn.df$X2, sort(unique(combn.df$X2)))
    combn.df$X3_mod <- match(combn.df$X3, sort(unique(combn.df$X3)))

    n1 <- length(unique(combn.df$X2))
    n2 <- length(unique(combn.df$X3))

    matrix.triangle  <- with(combn.df, {
      M <- matrix(nrow=n1, ncol=n2)
      M[cbind(X2_mod, X3_mod)] <- Likelihood
      M
    })

    #kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
    #persp(x =combn.df$X1_mod,  y =combn.df$X2_mod, z = M)
    output <- list()

    output[[1]] <- plotly::plot_ly(x = unique(combn.df$X2), y = unique(combn.df$X3), z = matrix.triangle) %>% add_surface() %>%
      layout(scene = list(xaxis=list(title = "Time of 2nd Breakpoint"),
                          yaxis=list(title = "Time of 1st Breakpoint"),
                          zaxis = list(title = "Log-Likelihood")))

    combn.df <- subset(combn.df, select = -c(X2_mod,X3_mod))
    output[[2]] <-  combn.df

    return(output)

  }

  return(combn.df)

}

##





nelson.aalen.haz <- function(time,status){
  #math.ucsd.edu/~rxu/math284/slect2.pdf
  Surv <-survfit(Surv(time, status) ~ 1,type = "fleming-harrington")$surv
  cum_haz <- -log(Surv)
  time.orderd <- survfit(Surv(time, status) ~ 1,type = "fleming-harrington")$time
  idx <- match(time,time.orderd)
  return(cum_haz[idx])
}


least.squares.method.mult.breaks <- function(breakpoints, df){
  #browser()
  breaks <- c(0,breakpoints, max(df$time))
  n.haz <- length(breakpoints)+1
  start.searching.here <- rbeta(n.haz,1,1)
  #start.searching.here <- rep(-1,n.haz)
  p <- start.searching.here

  #Because we use left breakpoints the last observation is NA, could avoid this
  #by adding an "epsilon" (small value to ensure it is above max value)

  df$timepoints <- cut(df$time, breaks = breaks, right = FALSE)
  df[nrow(df),"timepoints"] <- df[nrow(df)-1,"timepoints"]

  #Create a dataframe with 1's and 0's for each interval
  df <- cbind(df, model.matrix(~ timepoints + 0, data=df))
  #Drop redundant columns
  df <- df[,-which(colnames(df) == "timepoints")]
  #browser()
  break_2 <- c(0,breakpoints[-length(breakpoints)])
  break_diff <- breakpoints - break_2

  time_start <- min(grep("timepoints",colnames(df)))
  time_end <- max(grep("timepoints",colnames(df)))
  #df$diff <- rep(c(break_diff,0), colSums(df[time_start:time_end]))


  least.square.multi.func <-function(df, p, breakpoints, break_diff){
    #browser()
    if(length(which(p<0) != 0)){
      sum_total <- 1E10
      return(sum_total)
    }else{
      sum1    <- ((df$avg_haz-p[1])^2)%*%df[,(time_start)]
      sum2     <- ((df$avg_haz-p[2]-(p[1]-p[2])*(breakpoints[1]/df$time))^2)%*%df[,(time_start+1)]

      #sum1    <- ((df$avg_haz-p[1])^2)%*%df[,time_start]
      #part_2  <- p[1]*(breakpoints[1]/df$time)
      #sum2    <- ((y-p[2] + p[2]*(breakpoints[1]/df$time) - part_2)^2)%*%df[,(time_start+1)]
      breakpoint_par <- 1
      sum_total <- sum(sum1,sum2)

      if(length(breakpoints>1)){
        for(i in (time_start+2):time_end){

          prev_param <- i-time_start
          prev_beta <- p[1:prev_param]%*%break_diff[1:prev_param]
          #final part of equation
          prev_beta_time <- prev_beta/df$time

          breakpoint_par <- breakpoint_par + 1
          current_breakpoint <- breakpoints[breakpoint_par]

          current_par <- prev_param + 1
          sum_current<- ((df$avg_haz - p[current_par]+ p[current_par]*(current_breakpoint/df$time)-prev_beta_time)^2)%*%df[,i]

          sum_total <- sum_total + sum_current

        }

      }

      return(sum_total)

    }
  }

  return(optim(par = start.searching.here,
               fn = least.square.multi.func,
               df = df, breakpoints = breakpoints, break_diff = break_diff ))

}




grid.search.hazard_multi <- function(time,
                                     status,
                                     min.break,
                                     max.break,
                                     num.breaks,
                                     min.break.width,
                                     grid.width){

  Cum_haz  <- nelson.aalen.haz(time, status)
  df <- data.frame(time = time, status = status, Cum_haz = Cum_haz)
  df$avg_haz <- Cum_haz/time
  df <- df[-which(df$time == 0),]
  df <- df[order(df$time),]
  seq.new <- seq(min.break, max.break, by = grid.width)

  if(num.breaks == 1){
    combn.df <- data.frame(seq.new)
  }else{
    combn.df <- t(combn(seq.new, num.breaks))
    combn.df.diff <- t(apply(combn.df, 1, diff))
    drop.index <- which(apply(combn.df.diff, 1, function(x)(any(x < min.break.width))) == TRUE)
    if(length(drop.index) !=0){
      combn.df <- data.frame(combn.df[-drop.index,])
    }

  }


  Likelihood <- apply(combn.df, 1, least.squares.method.mult.breaks, df)
  least.sqrs <- c(NA, rep(nrow(combn.df)))
  for(i in 1:nrow(combn.df)){
    #Likelihood[[i]]<-  append(Likelihood[[i]], as.numeric(combn.df[i,]))
    Likelihood[[i]][["breakpoints"]] <-  as.numeric(combn.df[i,])
    least.sqrs[i] <- Likelihood[[i]][["value"]]
  }


  Likelihood["min_square"] <- which(least.sqrs==min(least.sqrs))

  if(num.breaks == 1){
    plot(y = unlist(Likelihood)[which(names(unlist(Likelihood)) == "value")],
         x = combn.df[,"seq.new"],
         ylab = "Least Squares",
         xlab = "Time of breakpoint")}

  return(Likelihood)

}


###Redundant function

piecewise.exponential.model <- function(time, status, breakpoints = numeric(0)){

  #stick the two vectors together to get a dataframe
  df <- data.frame(time= time,
                   status = status)


  # If user has specified breakpoints
  if(length(breakpoints)){

    #Because we use left breakpoints the last observation is NA, could avoid this
    #by adding an "epsilon" (small value to ensure it is above max value)
    df$timepoints <- cut(df$time, breaks=c(0,breakpoints, max(df$time)), right = FALSE)
    df[nrow(df),"timepoints"] <- df[nrow(df)-1,"timepoints"]

    #Create a dataframe with 1's and 0's for each interval
    df <- cbind(df, model.matrix(~ timepoints + 0, data=df))
    #Drop redundant columns
    df <- df[,-which(colnames(df) == "timepoints")]

    #Define some placeholder vectors
    ncol.df <- ncol(df)
    piecewise.haz <- rep(NA, length(breakpoints)+1)
    log.likelihood <- rep(NA, length(breakpoints)+1)
    nrows <- nrow(df)

    #Define some counters
    j <- 0
    t.prev <- 0
    #For 3rd column (first time interval) to the final time interval
    for(i in 3:ncol.df){
      j <- j + 1

      temp.df <- df[min(which(df[,i] == 1)):nrow(df),]
      nrows <- nrow(temp.df)
      end.point <- max(which(temp.df[,i] == 1))

      #If we are in an interval before the last one all observations
      # after the breakpoints are censored
      if(end.point < nrows){
        temp.df[end.point:nrows, "status"] <- 0
        temp.df[end.point:nrows, "time"] <- breakpoints[j]
      }
      # Center the time to be from the start of that interval
      temp.df$time <-  temp.df$time - t.prev

      # If all obserations are censored, MLE of hazard is 0 and Loglikelihood is 0
      if(all(temp.df[,"status"] ==0)){
        piecewise.haz[j] <- 0
        log.likelihood[j] <- 0
      }else{
        #Else run the optimization algorithm
        result <- optim(par=1, fn=jll.expo, method= "L-BFGS-B",
                        lower=c(0.000001), upper=c(100),
                        control=list(fnscale = -1),
                        tt=temp.df$time, status=temp.df[,"status"])

        piecewise.haz[j] <- result$par
        log.likelihood[j] <- result$value

      }
      #The current time is now the previous time for the next interval
      t.prev <- breakpoints[j]

    }
    #Output the results
    column.names <- colnames(df[,3:ncol.df])
    #Need to change the timepoint text to be "Inf" for technical correctness
    output <- data.frame(switchpoints = gsub(max(df$time), "Inf",column.names),
                         hazards  = piecewise.haz,
                         log.likelihood =log.likelihood )


  } else{
    #If no breakpoints provided by the use just fit a one parameter exponential model
    result <- optim(par=1, fn=jll.expo, method= "L-BFGS-B",
                    lower=c(0.000001), upper=c(100),
                    control=list(fnscale = -1),
                    tt=df$time, status=df[,"status"])

    output <-  data.frame(switchpoints = NA,
                          hazards  = result$par,
                          log.likelihood =result$value )


  }

  return(output)

}


jll.expo <- function(par, tt, status){

  dd <- sum(status)
  sum.t <- sum(tt)

  result <- dd*log(par) - par*sum.t

  result
}


nelder.mead.piecewise.pchreg <- function(time, status,breakpoints){

    #Errors can occur if the algorithm picks a value that is outside the time horizon
    #or if a smaller time is picked for a later interval (or if the same time is picked for both!)

    Surv.formula<- as.formula("Surv(time, status) ~ 1")
    breaks <- c(0,breakpoints, (max(time)+0.01))
    #browser()
    if(all(order(breaks) == 1:length(breaks)) &&
       all(duplicated(breaks)==FALSE) &&
       all(breaks >= 0) &&
       all(diff(breaks) >= 0.01)){

      log.likelihood <- pchreg(Surv.formula,breaks = breaks)$logLik


    }else{
      log.likelihood <- -Inf
    }
    return(log.likelihood)
}


##Need to fix this formula for censors
#Need to figure out why gompertz seems to have stable hazard

KM.hazards <- function(time.vec, status.vec, intervals =5){


  df.KM <-survfit(Surv(time.vec,status.vec) ~1)
  time <- df.KM$time
  tau <- c(df.KM$time[1],diff(df.KM$time))
  risk.time <- df.KM$n.risk*tau
  n.event <- df.KM$n.event

  #We need to find a sequence of 0's and add the risk time to the previous or next event time
  #by default we add to the previous one (execpt if the first observation is a censor)
  n.event.runs <- n.event
  n.event.runs[which(n.event.runs != 0)] <- which(n.event.runs != 0)
  runs <- rle(n.event.runs)
  myruns = which(runs$lengths > 1)

  #Find end of sequence
  runs.lengths.cumsum = cumsum(runs$lengths)
  ends = runs.lengths.cumsum[myruns]

  #Find Start of sequence
  newindex = ifelse(myruns>1, myruns-1, 0)
  starts = runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)

  #Need to sum up the censor risk times and add it to the next event !
  df <- data.frame(starts,ends)

  if(length(ends) !=0){
    drops <- c()
    for(i in 1:nrow(df)){
      if((df[i,2]) == length(risk.time)){
        #if last observations are censors, do nothing!

      }else{
        #else move to the next event
        risk.time[(df[i,2])+1] <- sum(risk.time[df[i,1]:df[i,2]]) + risk.time[(df[i,2])+1]

      }

      drops <- c(drops,df[i,1]:df[i,2])
    }

    risk.time <- risk.time[-drops]
    n.event <- n.event[-drops]
    time <- time[-drops]
    rm(drops)
  }

  #Need to do the same for single censors

  single.censor <- which(n.event == 0)

  drops2 <- c()

  if(length(single.censor)!=0){

    for(i in 1:length(single.censor)){

      if(single.censor[i]==length(risk.time)){


      }else{

        risk.time[single.censor[i]+1] <- risk.time[single.censor[i]] +   risk.time[single.censor[i]+1]
      }

      drops2 <- c(drops2,single.censor[i])


    }

    risk.time <- risk.time[-drops2]
    n.event <- n.event[-drops2]
    time <- time[-drops2]

  }

  df.new <- data.frame(n.event = n.event, risk.time=risk.time, time = time)

  num.rows <- nrow(df.new)
  obs.per.group<- round(num.rows/intervals)


  df.new <- df.new %>%
    group_by(grp = as.integer(gl(n(), obs.per.group, n()))) %>%
    #or with rep
    # group_by(grp = rep(row_number(), length.out = n(), each = 5))
    summarise(sum.event = sum(n.event),
              sum.risk = sum(risk.time ),
              median.time = median(time),
              first.time = first(time),
              last.time = last(time)) %>%
    mutate(avg.haz =sum.event/sum.risk )

lm.reg <-lm(log(df.new$sum.event/df.new$sum.risk)~ df.new$median.time)

 #https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph

  eq <- as.expression(substitute(~~italic(r)^2~"="~r2,
                    list(r2 = format(summary(lm.reg)$r.squared, digits = 3))))

  haz.plot <-ggplot(df.new,
                     aes(x = first.time, y = log(avg.haz) ))+
    geom_step(linetype = 3)+
    geom_point(inherit.aes = F, aes(x = median.time, y = log(avg.haz)))+
    xlab("time") + ylab("log(Average hazard)")+
    geom_text(x = sum(as.numeric(c(df.new[1,"first.time"],
                                   df.new[1,"last.time"])))/2, y = mean(log(df.new$avg.haz)),
              label = as.character(eq), parse = TRUE)

  return(haz.plot)

  #output <- list()
  #output[[1]] <- df.new
  #output[[2]] <- plot.haz

  #return(output)

}

lm_eqn <- function(y, x){
  #https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph

  m <- lm(y ~ x);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

rllogis2 <- function(n, lambda, alpha){
  theta <- log(lambda)
  surv.prob <-runif(n, 1, 99)


  surv.times <- (surv.prob*(exp(-theta))/(100-surv.prob))^(1/alpha)

  return(surv.times)
}


pwe <- function(data, timevar, deathvar, bounds) {
  # pwe: expands an S data frame for piece-wise exponential survival
  # G. Rodriguez, Nov 29, 1992
  #
  # Check arguments: time and death must be variables in the data frame
  # and boundaries must be non-negative and strictly increasing

  #source("http://data.princeton.edu/wws509/stata/pwe.txt")
  if(!is.data.frame(data)) stop("First argument must be a data frame")
  if(is.na(match(tn <- deparse(substitute(timevar)), names(data))))
    stop(paste("\n\tSurvival time", tn,
               "must be a variable in the data frame"))
  if(is.na(match(dn <- deparse(substitute(deathvar)), names(data))))
    stop(paste("\n\tDeath indicator", dn,
               "must be a variable in the data frame"))
  width <- diff(bounds)
  if(any(bounds < 0) | any(width <= 0)) stop(paste(
    "Invalid interval boundaries in", deparse(substitute(
      bounds))))      #
  # Expand the data frame creating one pseudo-observation for each
  # interval visited, add interval number, events and exposure time
  # (existing variables with these names will be overwriten)
  n <- cut(data[, tn], bounds)
  data <- data[rep(seq(along = n), n),  ]
  i <- NULL
  for(k in 1:length(n))
    i <- c(i, 1:n[k])
  data$events <- ifelse(data[, tn] > bounds[i + 1], 0, data[, dn])
  data$exposure <- ifelse(data[, tn] > bounds[i + 1], width[i], data[, tn
                                                                     ] - bounds[i])
  data$interval <- i
  attr(data$interval, "levels") <- attr(n, "levels")
  data
}



#New piewcewise exponentential function

piecewise.exponential2 <- function(breakpoints = NULL, time , status){

  time <- time
  status <- status
  df <- data.frame(time= time,
                   status = status)
  df <- df[order(time), ]

  df$id <- 1:nrow(df)

  ##Note some of the code is redundant and refers to when the Nelder-Mead algorithm was applied.
  #Errors can occur if the algorithm picks a value that is outside the time horizon
  #or if a smaller time is picked for a later interval (or if the same time is picked for both!)
  if(length(breakpoints)>0){

    breaks <- c(0,breakpoints, max(time)+.001)

    #Because we use left breakpoints the last observation is NA, could avoid this
    #by adding an "epsilon" (small value to ensure it is above max value)
    df$timepoints <- cut(df$time, breaks = breaks, right = FALSE)

    #Create a dataframe with 1's and 0's for each interval
    df <- cbind(df, model.matrix(~ timepoints + 0, data=df))
    #Drop redundant columns
    df <- df[,-which(colnames(df) == "timepoints")]

    df.likelihood <- df %>% select(time, status, id)

    #Make sure that middle intervals at least 1 event
    if(length(breakpoints)> 1 && any(apply(df[,3:(ncol(df)-1)],2, function(x)all(x== 0)))){

      likelihood <- NA
      return(likelihood)


    }else{

      #Define some placeholder vectors
      ncol.df <- ncol(df)
      piecewise.haz <- rep(NA, length(breakpoints)+1)

      nrows <- nrow(df)

      #Define some counters
      j <- 0
      t.prev <- 0
      #For 4rd column (first time interval) to the final time interval
      for(i in 4:ncol.df){

        j <- j + 1
        # browser()
        temp.df <- df[min(which(df[,i] == 1)):nrow(df),]
        nrows <- nrow(temp.df)
        end.point <- max(which(temp.df[,i] == 1))

        #If we are in an interval before the last one all observations
        # after the breakpoints are censored
        if(end.point < nrows){
          temp.df[end.point:nrows, "status"] <- 0
          temp.df[end.point:nrows, "time"] <- breakpoints[j]
        }
        # Center the time to be from the start of that interval
        temp.df$time <-  temp.df$time - t.prev

        haz <- sum(temp.df[,"status"])/sum(temp.df$time)

        temp.df <- temp.df %>% mutate(likelihood = (haz^status)*(exp(-haz*time))) %>%
          dplyr::select(id, likelihood)

        names(temp.df)[names(temp.df) == "likelihood"] <- paste0("likelihood",(i-3))

        df.likelihood <- df.likelihood %>% left_join(temp.df)
        #As defined in the chapter the log likelihood of a exponential constant hazard is
        #(Number of Events)*log(Number of Events/Exposure time in the interval)

        #log.likelihood.new <- sum(temp.df[,"status"])*log(sum(temp.df[,"status"])/sum(temp.df$time))


        #If no events in the final interval then the likelihood is zero?? NO this isn't true
        #if(is.nan(log.likelihood.new)){
        #log.likelihood.new <- 0
        #}


        #log.likelihood <- log.likelihood.new + log.likelihood
        #The current time is now the previous time for the next interval
        t.prev <- breakpoints[j]

      }

      df.likelihood[is.na(df.likelihood)] <- 1

      df.likelihood[,grep( "lik",colnames(df.likelihood))] <- log(df.likelihood[,grep( "lik",colnames(df.likelihood))])
      log.likelihood.value <- sum(rowSums(df.likelihood[,grep( "lik",colnames(df.likelihood))]))

      #likelihood <- exp(as.brob(log.likelihood.value))

      return(log.likelihood.value)

    }
  }else{
    #If no breakpoints provided by the use just fit a one parameter exponential model
    log.likelihood.value <- sum(log((haz^df$status)*(exp(-haz*df$time))))

    #likelihood <- exp(as.brob(log.likelihood.value))
    return(log.likelihood.value)

  }

}

#Hazard Plot Test
hazard.plot <- function(time.vec, cens.vec, lrg.haz.int = 1,  bw.method = "local"){

  #Plot the hazards
  max.time <- max(time.vec)


  result.hazard.pe.lrg <- pehaz(time.vec, cens.vec, width= lrg.haz.int, max.time=max.time)

  result.smooth <- muhaz(time.vec, cens.vec, b.cor="left",
                         max.time=max.time, bw.method = bw.method)

  #Plot the Survival function

  plot(result.hazard.pe.lrg,  col="black")
  lines(result.smooth, col = "red", lty= 2)

  #return(result.smooth)
}


#Gibbs functions

#Maximum Likelihood and likelihood with supplied hazards

piecewise_loglik <- function(df, changepoint, method = "ML", lambda = NULL) {

  surv.object <- with(df, Surv(enter, time, status))
  split <- SurvSplit(surv.object, changepoint)
  n.ivl <- length(changepoint) + 1

  T <- split$Y$exit - split$Y$enter
  d <- split$Y$event
  ivl <- split$ivl

  if(method == "ML"){

    alpha <- matrix(0, nrow = 1, ncol = n.ivl)
    Dtot <- sum(d)
    res <- -Dtot
    for (i in seq_len(n.ivl)) {
      indx <- ivl == i
      D <- sum(d[indx])
      if (D > 0) {
        sumT <- sum(T[indx])
        alpha[i] <- D/sumT
        res <- res + D * (log(D) - log(sumT))
      }
    }
    #return(list(loglik = res, hazards = alpha))
    return(res)
  }else{

    death.loglik <-0
    cum_haz.loglik <- 0

    for (i in seq_len(n.ivl)) {
      indx <- ivl == i
      death.loglik <-  sum(d[indx])*log(lambda[i]) + death.loglik
      cum_haz.loglik <- -sum(lambda[i]*T[indx]) +cum_haz.loglik

    }
    return(death.loglik +cum_haz.loglik)
  }
}


#Time and death intervals
exposure_death <- function(df, changepoint) {

  surv.object <- with(df, Surv(enter, time, status))
  split <- SurvSplit(surv.object, changepoint)
  n.ivl <- length(changepoint) + 1

  T <- split$Y$exit - split$Y$enter
  d <- split$Y$event
  ivl <- split$ivl

  result <- matrix(NA, ncol = 2, nrow = n.ivl)
  for (i in seq_len(n.ivl)) {
    indx <- ivl == i
    result[i,1] <- sum(d[indx])
    result[i,2] <- sum(T[indx])
  }
  return(result)
}

