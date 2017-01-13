CI.boot <- function(

      model,             # model with inputs x and the parameters separately
      errmod = 'add',    # error model
      init,              # initial values to start estimations with bootstrapped data

      data,              # data...
      xcol = 'x',
      ycol = 'y',

      xsupport,          # x values to evaluate CI
      CI.level = 0.9,    # confidence interval needs to be between 0 and 1
      nboot = 1000,      # number of simulations

      estmethod = 'nls', # choose estimation method 'nls' or 'nlm'

      ...
){
   
   ###########################################
   # Preparations

   # get parameters
   model.args <- names(formals(model))
   pars <- setdiff(model.args,'x')

   # number of observations
   n <- length(data$y)

   # make sure that parval is a list
   init <- as.list(init)

   # simulation for initials (assumed to be the point estimates)
   parstr <- paste(pars,pars,collapse=',',sep='=init$')

   ifelse(class(xsupport) == 'list', xsupport.vector <- xsupport[[xcol]], xsupport.vector <- xsupport)
   nsupp <- length(xsupport.vector)
   pred0 <- eval(parse(text=paste('data.frame(x=xsupport.vector, y = model(x=xsupport,',parstr,'))',sep='')))

   # choose estimation method
   switch(estmethod,
      nls = estfun <- est.nls,
      nlm = estfun <- est.nlm
   )
   if (!exists('estfun')) stop('Estimation method not defined.')

   # prepare simulation
   estparstr <- paste(pars,pars,collapse=',',sep='=est_i$')

   # prepare containers for results
   estimates <- data.frame()
   simulations <- predictions <- matrix(NA,nrow=nsupp,ncol=nboot)
   n.succ <- n.fail <- 0

   ###########################################

   for (i in seq(1,nboot)) {

   #cat("iteration = ", i, "\n")

   # * * * * * * * * * * * * #
   # 1) Bootstrap data

      # simulate data set
      if ('list' %in% class(data$x)) {
        bootdata_i <- data
        idx_smpl <- sample(n,rep = TRUE)
        bootdata_i$x[[xcol]] <- data$x[[xcol]][idx_smpl]
        bootdata_i$y <- data$y[idx_smpl] }
      else {
        bootdata_i <- data
        idx_smpl <- sample(n,rep = TRUE)
        bootdata_i$x <- data$x[idx_smpl]
        bootdata_i$y <- data$y[idx_smpl] }

   # * * * * * * * * * * * * #
   # 2) Estimate parameters with bootstrapped data

      if (i==1) {
         esttime_one <- system.time(est_i <- estfun(model, bootdata_i, init, errmod, ...))
         cat("\nEstimated time for bootstrapping:",format(nboot*esttime_one[['elapsed']]/60,digits=2),"min.\n")
      } else
         est_i <- estfun(model, bootdata_i, init, errmod, ...)

      if (!is.null(est_i)) {

         n.succ <- n.succ + 1
         estimates <- rbind(estimates,as.data.frame(est_i))

   # * * * * * * * * * * * * #
   # 3) Predict and simulate

         # predict -> confidence interval
         pred_i <- eval(parse(text=paste('y = model(x=xsupport,',estparstr,')',sep='')))
         predictions[,i] <- pred_i
         
         # simulate -> prediction interval (containing residual error)
         switch(errmod,
                add  = sim_i <- pred_i + rnorm(nsupp,sd=est_i$a),
                prop = sim_i <- pred_i + pred_i*rnorm(nsupp,sd=est_i$b),
                comb = sim_i <- pred_i + rnorm(nsupp,sd=est_i$a) + pred_i*rnorm(nsupp,sd=est_i$b),
                exp  = sim_i <- pred_i * rlnorm(nsupp,meanlog=0,sdlog=est_i$e)
         )
         simulations[,i] <- sim_i

      } else {
         n.fail <- n.fail + 1
      }

   rm(bootdata_i,est_i)
   }

   # * * * * * * * * * * * * #
   # 4) Confidence interval of parameter estimates
   #    i.e. percentiles of estimates from all bootstrapped data

   # determine quantiles corresponding to confidence level
   ci.quantiles <- c((1-CI.level)/2,1-(1-CI.level)/2)
   
   # get empirical quantiles for parameter estimates
   CI.estimates <- as.data.frame(t(apply(estimates,2,quantile,probs = ci.quantiles)))
   names(CI.estimates) <- c('CI.lb','CI.ub')
   CI.estimates$parameter <- rownames(CI.estimates)
   CI.estimates <- CI.estimates[,c('parameter','CI.lb','CI.ub')]

   # * * * * * * * * * * * * #
   # 5) Confidence interval for the prediction
   #    i.e. percentiles of predictions for estimates from all bootstrapped data

   # get empirical quantiles for predictions
   CI <- as.data.frame(cbind(xsupport.vector,t(apply(predictions,1,quantile,probs = ci.quantiles, na.rm=TRUE))))
   names(CI) <- c('x','CI.lb','CI.ub')

   # * * * * * * * * * * * * #
   # 6) Prediction interval of the simulations
   #    i.e. percentiles of simulations for estimates from all bootstrapped data

   # get empirical quantiles 
   PI <- as.data.frame(cbind(xsupport.vector,t(apply(simulations,1,quantile,probs = ci.quantiles, na.rm=TRUE))))
   names(PI) <- c('x','PI.lb','PI.ub')



# output
   OUT <- list(
               CImethod='bootstrap',
               CI.estimates=CI.estimates,
               pred=pred0,
               CI=CI,
               PI=PI,
               n = nboot,
               n.fail=n.fail
   )
   return(OUT)

}
