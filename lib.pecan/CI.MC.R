CI.mc <- function(

      model,     # model with inputs x and the parameters separately
      parval,    # named vector or list of parameter values

      errval = 0,     # residual error (need to be 2 values for comb error model, first one additive)
      errmod = 'add', #

      xsample,        # sampling timepoints
      xsupport,       # x values to evaluate CI
      xcol = 'time',  # name of column in x which is support dimension
      CI.level = 0.9,  # confidence interval needs to be between 0 and 1
      nmc = 1000,     # number of simulations

      estmethod = 'nls', # choose estimation method 'nls' or 'nlm'

      ...
){

   ###########################################
   # Preparations

   # get parameters
   model.args <- names(formals(model))
   pars <- setdiff(model.args,'x')

   # make sure that parval is a list
   parval <- as.list(parval)

   # choose estimation method
   switch(estmethod,
      nls = estfun <- est.nls,
      nlm = estfun <- est.nlm
   )
   if (!exists('estfun')) stop('Estimation method not defined.')

   # determine quantiles corresponding to confidence level
   ci.quantiles <- c((1-CI.level)/2,1-(1-CI.level)/2)

   # prepare x support vector
   ifelse(class(xsupport) == 'list', xsupport.vector <- xsupport[[xcol]], xsupport.vector <- xsupport)
   nsupp <- length(xsupport.vector)

   # prepare x sample vector
   ifelse(class(xsample) == 'list', xsample.vector <- xsample[[xcol]], xsample.vector <- xsample)
   nsmpl <- length(xsample.vector)

   # containers
   estimates <- data.frame()
   simulations <- predictions <- matrix(NA,nrow=nsupp,ncol=nmc)
   n.succ <- n.fail <- 0

   ###########################################

   # predictions for given parameter values
   predstr <- paste(pars,pars,collapse=',',sep='=parval$')

      # ... for x values for which to calculate predictions and their intervals
      pred0 <- eval(parse(text=paste('data.frame(x=xsupport.vector, y = model(x=xsupport,',predstr,'))',sep='')))

      # ... for x values of the data to prepare simulated data sets (for which only residual error is added)
      pred.sample <- eval(parse(text=paste('data.frame(x=xsample.vector, y = model(x=xsample,',predstr,'))',sep='')))


   # prepare simulation
   estparstr <- paste(pars,pars,collapse=',',sep='=est_i$')


   # * * * * * * * * * * * * #
   # 1) Create simulated data

      # simulate data set
      # (add random residual error)
         yhat <- matrix(pred.sample$y,nrow=nsmpl,ncol=nmc,byrow=FALSE)
         switch(errmod,
                add  = simdata <- yhat + matrix(rnorm(nsmpl*nmc,sd=errval),nrow=nsmpl,ncol=nmc),
                prop = simdata <- yhat + yhat*matrix(rnorm(nsmpl*nmc,sd=errval),nrow=nsmpl,ncol=nmc),
                comb = simdata <- yhat + rnorm(nsmpl*nmc,sd=errval[1]) + yhat*rnorm(nsmpl*nmc,sd=errval[2]),
                exp  = simdata <- yhat * rlnorm(nsmpl*nmc,sdlog=abs(errval))
         )

   for (i in seq(1,nmc)) {

   # * * * * * * * * * * * * #
   # 2) Estimate parameters for each simulated data set

      data_i <- list(x=xsample,y=simdata[,i])

      if (i==1) {
         esttime_one <- system.time(est_i <- estfun(model, data_i, parval, errmod, ...))
         cat("\nEstimated time for MC simulation/estimation:",format(nmc*esttime_one[['elapsed']]/60,digits=2),"min.\n",sep="")
      } else
         est_i <- estfun(model, data_i, parval, errmod, ...)

      if (!is.null(est_i)) {

         n.succ <- n.succ + 1
         estimates <- rbind(estimates,as.data.frame(est_i))

   # * * * * * * * * * * * * #
   # 3) Predict and simulate for estimated parameters on simulated dataset

         # predict -> confidence interval
         pred_i <- eval(parse(text=paste('y = model(x=xsupport,',estparstr,')',sep='')))
         predictions[,i] <- pred_i
         
         # simulate -> prediction interval
         switch(errmod,
                add  = sim_i <- pred_i + rnorm(nsupp,sd=est_i$a),
                prop = sim_i <- pred_i + pred_i*rnorm(nsupp,sd=est_i$b),
                comb = sim_i <- pred_i + rnorm(nsupp,sd=est_i$a) + pred_i*rnorm(nsupp,sd=est_i$b),
                exp  = sim_i <- pred_i * rlnorm(nsupp,meanlog=0,sdlog=est_i$e)
         )
         simulations[,i] <- sim_i
         rm(sim_i,est_i)

      } else {
         n.fail <- n.fail + 1
         rm(est_i)
      }

   }

   # * * * * * * * * * * * * #
   # 4) Confidence interval of parameter estimates
   #    i.e. percentiles of estimates from all simulated datasets

   # get empirical quantiles of parameter estimates
   CI.estimates <- as.data.frame(t(apply(estimates,2,quantile,probs = ci.quantiles)))
   names(CI.estimates) <- c('CI.lb','CI.ub')
   CI.estimates$parameter <- rownames(CI.estimates)
   CI.estimates <- CI.estimates[,c('parameter','CI.lb','CI.ub')]

   # * * * * * * * * * * * * #
   # 5) Confidence interval for the prediction
   #    i.e. percentiles of predictions for estimates from all bootstrapped data

   # get empirical quantiles of predictions
   CI <- as.data.frame(cbind(xsupport.vector,t(apply(predictions,1,quantile,probs = ci.quantiles, na.rm=TRUE))))
   names(CI) <- c('x','CI.lb','CI.ub')


   # * * * * * * * * * * * * #
   # 6) Prediction interval of the simulations
   #    i.e. percentiles of simulations for estimates from all simulated datasets

   # get empirical quantiles 
   PI <- as.data.frame(cbind(xsupport.vector,t(apply(simulations,1,quantile,probs = ci.quantiles, na.rm=TRUE))))
   names(PI) <- c('x','PI.lb','PI.ub')



# output
   OUT <- list(
               CImethod = 'MC simulation-estimation',
               CI.estimates=CI.estimates,
               pred = pred0,
               CI=CI,
               PI=PI,
               n = nmc,
               n.fail=n.fail
   )
   return(OUT)

}
