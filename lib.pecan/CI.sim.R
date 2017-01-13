CI.sim <- function(

      model,     # model with inputs x and the parameters separately
      parval,    # named vector or list of parameter values
      parvar,    # variance-covariance matrix of the parameters

      errval = 0,     # residual error (need to be 2 values for comb error model, first one additive)
      errmod = 'add', #

      xsupport,       # x values to evaluate CI
      xcol = 'time',  # name of column in x which is support dimension
      n = Inf,        # number of observations (Inf indicates non-estimated variance)
      CI.level = 0.9,  # confidence interval needs to be between 0 and 100
      nsim = 1000     # number of simulations

   
){

   ###########################################
   # Preparations

   # Load MASS package for sampling from multivariate normal distribution
   require(MASS)

   # get parameters
   model.args <- names(formals(model))
   pars <- setdiff(model.args,'x')

   # get parameter values in same order as given in model 
   parvalS <- parval[pars]

   # make sure that parval is a list
   parval <- as.list(parval)
   parvalS <- as.list(parvalS)

   # get inits in same order as given in model (needed since sampling from 
   parvarS <- parvar[pars,pars]

   # general definitions
         np <- length(pars)
         df <- n-np # degrees of freedom

   # determine quantiles corresponding to confidence level
   ci.quantiles <- c((1-CI.level)/2,1-(1-CI.level)/2)

   # translate quantiles to multiple of standard deviation
   # (automatically gets quantiles from normal distribution for df = Inf)
   fac <- qt(ci.quantiles, df = df) 

   ###########################################

   # * * * * * * * * * * * * #
   # 0) Determine confidence interval for parameters


      # get standard deviation for predictions from diagonal
      if (n == Inf){ # -> variance is not estimated, do not correct for bias
         sd <- sqrt(diag(parvar))
      } else { # correct for bias of estimated variance
         sd <- sqrt(diag(parvar) * n/df)
      }

      # calculate CI
      estimates <- unlist(parval[dimnames(parvar)[[1]]])
      CI.estimates <- cbind(dimnames(parvar)[[1]],as.data.frame(estimates + as.matrix(sd) %*% t(as.matrix(fac))))
      names(CI.estimates) <- c('parameter','CI.lb','CI.ub')

   # * * * * * * * * * * * * #
   # 1) Variance matrix to sample parameters from

#      if (n == Inf){ # -> variance is not estimated, do not correct for bias
         samplemat <- parvarS
#      } else { # correct for bias of estimated variance (OBS!! this does not work as sometimes non-positive matrices are generated)
#         samplemat <- parvarS * n/df
#      }  

   # * * * * * * * * * * * * #
   # 2) Sample parameters

      parsample <- mvrnorm(n=nsim, mu=unlist(parvalS), Sigma=samplemat)

   # * * * * * * * * * * * * #
   # 3) predict and simulate

      # model reformulation
      modelTP <- function(t,p){
         pars <- pars
         df <- data.frame(pars=pars,no=paste('p[',1:length(pars),']',sep=''))
         transstr <- paste(df$pars,df$no, sep = '=', collapse = ',')
         eval(parse(text=paste('y <- model(x=t,',transstr,' )',sep='')))
         return(y)
      }

      # predictions for point estimate
      ifelse(class(xsupport) == 'list', xsupport.vector <- xsupport[[xcol]], xsupport.vector <- xsupport)
      nsupp <- length(xsupport.vector)
      pred0 <- data.frame(x=xsupport.vector, y=modelTP(t=xsupport,p=unlist(parvalS)))

      # predictions (w/o residual error)
      pred <- matrix(NA,nsupp,nsim)
      for (i in 1:nsim){
         pred[,i] <- modelTP(t=xsupport,p=parsample[i,])
      }

      # simulations (w residual error)
      if (all(errval == 0)) {
         sim <- NULL
      } else {

         # sample residual error
         eps <- matrix(rnorm(nsupp*nsim),nsupp,nsim)

         # get sigma per prediction
         switch(errmod,
            add = g <- matrix(errval,nsupp,nsim),
            prop = g <- pred * matrix(errval,nsupp,nsim),
            comb = g <- matrix(errval[1],nsupp,nsim) + pred * matrix(errval[2],nsupp,nsim),
            exp = g <- matrix(errval,nsupp,nsim)
         )
         
         # add residual error on predictions
         if (errmod == 'exp'){
            sim <- exp( log(pred) + g*eps )
         } else {
            sim <- pred + g*eps
         }

      }
      n.fail <- sum(apply(is.na(pred),2,all))

   # * * * * * * * * * * * * #
   # 4) Confidence interval for the prediction
   #    i.e. percentiles of predictions for estimates from all simulated datasets

   # get empirical percentiles 
   CI <- as.data.frame(cbind(xsupport.vector,t(apply(pred,1,quantile,probs = ci.quantiles, na.rm=TRUE))))
   names(CI) <- c('x','CI.lb','CI.ub')


   # * * * * * * * * * * * * #
   # 5) Prediction interval of the simulations
   #    i.e. percentiles of simulations for estimates from all bootstrapped data

   # get empirical percentiles
   if (is.null(sim)){
      PI <- NULL
   } else {
      PI <- as.data.frame(cbind(xsupport.vector,t(apply(sim,1,quantile,probs = ci.quantiles, na.rm=TRUE))))
      names(PI) <- c('x','PI.lb','PI.ub')
   } 


   # output
   OUT <- list(
               CImethod='simulation',
               CI.estimates=CI.estimates,
               pred = pred0,
               CI=CI,
               PI=PI,
               n = nsim,
               n.fail=n.fail
   )
   return(OUT)

}
