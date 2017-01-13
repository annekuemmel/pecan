parameter.estimation.nls <- function(

# Parameter estimation and function CI assessment based on delta method
# inputs: 
# - model: model function with x as first input argument 
#          followed by the model parameters
#          (required)
# - init: named list with initial guesses of parameter values
#          (required)
# - data: data.frame or list containing independent (x) and
#          dependent variable (y) respectively
#          (required), x may contain additional design variables
#
# Anne Kuemmel, Oct 2015
 
   model = function(
         x=seq(0,1,0.1),
         E0 = 0.2,
         Emax = 1,
         EC50 = 0.5
      ){
         y <- (Emax - E0) * x/(EC50 + x) + E0
         return(y)
      },
   init = list(E0=0.5,Emax=1.5,EC50=0.4),
   data = list(x=10^seq(-2,1,0.5),
                     y=model(x=10^seq(-2,1,0.5)*exp(rnorm(n=length(seq(-2,1,0.5)),sd=0.1)))),

   ...

){

   ###########################################
   # Functionality
   require(plyr)

   ###########################################
   # Preparations
   model.args <- names(formals(model))
   pars <- setdiff(model.args,'x')

   # get inits in same order as given in model
   init <- init[pars]
               
   ###########################################
   # Definitions
   m <- length(pars)   

   ###########################################
   # Parameter estimation

   # define strings to construct expression to estimate parameters given the initials
   parstr <- paste(pars,collapse=',')
   initstr <- paste(pars,unlist(init),collapse=',',sep='=')

   # do estimation with nls
   stime <- system.time(
      suppressWarnings(
           estres <- eval(parse(text=paste('nls(y ~ model(x,',parstr,'), data = data, start = list(',initstr,'), ...)',sep='')))
      )
   )

   # extract results

     # degrees of freedom
     DF <- df.residual(estres)

     # parameter estimates
     estimates <- coef(estres)
     # add residual error estimate
     resSD <- sqrt(sum(residuals(estres)^2)/DF) # OR sd(residuals(estres)) ?
     estimates <- c(coef(estres),resSD)
     names(estimates) <- c(pars,'a')
   
     # variance-covariance matrix of estimates
     estvar    <- vcov(estres)

     # standard error of estimates
     SE <- sqrt(diag(estvar))
     SEresSD <- resSD / sqrt(2*DF) # this is an approximation for sufficiently large DF!!
     SE <- c(SE,"a"=SEresSD)



   # predicted values
   ypred <- fitted(estres)
   # residuals
   res <- residuals(estres)
   # weighted residuals
   wres <- res / sd(residuals(estres))

   # assemble fit diagnostics
   if ( is.list(data$x) ) { diagn <- data.frame(t=data$x[[1]],y=data$y,ypred=ypred,res=res,wres=wres)
   }  else  { diagn <- data.frame(t=data$x,y=data$y,ypred=ypred,res=res,wres=wres) }

   ###########################################
   # Output
   OUT <- list(
            estmethod = 'nls',
            errmod = 'add',
            estimates = estimates,
            errval = sd(residuals(estres)),
            SE = SE,
            estvar = estvar,
            diagn=diagn,
            estimation = estres
          )

   cat('Time elapsed for estimation is',stime[['elapsed']],'sec.\n')

   return(OUT)

}
