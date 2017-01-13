parameter.estimation.nlm <- function(

# Parameter estimation using nonlinear likelihood maximation: 
# - model: model function with scalar output and input arguments
#          including x and the model parameters
# - errmod: error model ('add', 'prop', 'comb', 'exp' possible)
# - init: named list with initial guesses of parameter values
# - data: data.frame containing independent (xcol) and
#          dependent variable (ycol) respectively

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
   errmod = 'add',

   init = list(E0=0.5,Emax=1.5,EC50=0.4),

   data = list(x=10^seq(-2,1,0.5),y=model(x=10^seq(-2,1,0.5))*exp(rnorm(n=length(seq(-2,1,0.5)),sd=0.1))),

   init.err=NULL,
   ...

){

   ###########################################
   # Functionality
   require(plyr)

   ###########################################
   # Preparations

   # get parameter names from model function
   model.args <- names(formals(model))
   pars <- setdiff(model.args,'x')

   # get inits in same order as given in model
   init <- init[pars]

   # check whether initial for residuals error is given
   if (!is.null(init.err)) {
      cat('\nUser-specified initial value for residual error detected:\n')
      if (errmod == 'comb') {
        if (length(init.err) != 2) {
           cat('  Sorry, you need to specify TWO initial values for residual errors for combination error model.\n  Default (c(1,1)) set as initials.\n')
           init.err <- c(1,1)
        } else {cat('  Accepted.\n')}
      } else {
        if (length(init.err) != 1) {
           cat('  Sorry, you need to specify ONE initial value for residual errors for chosen error model.\n  First set as initials.\n')
           init.err <- init.err[1]
        } else {cat('  Accepted.\n')}
      }
   } else {
      if (errmod == 'comb') {
           init.err <- c(1,1)
      } else {
           init.err <- 1
      }
   }

   # check whether dependent variables are positive if the error model is exponential
   if (errmod == 'exp' & any(data$y <= 0)) stop('\nThis does not work: Negative or zero observations but exponential error model.\n', call.=FALSE)
               
   ###########################################
   # Definitions
   n <- length(data$y)    # number of observations
   nphi <- length(pars)   # number of parameter w/o error model

   ###########################################
   # Model reformulation

   # (need to pack all parameters into one variable p when using nlm function) 
   modelTP <- function(t,p){
       pars <- setdiff(names(formals(model)),'x')
       df <- data.frame(pars=pars,no=paste('p[',1:length(pars),']',sep=''))
       transstr <- paste(df$pars,df$no, sep = '=', collapse = ',')
       eval(parse(text=paste('y <- model(x=t,',transstr,' )',sep='')))
       return(y)
   }

   ###########################################
   # Parameter estimation

   # define score function (-2LL) depending on the error model
   if (errmod == 'add'){
      fmin <- function(p,y,t){
         nphi <- length(p)-1
         yhat <- modelTP(t, p)
         g <- p[nphi+1]
         e <- sum( ((y-yhat)/g)^2 + log(g^2) )
         return(e)
      }
      start <- c(unlist(init),init.err)
   }
   if (errmod == 'prop') {
      fmin <- function(p,y,t){
         nphi <- length(p)-1
         yhat <- modelTP(t, p)
         g <- p[nphi+1]*yhat
         e <- sum( ((y-yhat)/g)^2 + log(g^2) )
         return(e)
      }
      start <- c(unlist(init),init.err)
   }
   if (errmod == 'comb') {
      fmin <- function(p,y,t){
         nphi <- length(p)-2
         yhat <- modelTP(t, p)
         g <- abs(p[nphi+1]) + abs(p[nphi+2])*yhat
         e <- sum( ((y-yhat)/g)^2 + log(g^2) )
         return(e)
      }
      start <- c(unlist(init),init.err)
   }
   if (errmod == 'exp') {
      fmin <- function(p,y,t){
         nphi <- length(p)-1
         yhat <- modelTP(t, p)
         g <- p[nphi+1]
         e <- sum( ((log(y)-log(yhat))/g)^2 + log(g^2) )
         return(e)
      }
      start <- c(unlist(init),init.err)
   }
   if (!exists('fmin')) stop('Score function for MLE not defined. Given error model unknown.')
   npsi <- length(start) # number of parameters including error model

   # do estimation with nlm
   stime <- system.time(
      suppressWarnings(
         estres <- nlm(fmin, start, y=data$y,t=data$x, hessian = TRUE, ...)
      )
   )

   # get -2LL
   if (errmod == 'exp')
   {
     estres$neg2LL <- estres$minimum + n*log(2*pi) + 2*sum(log(data$y))
   } else {
     estres$neg2LL <- estres$minimum + n*log(2*pi)
   }
   # get BIC
   if (errmod == 'comb')
   {
     estres$BIC <- estres$neg2LL + log(n)*(nphi+2)
   } else {
     estres$BIC <- estres$neg2LL + log(n)*(nphi+1)
   }

   # make sure that standard deviations for errors are positive
   if (errmod == 'comb') estres$estimate[npsi-c(1,0)] <- abs(estres$estimate[npsi-c(1,0)])
   if (errmod != 'comb') estres$estimate[npsi] <- abs(estres$estimate[npsi])

   # parameter estimates
   estimates <- estres$estimate  # parameter estimates
   if (errmod == 'add') {names(estimates) <- c(pars,'a'); errval <- estimates['a']}
   if (errmod == 'prop') {names(estimates) <- c(pars,'b'); errval <- estimates['b']}
   if (errmod == 'comb') {names(estimates) <- c(pars,'a','b'); errval <- estimates[c('a','b')]}
   if (errmod == 'exp') {names(estimates) <- c(pars,'e'); errval <- estimates['e']}

   # variance-covariance matrix of estimates
   FIM <- estres$hessian/2  # need to devide by 2 as 2LL was maximized
   estvar <- solve(FIM)
   dimnames(estvar) <- list(names(estimates),names(estimates))

   # standard error of estimates
   SE <- sqrt( diag(estvar) * n/(n-npsi) )

   # predicted values
   ypred <- modelTP(t=data$x, p=estres$estimate)

   # residuals and weighted residuals
   if (errmod == 'add') {res <- data$y-ypred; wres  <- res / (estimates[['a']])}
   if (errmod == 'prop') {res <- data$y-ypred; wres  <- res / sqrt( (ypred*estimates[['b']])^(2) )} 
   if (errmod == 'comb') {res <- data$y-ypred; wres  <- res / sqrt( (estimates[['a']] + estimates[['b']]*ypred)^(2) )} 
   if (errmod == 'exp') {res <- log(data$y)-log(ypred); wres  <- res / estimates[['e']] }

   # assemble fit diagnostics to data frame
   if ( is.list(data$x) ) { diagn <- data.frame(t=data$x[[1]],y=data$y,ypred=ypred,res=res,wres=wres)
   }  else  { diagn <- data.frame(t=data$x,y=data$y,ypred=ypred,res=res,wres=wres) }

   ###########################################
   # Output
   OUT <- list(
            estmethod = 'nlm',
            errmod = errmod,
            estimates = estimates,
            errval = errval,
            SE = SE,
            estvar = estvar,
            diagn=diagn,
            estimation = estres
          )

   cat('\nTime elapsed for estimation is',stime[['elapsed']],'sec.\n')

   return(OUT)

}
