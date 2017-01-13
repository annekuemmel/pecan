CI.delta <- function(

      model,              # model with inputs x and the parameters separately
      model.deriv = NULL, # partial derivatives of model wrt parameters
      errmod = "add",     # error model
      parval,             # named vector or list of parameter values
      parvar,             # variance-covariance matrix of the parameters
      xsupport,           # x values to evaluate CI
      n = Inf,            # number of observations (Inf indicates non-estimated variance)
      CI.level = 0.9,     # confidence interval needs to be between 0 and 1
      xcol = 'time'       # name of x variable within xsupport (others are design variables)

){

cat('\nINFO for delta method: Confidence of parameters and function calculated\nassuming normal or t distribution.\n')

   ###########################################
   # Preparations

   # get parameters
   model.args <- names(formals(model))
   pars <- setdiff(model.args,'x')
 
   # make sure that parval is a list
   parval <- as.list(parval)

   # general definitions
   np <- length(pars)
   df <- n-np # degrees of freedom

   # determine quantiles corresponding to confidence level
   ci.quantiles <- c((1-CI.level)/2,1-(1-CI.level)/2)

   # translate quantiles to multiple of standard deviation
   # (automatically gets quantiles from normal distribution for df = Inf)
   fac <- qt(ci.quantiles, df = df) 

   # string to assign given parameter values to parameter variables
   parstr <- paste(pars,pars,collapse=',',sep='=parval$')

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
   # 1) derivative of the functions wrt parameter

   if (is.function(model.deriv)) { 
      # -> function for derivative of model wrt to parameters given

      # check whether arguments of function and derivative function agree
      model.deriv.args <- names(formals(model.deriv))
      if (!identical(model.args,model.deriv.args)) stop('Arguments of given derivative are not same as for model.')

      # calculate Jacobian using given derivative
      jacob <- eval(parse(text=paste('model.deriv(x=xsupport,',parstr,')',sep='')))
      dimnames(jacob)[[2]] <- aaply(dimnames(jacob)[[2]],1, function(x) substr(x, 2, nchar(x)))

   } else {

      # -> function for derivative not given, use numerical derivative
      jacob <- calc.jacobian(model, parval[pars], xsupport)

   }

   # * * * * * * * * * * * * #
   # 2) calculate variance-covariance matrix of the predictions

      VarFun <- calc.VarFun(parvar[pars,pars], jacob)


   # * * * * * * * * * * * * #
   # 3) calculate CI of predictions

      # support vector w/o additional design variables
      if (class(xsupport) %in% c('list')) {
        xsupport.x <- xsupport[[xcol]]
      } else {
        xsupport.x <- xsupport
      }

      # Predictions
      yhat <- eval(parse(text=paste('model(x=xsupport,',parstr,')',sep='')))
      pred0 <- data.frame(x=xsupport.x, y=yhat)

      # get variance for predictions from diagonal
      if (n == Inf){ # -> variance is not estimated, do not correct for bias
         varconf <- diag(VarFun)
      } else { # correct for bias of estimated variance
         varconf <- diag(VarFun) * n/df
      }
      sd <- sqrt(varconf)
 
      # calculate CI for predictions
      CI <- cbind(pred0$x,as.data.frame(yhat + as.matrix(sd) %*% t(as.matrix(fac))))
      names(CI) <- c('x','CI.lb','CI.ub')

   # * * * * * * * * * * * * #
   # 3) calculate PI

     varpred <- switch(errmod,
       add =  rep(parval$a^2, length(xsupport.x)),
       prop = (parval$b * xsupport.x)^2,
       comb = (parval$a + parval$b * xsupport.x)^2,
       exp = (parval$e * xsupport.x)^2
     )
    sdpred <- sqrt(varpred+varconf)

      # calculate PI
      PI <- cbind(pred0$x,as.data.frame(yhat + as.matrix(sdpred) %*% t(as.matrix(fac))))
      names(PI) <- c('x','PI.lb','PI.ub')

   ###########################################
   # OUTPUT
     return(list(CImethod='delta',CI.estimates=CI.estimates,pred = pred0,CI=CI,PI=PI))

}
