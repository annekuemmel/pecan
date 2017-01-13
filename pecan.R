pecan <- function(

  # Required input arguments are data, init and model.
  #  (However, running without any of the required inputs will run an example.)
  # Model parameter names a, b and e not allowed. They are reserved for error model parameters.
  # When using ADVAN models, 
  #   - the design input argument has to be a list with either
  #      o 'dostime' and 'dose' specifying oral or iv dosing or
  #      o 'infstart','infend','rate' specifying an infusion
  #   - the xcol has to be "time"

  # OUTPUT:
  #  list with following entities:
  #    estres: list with information on and results for the parameter estimation
  #    pred:   list with information on and results for the confidence and prediction intervals
  #    gr:     graph displaying predictions, confidence overlaid with data and inforamtion on the estimation
  #    diagnostics:  list of diagnostic plots (observed vs. predicted and weighted residual plots)

   data = Indometh,       # (required) named data frame
   xcol = 'time',         # column name of independent variable
   ycol = 'conc',         # column name of dependent variable
   design = NULL,         # additional design variables needed in model (e.g. dosing specification)

   init = list(ke=1,V=0.5,k12=0.5,k21=0.5),     # (required) named list of initial parameter guesses, example: init = list(V=1,ke=3)
   par.trans = NULL,                            # parameter transformation, either NULL, single string or character vector with length(init).
                                                # Allowed strings: "log", "logit", "c"

   model = function(x,ke,V,k12,k21){
     s1 <- k12+k21+ke
     s2 <- k21*ke
     beta <- 0.5*(s1-sqrt(s1^2 - 4*s2))
     alpha <- s2 / beta
     A <- 1/V * (alpha-k21)/(alpha-beta)
     B <- 1/V * (beta-k21)/(beta-alpha)

     y <- A*exp(-alpha*x) + B*exp(-beta*x)
     return(y)
   },                  # (required) model function with 'x' as first input argument
                       # followed by the model parameters: function(x, p1, p2, p3 , ...) 

   model.deriv = NULL, # optional function computing the derivatives of the model wrt the parameters,
                       # if not given numerical derivaties used
   errmod = 'exp',     # error model: 'add', 'prop', 'comb', 'exp'; when using 'nls' automatically add is assumed
   init.err = NULL,    # initial estimates for error model parameters,
                       #    specify as numeric vector,
                       #    order for combined error model: first additive then proportional term
   estmethod = 'nlm',  # estimation method: 'nls' (least-squares) or 'nlm' (MLE)

   CI.level = 0.9,     # confidence/coverage interval
   CImethod = 'sim', # method to determine confidence intervals, 'delta', 'sim', 'boot', 'mc'
   xsupport = NULL,    # vector of x-values at which to calculate predictions and confidence intervals
   nsmpl = 100,       # number of samples for simulation/estimation methods

   # graphic settings:
   log.x = FALSE,      # logarithmic x-scale?
   log.y = TRUE,       # Note: Only for plotting, original scale used for estimation

   col.dat='steelblue',   # color for observations
   shape.dat=1,           # shape for observations
   col.pred='olivedrab',  # color for confidence interval
   col.sim='olivedrab',   # color for prediction/coverage interval
   xlab='time',ylab='concentration',     # axis labels

   xbreaks=waiver(), ybreaks=waiver(),           # breaks/major grid for axes
   xbreaksminor=waiver(), ybreaksminor=waiver(), # minor grid for axes

   ... # additional arguments passed to estimation functions (nls and nlm),
       # e.g., algorithms or tolerance levels.

   ){





   ###########################################
   # Input checks

   # check whether all or none example inputs are given ("required inputs")
   ex_input_var <- c("data","init","model")
   ex_input_var_flag <- c()
   for (iex in ex_input_var) ex_input_var_flag[iex] <- eval(parse(text=paste("hasArg(",iex,")",sep="")))
   if (!(all(ex_input_var_flag) | all(!ex_input_var_flag))) stop("Give all required inputs (data, model, init) or none to run example.\n", call.=FALSE)

   # check whether model is function and exists
   if (!is.function(model)) stop('Oops, given model does not exist or is no function.\n', .call=FALSE)

   # check whether first argument in independent variable 'x'
   model.args <- names(formals(model))
   if (model.args[1] != 'x') stop('First input argument of model function needs to be "x".\n', call.=FALSE)

   # check whether given initials are according to parameter input argument of model function
   pars <- setdiff(model.args,'x')
   init.pars <- names(init)
   if (!setequal(pars,init.pars)) stop("Oops, given initials not corresponding to model parameters.\n", call.=FALSE)

   # check whether model has not-allowed parameter names
   parnames_nogo <- c('a','b','e')
   if (any(parnames_nogo %in% pars)) stop('Please do not use ',paste(intersect(parnames_nogo,pars),collapse=', '),' as model parameter names.\n', call.=FALSE)

   # check whether data is data frame
   if (!is.data.frame(data)) stop('Please give data as data frame.\n', call.=FALSE)

   # check whether data has defined columns defined
   if (!all(c(xcol,ycol) %in% names(data))) stop('Non-existing colums defined as x and/or y.\n', call.=FALSE)

   # check missing values in xcol and ycol
   miss.rows <- apply(is.na(data[,c(xcol,ycol)]),1,any)
   if ( any(miss.rows) )
   {
      cat('Remove', sum(miss.rows), 'observations due to missingness.\n')
      data <- data[!miss.rows,]
   }

   # check whether CI level is between 0 and 1
   if (CI.level > 1 | CI.level < 0) stop('Invalid CI level (not between 0 and 1).\n', call.=FALSE)

   # check whether number of sampling for estimation/simulations meaningful
   if (CImethod %in% c('sim','mc','boot'))
   {
      if (nsmpl <= 1)  stop('\nNumber of replicates for simulations/estimations cannot be smaller than 2.\n',call.=FALSE)
      if (nsmpl <= 20) cat('\n!!! Number of replicates for simulations/estimations quite small. Results may not be reliable.\n')
   }

   # check on parameter transformation
   if (!is.null(par.trans))
   {
      # check whether only allowed transforamtions given
      if (!all(par.trans %in% c("c","log"))) stop("\n Unknown parameter transformation given. Currently, only log or non-tranformed allowed. \n")
      # if only one transformation is given, assume that it applies to all
      cat("\n Unique parameter transformation given applied to all parameters !!! \n")
      if (length(par.trans) == 1) par.trans <- rep(par.trans, length(init))
      else
      {
         if (length(par.trans != length(init))) stop("\n Number of parameter transformation not equal to number of parameters! \n")
      }
      # give names (OBS! this means that input order matters)
      names(par.trans) <- names(init)

      # add residual parameter estimates
      switch(errmod,
        add = par.trans <- c(par.trans,"a"="c"),
        comb = par.trans <- c(par.trans,"a"="c","b"="c"),
        prop = par.trans <- c(par.trans,"b"="c"),
        exp = par.trans <- c(par.trans,"e"="c"),
      )
   }

   ##########################
   #      PREPARATIONS

   # get number of observations
   n <-  length(data[[xcol]])

   # check whether independent variable vector used for predictions exist
   # if not, create one:
   # support vector has equal steps on scale used for plotting
   if (is.null(xsupport)) 
   {
      if (log.x) {
         # check whether zero or negative x values exist
         posX <- data[[xcol]] > 0
         if (all(posX)) {
            # all positive, support vector covering full x range
            xsupport <- 10^seq(min(log10(data[[xcol]])),max(log10(data[[xcol]])),length.out = 20)
         } else {
            # support vector to cover positive range
            # (starting one log10 lower than smallest positive value)
            xsupport <- 10^seq(min(log10(data[[xcol]][posX])-1),max(log10(data[[xcol]][posX])),length.out = 20)
         }
      } else {
         xsupport <- seq(min(data[[xcol]]),max(data[[xcol]]),length.out = 20)
      }
   }

   # Reformat data to have design columns within x:
   # The data and pred internal variables are lists with x and y entries
   # The y entry is always the numeric vector "ycol" 
   # The x entry is the numeric vector "xcol" if no design variables given/required and
   #                a list with the xcol and the design variables with xcol as first entry
   #                (collection of diagnostics relies when estimating parameters relies on this!!!)
   if (!is.null(design)) {

   # design information needed in model -> pack design information to x in internal data container
   cat('\nDesign information found. Make sure that given design variable is properly handled by model function.\n')
     dataint <- list(x=c(list(data[[xcol]]),design), y=data[[ycol]])
     names(dataint$x) <- c(xcol,names(design))
     predint <- dataint; predint$x[[xcol]] <- xsupport; predint$y <- NULL }
   else {

   # assume no further information needed and only regressor is the x vector
     dataint <- list(x=data[[xcol]],y=data[[ycol]])
     predint <- list(x=xsupport) }

   # store original values
   model.orig <- model
   init.orig <- init

   # Parameter transformation if required
   if (!is.null(par.trans)) {

     # transform parameter initial values
     for (i in 1:length(init))
       init[[i]] <- trans.par(init[[i]],par.trans[i])

     # transform model inputs
     df <- data.frame(trans=NA,par=pars, stringsAsFactors=FALSE)
     for (i in 1:length(init)) # ! OPS ! assumes that residual parameters are at the end....
       df$trans[i] <- trans.mod(par.trans[i])
     strmodtrans <- paste(df$trans,df$par,sep="(", collapse = "), ")

     eval(parse(text=paste(
       "model <- function( x,",
       paste(pars,collapse=", ",sep=""),
       ") { y <- model.orig(x, ",
       strmodtrans,
       ") ); return(y) }",
     sep="")))
   }



   ##########################
   #       ESTIMATION
   # estimate parameters using either NLS or MLE

   flagestmet <- FALSE   
   if (estmethod == 'nls') {

      flagestmet <-TRUE
      if (errmod != 'add') cat('\nWhen estimating parameters using NLS, error model is additive!\n')
      errmod = 'add'
      estres <- parameter.estimation.nls(
          model = model,
          init=init, data = dataint,
          ...
          )

      # back transformation of parameters for output
      estres.out <- trans.output.est(estres, par.trans)
      display.summary.nls(estres.out)

   }
   if (estmethod == 'nlm') {

      flagestmet <- TRUE
      estres <- parameter.estimation.nlm(
          model = model, errmod=errmod,
          init=init, data = dataint,
          init.err = init.err,
          ...
          )

      # back transformation of parameters for output
      estres.out <- trans.output.est(estres, par.trans)
      display.summary.nlm(estres.out)

   }
   if (!flagestmet) stop('Unknown estimation method.\n', call.=FALSE)


   #########################
   #      DIAGNOSTICS
   # generates diagnostic plots
   # (observed/predicted, residual histogram,
   # residuals against predicted and residuals against independent variable)

   gr_diagn <- plot.diagnostics(estres$diagn,log.x,log.y)


   #########################
   #       CONFIDENCE
   # determine confidences for parameters and predictions using method of choice
   flagCImet <- FALSE

   # delta method
   if (CImethod == 'delta'){
       flagCImet <- TRUE
       pred <- CI.delta(
          model = model,
          model.deriv = model.deriv,
          errmod = errmod,
          parval = estres$estimates,
          parvar = estres$estvar,
          xsupport = predint$x,
          n=n,
          CI.level=CI.level,
          xcol=xcol
       )
   }

   # based on simulations
   # sampling parameters from their "error distribution"
   if (CImethod == 'sim'){
       flagCImet <- TRUE
       pred <- CI.sim(
          model = model,
          errmod = errmod,
          errval = estres$errval,
          parval = estres$estimates,
          parvar = estres$estvar,
          xsupport = predint$x,
          n=n,
          xcol=xcol, 
          CI.level=CI.level,
          nsim = nsmpl
       )
   }

   # based on MC simulations of datasets
   # aka "parametric bootstrap"
   if (CImethod == 'mc'){
       flagCImet <- TRUE
       pred <- CI.mc(
          model = model,
          errmod = errmod,
          errval = estres$errval,
          parval = estres$estimates,
          xsample = dataint$x,
          xsupport = predint$x,
          xcol = xcol,
          CI.level=CI.level,
          nmc = nsmpl,
          estmethod = estmethod,
          ...
       )
   }

   # based on bootstrapping the data
   if (CImethod == 'boot'){
       flagCImet <- TRUE
       pred <- CI.boot(
          model = model,
          errmod = errmod,
          init = estres$estimates,
          data = dataint,
          xcol=xcol,
          xsupport = predint$x,
          CI.level=CI.level,
          nboot = nsmpl,
          estmethod = estmethod,
          ...
       )
   }
   if (!flagCImet) {cat('\nSorry, unknown method to determine confidences!\n'); pred=NA }


   # potential back-transformation of CI limits of parameter estimates
   pred.out <- trans.output.pred(pred, par.trans)


   #########################
   #        GRAPH
   # Predictions and intervals overlaid with the data
   # and parameter estimates with some information on the calculations
   gr <- plot.data.CI(
           data=data,
           xcol=xcol,ycol=ycol,
           est=estres.out,
           pred=pred.out,
           par.trans=par.trans,
           CI.level=CI.level,
           log.x=log.x,log.y=log.y,
           col.dat=col.dat,col.pred=col.pred,col.sim=col.sim,
           xlab=xlab,ylab=ylab,
           xbreaks = xbreaks, ybreaks = ybreaks,
           xbreaksminor = xbreaksminor, ybreaksminor = ybreaksminor,
         )


   #########################
   #        OUTPUT
   return(
     list(
       estres=estres.out,
       pred=pred.out,
       gr=gr,
       diagnostics=gr_diagn
     )
   )

}