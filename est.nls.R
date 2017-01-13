est.nls <- function(model, data, init, errmod, ...){

# estimation function without much diagnostics and checking to be used during MC or bootrapping

   # get parameters of the model
   model.args <- names(formals(model))
   pars <- setdiff(model.args,'x')

   # define strings to construct expression to estimate parameters given the initials
   parstr <- paste(pars,collapse=',')
   initstr <- paste(pars,unlist(init[pars]),collapse=',',sep='=')

   # do estimation with nls
   estres <- try(suppressWarnings(eval(parse(text=paste('nls(y ~ model(x,',parstr,'), data = data, start = list(',initstr,'),...)',sep='')))), silent = TRUE)

   # extract results
   if (class(estres) != 'try-error') {
      estimates <- as.list(c(coef(estres),sd(residuals(estres))))  # parameter estimates plus residual error
      names(estimates) <- c(pars,'a')
   } else {
      estimates <- NULL
   }
   
   return(estimates)
}