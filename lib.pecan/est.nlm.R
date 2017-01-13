est.nlm <- function(model, data, init, errmod, ...){

# estimation function without much diagnostics and checking to be used during MC or bootrapping

   # get model parameters
   pars <- setdiff(names(formals(model)),'x')

   # number of parameter w/o error model
   nphi <- length(pars)

   # model reformulation
   modelTP <- function(t,p){
       pars <- pars
       df <- data.frame(pars=pars,no=paste('p[',1:length(pars),']',sep=''))
       transstr <- paste(df$pars,df$no, sep = '=', collapse = ',')
       eval(parse(text=paste('y <- model(x=t,',transstr,' )',sep='')))
       return(y)
   }

   # define score function (-2LL) depending on error model
   if (errmod == 'add'){
      fmin <- function(p,y,t){
         nphi <- length(p)-1
         yhat <- modelTP(t, p)
         g <- p[nphi+1]
         e <- sum( ((y-yhat)/g)^2 + log(g^2) )
         return(e)
      }
      start <- c(unlist(init))
   }
   if (errmod == 'prop') {
      fmin <- function(p,y,t){
         nphi <- length(p)-1
         yhat <- modelTP(t, p)
         g <- p[nphi+1]*yhat
         e <- sum( ((y-yhat)/g)^2 + log(g^2) )
         return(e)
      }
      start <- c(unlist(init))
   }
   if (errmod == 'comb') {
      fmin <- function(p,y,t){
         nphi <- length(p)-2
         yhat <- modelTP(t, p)
         g <- abs(p[nphi+1]) + abs(p[nphi+2])*yhat
         e <- sum( ((y-yhat)/g)^2 + log(g^2) )
         return(e)
      }
      start <- c(unlist(init))
   }
   if (errmod == 'exp') {
      fmin <- function(p,y,t){
         nphi <- length(p)-1
         yhat <- modelTP(t, p)
         g <- p[nphi+1]
         e <- sum( ((log(y)-log(yhat))/g)^2 + log(g^2) )
         return(e)
      }
      start <- c(unlist(init))
   }
   if (!exists('fmin')) stop('Score function for MLE not defined. Given error model unknown.')
   npsi <- length(start) # number of parameters including error model


   # do estimation
   estres <- try(suppressWarnings(nlm(fmin, start, y=data$y,t=data$x,...)), silent = TRUE)


   # extract results
   if (class(estres) != 'try-error') {

      # make sure that estimates for residual error standard deviation are positive
      if (errmod == 'comb') estres$estimate[npsi-c(1,0)] <- abs(estres$estimate[npsi-c(1,0)])
      if (errmod != 'comb') estres$estimate[npsi] <- abs(estres$estimate[npsi])

      estimates <- as.list(estres$estimate)  # parameter estimates plus residual error
      switch(errmod,
         add  = names(estimates) <- c(pars,'a'),
         prop = names(estimates) <- c(pars,'b'),
         comb = names(estimates) <- c(pars,'a','b'),
         exp  = names(estimates) <- c(pars,'e')
      )
   } else {
      estimates <- NULL
   }
   
   return(estimates)
}
