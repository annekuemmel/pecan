#           function to summarize MLE estimation
# using the output list from parameter estimation function

   display.summary.nlm <- function(OUT){
     with(OUT,{

       switch(errmod,
          add = errmod.name <- 'additive',
          prop = errmod.name <- 'proportional',
          comb = errmod.name <- 'combined (additive + proportional)',
          exp = errmod.name <- 'exponential'
       )
       cat('\nLog-likelihood maximization: y ~ model(x,',paste(setdiff(names(estimates),c('a','b','e')),collapse=','),') with ',errmod.name,' error\n',sep="")

       cat('\nEstimates:\n')
       resmat <- matrix(c(estimates, SE, SE/estimates*100),ncol=3)
       dimnames(resmat) <- list(names(estimates),c('Estimate', 'SE','RSE (%)'))
       print(resmat, digits = 4)

       cat('\nCorrelation of estimates:\n')
       estcor <- cov2cor(estvar)
       estcor[upper.tri(estcor,diag=TRUE)] <- NA
       estcor <- estcor[-1,-dim(estcor)[1]]
       print(estcor, digits = 3, na.print="")

       cat('\nNumber of iterations:',estimation$iterations,'\n')
       switch(estimation$code,
              "1"=code.term <- 'Relative gradient is close to zero. OK!',
              "2"=code.term <- 'Successive iterates within tolerance. OK!',
              "3"=code.term <- 'Last global step failed to locate a point lower than estimate. OK? steptol too small?',
              "4"=code.term <- 'Iteration limit exceeded.',
              "5"=code.term <- 'Maximum step size stepmax exceeded five consecutive times.\n    Function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small?'
       )
       cat('\nOptimization termination reason:\n   ',code.term,'\n')

       cat('\nEstimated maximum likelihood (-2LL) :', estimation$neg2LL)
       cat('\nBayesian Information Criterion (BIC):', estimation$BIC,'\n')
     })
   }