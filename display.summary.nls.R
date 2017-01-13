#           function to summarize NLS estimation
# using the output list from parameter estimation function

   display.summary.nls <- function(OUT){
     with(OUT,{

       cat('\nNonlinear least-squares estimation: y ~ model(x,',paste(setdiff(names(estimates),c('a','b','e')),collapse=','),')\n',sep="")

       cat('\nEstimates:\n')
       resmat <- matrix(c(estimates, SE, SE/estimates*100),ncol=3)
       dimnames(resmat) <- list(names(estimates),c('Estimate', 'SE','RSE (%)'))
       print(resmat, digits = 4)

       cat('\nCorrelation of estimates:\n')
       estcor <- cov2cor(estvar)
       estcor[upper.tri(estcor,diag=TRUE)] <- NA
       estcor <- estcor[-1,-dim(estcor)[1]]
       print(estcor, digits = 3, na.print="")

       DF <- df.residual(estimation)
       resSE <- sqrt(sum(residuals(estimation)^2)/DF)
       cat('\nResidual standard error:',resSE,'on',DF,'degrees of freedom\n')

       if (estimation$convInfo$isConv) code.term <- '' else code.term <- 'NOT '

       cat('\nConverence ',code.term,'achieved.\n',sep='')
       cat('Number of iterations:',estimation$convInfo$finIter,'\n')
       cat('Achieved convergence tolerance:',estimation$convInfo$finTol,'\n')
     })
   }