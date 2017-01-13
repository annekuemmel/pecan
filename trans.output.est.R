trans.output.est <- function(estres, par.trans)
{
   estres.out <- estres

   if (!is.null(par.trans))
   {

      # parameter back transformation
      for (i in 1:length(estres.out$estimates))
         estres.out$estimates[[i]] <- transinv.par(estres.out$estimates[[i]],par.trans[i])

      # SE back transformation
      # relies on assumption that parameters are back-transformed
      for (i in 1:length(estres.out$SE))
         estres.out$SE[[i]] <- transinv.se(estres.out$SE[[i]],estres.out$estimates[[i]],par.trans[i])

   }

   return(estres.out)
}