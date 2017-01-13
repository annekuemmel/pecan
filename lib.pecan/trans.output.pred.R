trans.output.pred <- function(pred, par.trans)
{
   pred.out <- pred

   if (!is.null(par.trans))
   {

      # parameter back transformation
      for (i in 1:dim(pred.out$CI.estimates)[1])
      {
         pred.out$CI.estimates[i,"CI.lb"] <- transinv.par(pred.out$CI.estimates[i,"CI.lb"],par.trans[i])
         pred.out$CI.estimates[i,"CI.ub"] <- transinv.par(pred.out$CI.estimates[i,"CI.ub"],par.trans[i])
      }
   }

   return(pred.out)
}


