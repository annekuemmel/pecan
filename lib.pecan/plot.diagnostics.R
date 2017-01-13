plot.diagnostics <- function(diagn,log.x,log.y){

   # Observed versus predicted
   yrange <- range(c(diagn$y,diagn$ypred))
   gr_obsVpred <- ggplot(diagn, aes(ypred,y)) +
                     geom_abline(color='firebrick') +
                     geom_point(shape=1) +
                     geom_smooth(color='steelblue', se = FALSE, method = 'loess',size=0.8) +
                     coord_fixed(xlim=yrange,ylim=yrange) +
                     labs(x='Predictions',y='Observations')
   if (log.y) gr_obsVpred <- gr_obsVpred + scale_x_log10() + scale_y_log10()

   # Residuals against predictions
   gr_res_ypred <- ggplot(diagn, aes(ypred,wres)) +
                     geom_hline(color='firebrick') +
                     geom_point(shape=1) +
                     geom_smooth(color='steelblue', se = FALSE, method = 'loess',size=0.8) +
                     labs(x='Predictions',y='Weighted residuals')
   if (log.y) gr_res_ypred <- gr_res_ypred + scale_x_log10()

   # Residuals against independent variable (time)
   diagn$cens <- FALSE
   if (log.x) 
   {
     diagn$cens[diagn$t<=0] <- TRUE
     diagn$t[diagn$cens] <- exp(min(log10(diagn$t[!diagn$cens]))-0.05*diff(range(log10(diagn$t[!diagn$cens]))))
   }
   gr_res_t <- ggplot(diagn, aes(t,wres)) +
                     geom_hline(color='firebrick') +
                     geom_point(shape=21,aes(fill=cens)) +
                     geom_smooth(color='steelblue', se = FALSE, method = 'loess',size=0.8) +
                     labs(x='Independent variable',y='Weighted residuals') +
                     scale_fill_manual(values=c("transparent","tomato"),guide=FALSE)
   if (log.x) gr_res_t <- gr_res_t + scale_x_log10()

   # Histogram of residuals
   gr_res_hist <- ggplot(diagn,aes(wres)) +
                     geom_histogram(aes(y = ..density..),binwidth = diff(range(diagn$wres))/30,
                                    fill='grey',color='black') +
                     geom_density(color='steelblue', size=0.8) +
                     stat_function(fun = dnorm, colour = "red") +
                     coord_cartesian(xlim=c(-1.1,1.1)*max(abs(diagn$wres))) +
                     labs(x='Weighted residuals', y='Density')

   OUT <- list(obsVpred=gr_obsVpred,res_ypred=gr_res_ypred,res_t=gr_res_t,res_hist=gr_res_hist)


return(OUT)

}
