# function to plot data points / sampling points and predictions with CI interval
plot.data.CI <- function(

   data, # data.frame
   xcol = 'x',
   ycol = 'y',

   est, # estimation results; !obs! parameter confidence taken from pred!

   pred, # list with predictions, confidence and prediction intervals

   par.trans = NULL,

   col.dat = 'navyblue',
   shape.dat = 1,
   col.pred = 'firebrick',
   col.sim = 'firebrick',

   xlab = 'x',
   ylab = 'y',
   log.x = FALSE,
   log.y = FALSE,

   xbreaks = waiver(),
   ybreaks = waiver(),
   xbreaksminor = waiver(),
   ybreaksminor = waiver(),

   CI.level
){

  ################
  # input checks #

   # check whether data is data.frame and has column defined as x and y
   if (!('data.frame' %in% class(data))) data <- as.data.frame(data[c(xcol,ycol)])

   # keep only xcol and ycol and rename with 'x' and 'y'
   data <- data[,c(xcol,ycol)]
   names(data)  <- c('x','y')

   # check whether CI available
   if (!is.null(pred$CI)) flag.CI <- TRUE else flag.CI <- FALSE
   # check whether PI available
   if (!is.null(pred$PI)) flag.PI <- TRUE else flag.PI <- FALSE
   # check whether n available
   if (!is.null(pred$n)) flag.n <- TRUE else flag.n <- FALSE

   # check whether CI.level given and define percentage
   if (!hasArg(CI.level)) cat('PLOT: confidence level not given.\n')
   CI.perc = CI.level*100
   
  #############
  # functions #
  # get required packages
  require(ggplot2)
  theme_set(theme_bw())
  theme_update(legend.key = element_blank(), panel.background  = element_blank(), plot.background = element_blank())
  require(grid)
  require(gridExtra)
  # auxiliary functions
  number_ticks <- function(n) {function(limits) pretty(limits, n)}

   ##############
   # build plot #

   # define color scale
   color.def <- c(col.dat, col.pred, col.pred, col.sim)
   # define linetype scale
   linetype.def <- c(0, 1, 2, 3)
   # define shape scale
   shape.def <- c(shape.dat, NA, NA, NA)
   # name scales
   names(linetype.def) <- names(color.def) <- names(shape.def) <- c('Observations', 'Prediction', paste(CI.perc,'% CI',sep=''), paste(CI.perc,'% PI',sep=''))

   # create data frame for plotting
   df.plot <- rbind(within(data,{what='Observations';which='a'}),
                    within(pred$pred, {what='Prediction';which='b'}))
   whatlevels <- c('Observations','Prediction') # used to get legend in proper order
   if (flag.CI){
      df.plot <- rbind(df.plot,
                    within(pred$CI, {y=CI.lb;what=paste(CI.perc,'% CI',sep='');which='c'})[,c('x','y','which','what')],
                    within(pred$CI, {y=CI.ub;what=paste(CI.perc,'% CI',sep='');which='d'})[,c('x','y','which','what')])
         whatlevels <- c(whatlevels,paste(CI.perc,'% CI',sep=''))
   }
   if (flag.PI){
      df.plot <- rbind(df.plot,
                    within(pred$PI, {y=PI.lb;what=paste(CI.perc,'% PI',sep='');which='e'})[,c('x','y','which','what')],
                    within(pred$PI, {y=PI.ub;what=paste(CI.perc,'% PI',sep='');which='f'})[,c('x','y','which','what')])
         whatlevels <- c(whatlevels,paste(CI.perc,'% PI',sep=''))
   }
   df.plot$what <- factor(df.plot$what, levels = whatlevels)

   # build plot
   gr <- ggplot(data=subset(df.plot)) +
                geom_line(aes(x=x,y=y, linetype = what, color = what, group = which), na.rm=TRUE, size=0.8) +
                geom_point(aes(x=x,y=y, color = what, shape = what), na.rm=TRUE) +
                scale_color_manual('', values=color.def) +
                scale_linetype_manual('', values=linetype.def) +
                scale_shape_manual('', values=shape.def) +
                labs(x=xlab,y=ylab)

   # do logarithmic axis if required
   logstr <- ''
   if (log.x) {
      gr <- gr + scale_x_log10(breaks=xbreaks, minor_breaks=xbreaksminor)
      logstr <- paste(logstr,'b',sep='')
   } else {
      if (class(xbreaks) == 'waiver') xbreaks = number_ticks(n=6)
      gr <- gr + scale_x_continuous(breaks=xbreaks, minor_breaks=xbreaksminor)
   }
   if (log.y) {
      gr <- gr + scale_y_log10(breaks=ybreaks, minor_breaks=ybreaksminor)
      logstr <- paste(logstr,'l',sep='')
   } else {
      if (class(xbreaks) == 'waiver') ybreaks = number_ticks(n=6)
      gr <- gr + scale_y_continuous(breaks=ybreaks)
   }
   gr <- gr + annotation_logticks(sides=logstr)


   # Info for iterative procedures
    if (flag.n){
     build <- ggplot_build(gr)
     get.textcoord <- function(range, frac, log){
      if (log) {
         coord <- 10^(range[1] + frac*diff(range))
      } else {
         coord <- range[1] + frac*diff(range)
      }
     }
     info <- data.frame(x=get.textcoord(build$panel$ranges[[1]]$x.range, 0.05, log.x),
                        y=get.textcoord(build$panel$ranges[[1]]$y.range, 0.05, log.y),
                        label=paste(pred$n.fail,' out of ',pred$n,' runs failed.',sep=''))
     gr <- gr +  geom_text(data=info, aes(x=x,y=y,label=label), color = 'black', size = 3, hjust=0,vjust=0)
    }


   # Parameter estimate display

     # information on parameter transformation
     if (is.null(par.trans))
     {
       par.trans <- rep("",length(est$estimates))
       names(par.trans) <- names(est$estimates)
     } else {
       par.trans <- gsub(" (c)","",paste(" (",par.trans,")",sep=""), fixed=TRUE)
       names(par.trans) <- names(est$estimates)
       if (est$estmethod == "nls") par.trans <- c(par.trans,"a"="") # residual error parameter for nls
     }

     # construction of parameter estimate table
       # confidence interval not calculated for all error parameters depending on estimation method and CI calculation method used
       # need to fill up output table here

       CI <- paste(
          ' (',
          formatC(pred$CI.estimates$CI.lb[pred$CI.estimates$parameter %in% names(est$estimates)]),' - ',
          formatC(pred$CI.estimates$CI.ub[pred$CI.estimates$parameter %in% names(est$estimates)]),')',
          sep=''
       )
       if (length(pred$CI.estimates$parameter) < length(est$estimates))
       {
         ndiff <- length(est$estimates) - length(pred$CI.estimates$parameter)
         CI <- c(CI,rep("-",ndiff))
       }
       resulttable <- data.frame(
         Parameter = paste(names(est$estimates),par.trans, sep=""),
         Estimate  = formatC(est$estimates),
         CI = CI
       )
       names(resulttable) <- gsub("CI",paste(CI.perc,"% CI",sep=""),names(resulttable))

       # build grob
       gr.resinfo <- tableGrob(
         resulttable,
         gp = gpar(cex=0.75), padding.v = unit(3, "mm"),
         show.rownames = FALSE
       )
       h1 <- grobHeight(gr.resinfo)
       gr.resinfo$vp <- viewport(x=0, y=unit(1,"npc")-h1*0.5-unit(1,"line"), just=c(0,0.5))

   # Constrauction of method information display
   # error model information
     switch(est$errmod,
       add  = errmodstr <- 'additive (a)',
       prop = errmodstr <- 'proportional (b)',
       comb = errmodstr <- 'additive + proportional (a + b)',
       exp  = errmodstr <- 'exponential (e)' )
   
     if (est$estmethod == 'nlm') estmeth <- 'MLE' else estmeth <- est$estmethod

     methinfo <- data.frame(
       Setting = c("Estimation method", "CI calculation method", "Error model"),
       Value = c(estmeth, pred$CImethod, errmodstr)
     )
     gr.methinfo <- tableGrob(
       methinfo,
       gp = gpar(cex=0.75), padding.v = unit(3, "mm"),
       show.rownames = FALSE, show.colnames = TRUE
     )
     h2 <- grobHeight(gr.methinfo)
     gr.methinfo$vp <- viewport(x=0.2,y=unit(1,"npc")-h2*0.5-unit(1,"line"), just=c(0,0.5), clip = "off")

   # build graphical output
   gr.info <- arrangeGrob(gr.methinfo,gr.resinfo, ncol=2, nrow = 1, widths = c(1,1.5))
   gr.all <- arrangeGrob(gr.info,gr, ncol=1, nrow = 2, heights = c(1,2.5))

   return(gr.all)
}



  
