# # # # # # # # # # # # # # #
#      PREPARATIONS

# libraries
require(ggplot2)
theme_set(theme_bw())

# get pecan calculation functions
for (ifi in list.files('./lib.pecan', pattern=".R$")) source(file.path('./lib.pecan',ifi))
rm(ifi)


# get model functions
for (ifi in list.files('./lib.models', pattern=".R$")) source(file.path('./lib.models',ifi))
rm(ifi)


# # # # # # # # # # # # # # # # # # # # #
#              EXAMPLE RUN

# data:
data = Indometh
gr_data <- ggplot(data, aes(time,conc)) + geom_point() + scale_y_log10()
gr_data

model = function(x,ke,V,k12,k21){
     s1 <- k12+k21+ke
     s2 <- k21*ke
     beta <- 0.5*(s1-sqrt(s1^2 - 4*s2))
     alpha <- s2 / beta
     A <- 1/V * (alpha-k21)/(alpha-beta)
     B <- 1/V * (beta-k21)/(beta-alpha)

     y <- A*exp(-alpha*x) + B*exp(-beta*x)
     return(y)
   }

init = list(ke=1,V=0.5,k12=0.5,k21=0.5)

gr_data_init <- gr_data + stat_function(fun=model, args=init)
gr_data_init

# "pecan" estimates model parameters for regression models (naive pool)
# Example with Indometh dataset and a two compartment model:
res.example <- pecan(errmod="exp")

# Results are displayed in a graph:
grid.arrange(res.example$gr)


# Diagnostic plots for the estimation:
res.example$diagnostics[[1]] # observations/predictions
res.example$diagnostics[[2]] # weighted residuals against predictions
res.example$diagnostics[[3]] # weighted residuals against independent variable
res.example$diagnostics[[4]] # histogram of weighted residuals


# Results stored in other list elements:
# - estres: parameter estimation information and results
# - pred: confidence intervals for parameters and confidence and prediction intervals for dependent variable



# # # # # # # # # # # # # # # # # # # # #
#         EXAMPLE: MLE / NLS

# example for emax model with simulated data

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# simulate AUC-response data
  set.seed(770718)

# set doses
  doses <- c(0, 3, 10, 30, 100, 500)

# set parameters for Emax model
  Emax = 100   # maximum effect
  E0 = 10      # baseline of effect
  EC50 = 45    # AUC resulting in 50% of maximum effect
  a = 5        # additive residual error on effect
  b = 0.2      # proportional residual error on effect
  iivCL = 0.3  # inter-individual variability of clearance
  dose = 1     # dose

# simulate data
  data.emax <- within(
    data.frame(dose=rep(doses,each=8)),  # simulate 8 observations per dose
    {
      auc <- dose * exp(rnorm(length(dose),sd=iivCL))  # simulate AUCs
      effect <- emax(x=auc, E0=E0, Emax=Emax, EC50 = EC50)  # predict effect
      effect <- effect + (a+effect*b)*rnorm(length(dose))   # add residual error
    }
  )

# set initials
  init.emax <- list(Emax=90, E0 = 20, EC50 = 100)

# plot simulated data:
  ggplot(data.emax, aes(auc,effect)) + geom_point() + scale_x_log10()


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# estimate parameters using NLS and CI by simulation
res.emax.nls <- pecan(
   
   # required parameters
   data = data.emax,
   init = init.emax,
   model = emax,  # chosen from model library

   # define x and y columns
   xcol = "auc",
   ycol = "effect",

   # choose estimation method
   estmethod = 'nls',

   # axis settings for plotting
   log.x = TRUE,
   log.y = FALSE,
   xlab = "AUC",
   ylab = "effect"

)
grid.arrange(res.emax.nls$gr)
res.emax.nls$diagnostics[[2]]


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# estimate parameters using maximum likelihood and CI by simulation
res.emax.nlm <- pecan(
   
   # required parameters
   data = data.emax,
   init = init.emax,
   model = emax,

   # error model
   errmod = "comb",
   init.err = c(10,1),

   # define x and y columns
   xcol = "auc",
   ycol = "effect",

   # choose estimation method
   estmethod = 'nlm',

   # choose CI calculation method
   CImethod = 'sim',

   # axis settings for plotting
   log.x = TRUE,
   log.y = FALSE,
   xlab = "AUC",
   ylab = "effect"
)
grid.arrange(res.emax.nlm$gr)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#         EXAMPLE: DELTA / SIM / BOOT / MC



# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# estimate parameters using maximum likelihood and CI by delta method
res.emax.nlm.delta <- pecan(
   
   # required parameters:
   data = data.emax,
   init = init.emax,
   model = emax,

   # error model
   errmod = "comb",
   init.err = c(10,1),

   # define x and y columns:
   xcol = "auc",
   ycol = "effect",

   # choose estimation method:
   estmethod = 'nlm',

   # choose method to define confidence and prediction interval:
   CImethod = "delta",

   # axis settings for plotting:
   log.x = TRUE,
   log.y = FALSE,
   xlab = "AUC",
   ylab = "effect"

)
grid.arrange(res.emax.nlm.delta$gr)




# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# estimate parameters using maximum likelihood and CI by MC
res.emax.nlm.mc <- pecan(
   
   # required parameters:
   data = data.emax,
   init = init.emax,
   model = emax,

   # error model
   errmod = "comb",
   init.err = c(10,1),

   # define x and y columns:
   xcol = "auc",
   ycol = "effect",

   # choose estimation method:
   estmethod = 'nlm',

   # choose method to define confidence and prediction interval:
   CImethod = "mc",
   nsmpl = 500,

   # axis settings for plotting:
   log.x = TRUE,
   log.y = FALSE,
   xlab = "AUC",
   ylab = "effect"

)
grid.arrange(res.emax.nlm.mc$gr)




# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# estimate parameters using maximum likelihood and CI by bootstrapping
res.emax.nlm.boot <- pecan(
   
   # required parameters:
   data = subset(data.emax, sample(c(TRUE,FALSE),dim(data)[1],replace=TRUE,prob=c(0.3,0.7))),
   init = init.emax,
   model = emax,

   # error model
   errmod = "comb",
   init.err = c(10,1),

   # define x and y columns:
   xcol = "auc",
   ycol = "effect",

   # choose estimation method:
   estmethod = 'nlm',

   # choose method to define confidence and prediction interval:
   CImethod = "boot",
   nsmpl = 500,

   # axis settings for plotting:
   log.x = TRUE,
   log.y = FALSE,
   xlab = "AUC",
   ylab = "effect"

)
grid.arrange(res.emax.nlm.boot$gr)



# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #
# The follwing code needs ADVAN-style implementation of PK models
# (only run if folder exists that should contain these functions)
# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

if (file.exists("lib.advan"))
{

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#     EXAMPLE: PK MODEL USING ADVAN IMPLEMENTATION
# data for oral dosing o.d. for 7 days is simulated with
# sampling trough values and profile after last dose

# get Advan library
source('./lib.advan/Wrap_ADVANfunctions.R')

  # When using ADVAN models, 
  #   - the design input argument has to be a list with either
  #      o 'dostime' and 'dose' specifying oral or iv dosing or
  #      o 'infstart','infend','rate' specifying an infusion
  #   - the xcol has to be "time"
  #   - for oral dosing F1 (bioavailability) is set to 1, hence, volumes are apparent



# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# simulate data:
# 2-cpt model with multiple oral doses
set.seed(770719)
# model parameters
CL = 0.1
V2 = 5
Q = 0.1
V3 = 10
KA = 0.6
F1 = 1

# define dosing
dosing = list(dostime=24*0:14,dose=rep(1,length(0:14)))

# prediction
data <- data.frame(time=rep(c((24*1:14)-0.01,24*14+c(0.1,0.25,0.5,1,2,4,8,12,24,36,48,72,96,120)),each=3))
data$conc <- TwoCompOral.fun(c(data,dosing),CL = CL, V2 = V2, Q = Q, V3 = V3, KA = KA)

# add exponential residual error
data$conc <- data$conc * rlnorm(length(data$conc),sd = 0.1)

ggplot(data, aes(time, conc)) + geom_point()


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# parameter estimation using MLE and delta method for confidences
res.advan <- pecan(

   data = data,               # data frame
   xcol = 'time',             # column name of independent variable
   ycol = 'conc',             # column name of dependent variable
   design = dosing,           # list with design definition (user needs to make sure that requirement by ADVAN function is satisfied)
   init = list(CL=1,V2=1,V3=1,Q=1,KA=1),   # named list of initial parameter guesses, example: init = list(V=1,ke=3)

   model = TwoCompOral.fun,

   errmod = 'exp',     # error model: 'add', 'prop', 'comb', 'exp'; when using 'nls' automatically add is assumed
   estmethod = 'nlm',  # estimation method: 'nls' (least-squares) or 'nlm' (MLE)
   init.err = c(0.2),

   CImethod = 'delta',            # method to determine confidence intervals, 'delta', 'sim', 'boot', 'mc'
   xsupport = c(0:14*24,14*24+seq(0.1,96,0.1)),    # vector of x-values to calculate predictions and confidence intervals
   nsmpl = 100,                 # number of samples for simulation/estimation methods

   log.x = FALSE,            # graphical setting only!
   log.y = FALSE            # graphical setting only!

)
# model fit:
grid.arrange(res.advan$gr)


# diagnostics:
res.advan$diagnostics[[1]]
res.advan$diagnostics[[2]]
res.advan$diagnostics[[3]]
res.advan$diagnostics[[4]]

}
