# x: parameter on original scale
# p: parameter on transformed scale (used for estimation)
# se: standard error for p (on transformed scale)


     trans.par <- function(x,trans)
        switch(trans,
               log = log(x),
               logit = logit(x),
               x)

     transinv.par <- function(p,trans)
        switch(trans,
               log = exp(p),
               logit = logitinv(p),
               p)

     transinv.se <- function(se,x,trans)
        switch(trans,
               log = se*x,
               logit = x*(1-x)*se,
               se)

     trans.mod <- function(trans)
        switch(trans,
               log = "exp",
               logit = "logitinv",
               "c")

# auxiliary transformation functions
logit <- function(x) log(x/1-x) 
logitinv <- function(p) 1/(1+exp(-p))
