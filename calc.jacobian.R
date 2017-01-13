calc.jacobian <- function(fun, pars, xsupport)
# numerical calculation of partial derivative of a function wrt its parameters
# - fun: function (needs to accept the parameters as names arguments
# - pars: named list of parameters
# - xsupport: vector of independent variable values at which to calculate the derivative
{
      require(numDeriv)
   
	m <- length(pars)
	parnames <- names(pars)
	parstrdf <- data.frame(names=parnames, pars=paste('pars$',parnames,sep=''))
	exprstr <- paste('fun_i <- function(x) fun(x=xsupport, ', paste(parstrdf$names, parstrdf$pars, collapse = ', ', sep = '=') ,' )')
	
	for (i in 1:m){
		exprstr_i <- sub(parstrdf$pars[i],'x',exprstr, fixed = TRUE)
		eval(parse(text=exprstr_i))
		numderiv_i <- jacobian(fun_i, pars[[parnames[i]]])

		if (i==1) numderiv <- numderiv_i else numderiv <- cbind(numderiv, numderiv_i)
	}
	dimnames(numderiv) <- list(NULL, parnames)
	return(numderiv)
}