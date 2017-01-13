calc.VarFun <- function(parvar, jac)
# calculates the variance of the function at given values
# of the indendent variable based on the delta method
# - parvar: variance of the parameters [m,m]
# - jac: jacobian wrt to the parameters at the indendent variable values [n,m]
# OBS! Jacobian has to have the first dimension n (number of supports)
# and the second dimension m (number of parameters)
{	
	# get dimensions
	n <- dim(jac)[1]
	m <- dim(jac)[2]

	# check whether parameter covariance matrix is symmetric
	if (!isSymmetric(parvar)) stop('Parameter covariance matrix is not symmetric.')

	# check whether parameter covariance matrix has the same dimension as the parameter dimension of the jacobian
      if (dim(parvar)[1] != m) stop('Number of parameters in Jacobian and Variance matrix not the same.')

      # check whether parameter names are same
      pars.j <- dimnames(jac)[[2]]
      pars.v <- dimnames(parvar)[[1]]
      if (!setequal(pars.j,pars.v)) stop('Parameter names for Jacobian and Variance matrix not the same.')

      ind_p <- match(pars.j,pars.v)
	Var <- (jac %*% parvar[ind_p,ind_p]) %*% t(jac)

	return(Var)
}