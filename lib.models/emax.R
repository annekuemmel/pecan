emax <- function(
	x=seq(0,1,0.1),
	E0 = 0.2,
	Emax = 1,
	EC50 = 0.5
){
	y <- (Emax - E0) * x/(EC50 + x) + E0
	return(y)
}

emax.deriv <- function(
	x=seq(0,1,0.1),
	E0 = 0.2,
	Emax = 1,
	EC50 = 0.5
){

	dE0 <- 1 - x / (EC50 + x)
	dEmax <- x / (EC50 + x)
	dEC50 <- x * (E0-Emax) * (EC50 + x)^(-2)

	y <- cbind(dE0, dEmax, dEC50)		
	return(y)

}