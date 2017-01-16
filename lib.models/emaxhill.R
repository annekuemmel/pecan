emaxhill <- function(
	x=seq(0,1,0.1),
	E0 = 0.2,
	Emax = 1,
	EC50 = 0.5,
	g = 1.5
){
	y <- (Emax - E0) * x^g/(EC50^g + x^g) + E0
	return(y)
}

emaxhill.deriv <- function(
	x=seq(0,1,0.1),
	E0 = 0.2,
	Emax = 1,
	EC50 = 0.5,
	g = 1.5
){

	dE0 <- 1 - x^g / (EC50^g + x^g)
	dEmax <- x^g / (EC50^g + x^g)
	dEC50 <- x^g * (E0-Emax) * (EC50^g + x^g)^(-2) * g * EC50^(g-1)
	dg <- (E0-Emax) * (1 + (EC50/x)^g)^(-2) * (EC50/x)^g * log(EC50/x)

	y <- cbind(dE0, dEmax, dEC50, dg)		
	return(y)

}