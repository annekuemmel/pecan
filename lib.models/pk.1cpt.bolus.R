# pk.1cpt.bolus.R


pk.1cpt.bolus <- function(
	x=seq(0,1,0.1),
	ke = 2.5,
	V = 2
){
	y <- 1 / V * exp(-ke*x)
	return(y)
}

pk.1cpt.bolus.deriv <- function(
	x=seq(0,1,0.1),
	ke = 2.5,
	V = 2
){


	# df/dke
	dke <- 1/V * - x * exp(-ke*x)
	
	# df/dV
	dV <- - 1/(V^2) * exp(-ke*x)
	

	y <- cbind(dke, dV)		
	return(y)

}