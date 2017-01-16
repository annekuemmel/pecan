# pk.1cpt.abs.R


pk.1cpt.oral <- function(
	x=seq(0,1,0.1),
	ka = 5,
	ke = 2.5,
	V = 2
){
	y <- ka / (V * (ka-ke)) * (exp(-ke*x) - exp(-ka*x))
	return(y)
}

pk.1cpt.oral.deriv <- function(
	x=seq(0,1,0.1),
	ka = 5,
	ke = 2.5,
	V = 2
){

	a <- ka / (V * (ka-ke))
	b <- exp(-ke*x) - exp(-ka*x)

	# df/dka
	DaDka <- 1/V * ( 1/(ka-ke) - ka / (ka-ke)^2 )
	DbDka <- x * exp(-ka*x)
	dka <- DaDka*b + a*DbDka

	# df/dke
	DaDke <- ka/V * 1/(ka-ke)^2
	DbDke <- - x * exp(-ke*x)
	dke <- DaDke*b + a*DbDke
	
	# df/dV
	DaDV <- - ka / (ka-ke) * 1/V^2
	dV <- DaDV*b
	

	y <- cbind(dka, dke, dV)		
	return(y)

}