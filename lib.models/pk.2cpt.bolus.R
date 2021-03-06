pk.2cpt.bolus <- function(
	x=seq(0,1,0.1),
	ke = 2.5,
	V = 2,
	k12 = 2,
	k21 = 2
){
     s1 <- k12+k21+ke
     s2 <- k21*ke
     beta <- 0.5*(s1-sqrt(s1^2 - 4*s2))
     alpha <- s2 / beta
     A <- 1/V * (alpha-k21)/(alpha-beta)
     B <- 1/V * (beta-k21)/(beta-alpha)

	y <- A*exp(-alpha*x) + B*exp(-beta*x)
	return(y)
}

