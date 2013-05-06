pin_likelihood <- function(params, n_trades){ 
	
	#Assign the number of buy- and sell-trades 
	B <- n_trades[,1]; 
	S <- n_trades[,2];
	
	#Assign number of trading days
	trad_days <- length( B );
	
	#Initialize parameter values
	epsi <- params[1];
	miu  <- params[2];
	alph <- params[3];
	delt <- params[4];
	
	#Initialize
	likel <- c(0);
	
	
	for (j in 1:trad_days){
		
		#Choose number of buy- and sell-trades for each trading day
		buy_s     <- B[j];
		sell_s    <- S[j];
		
		#Compute values of interest for the log-likelihood function
		M <- min( buy_s,sell_s ) + max( buy_s,sell_s )/2;
		x <- epsi/( miu+epsi );
		
		a1 <- exp(-miu);
		a2 <- x^(sell_s-M);
		a3 <- x^(buy_s-M);
		a4 <- x^(buy_s+sell_s-M);
		
		
		#Split the log-likelihood in two parts and compute the relevant terms
		part1 <- -2*epsi + M*log(x) + ( buy_s+sell_s )*log( miu+epsi );
		part2 <- log( alph*(1-alph)*( a1 )*( a2 ) + alph*delt*( a1 )*( a3 ) + (1-alph)*( a4 ) );
		
		#Compute the sum of log-likehoods over the window of trading days
		likel <- likel + ( part1+part2 )
	    
		#Delete past term for safety
		rm( part1, part2 ) 
	}
	
	#Impose boundary conditions and obtain the negative of the log-likelihood (for use with numerical minimization routines)
	if ((epsi>=0) && (miu>=0) && (alph>=0) && (delt>=0) && (epsi<=1) && (miu<=1) && (alph<=1) && (delt<=1)){
		likel_final <- -likel;
	}else{
		likel_final <- -Inf;
	}
	
	#Return the value of the approximated log-likelihood function 
	return( llikeval=likel_final )
	
}