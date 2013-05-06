\name{pin_likelihood}
\alias{pin_likelihood}

\title{
Compute value of log-likelihood function for the estimation of the structural model of Easley and O'Hara (1987).
}
\description{
This function computes the value of the log-likelihood function (negative of) for the model of strategic trading with asymmetric information proposed by Easley and O'Hara (1987). This is an application of the sequential trade model originally introduced by Kyle (1985). To deal with possible convergence issues in numerical optimization, the code implements the modified likelihood proposed by Aktasa et al. (2007).
}
\usage{
pin_likelihood(params, n_trades)
}
\arguments{
  \item{params}{
Parameters of the structural model in the following order: epsilon, mu, alpha, delta.
}
  \item{n_trades}{
A data matrix with a number of rows equal to the trading days, and two columns. Buy-initiated trades are stored in the first column. Sell-initiated trades are indexed by the second column. 
}
}
\details{
A data matrix with a number of rows equal to the trading days, and two columns. Buy-initiated trades are stored in the first column. Sell-initiated trades are indexed by the second column.  
}
\value{
\item{likel_final }{Value of the modified log-likelihood function.}
}
\references{
N. Aktasa, E. de Bodta, F. Declerckc and H. Van Oppensa. The PIN Anomaly Around M&A Announcements, Journal of Financial Markets, 10(2), 169-191, 2007.

D. Easley and M. O'Hara. Price, Trade Size, and Information in Securities Markets. Journal of Financial Economics, 19(1), 69-90, 1987.

D. Easley, R. F. Engle, M. O'Hara and L. Wu. Time-Varying Arrival Rates of Informed and Uninformed Traders. Journal of Financial Econometrics, 6(2), 171-207, 2008.
D. Easley, N. M. Kiefer, M. O'Hara and J. B. Paperman. Liquidity, Information, and Infrequently Traded Stocks. Journal of Finance, 51(4), 1405-1436, 1996.

L. R. Glosten and L. E. Harris. Estimating the Components of the Bid/Ask Spread. Journal of Financial Economics, 21(1), 123-142, 1988.

A. S. Kyle. Continuous Auctions and Insider Trading. Econometrica, 53(6), 1315-1335, 1985.
}
\author{
Paolo Zagaglia <paolo.zagaglia@gmail.com>
}



\examples{

#----------------------------------------------------------------------
# Preliminaries
#Load relevant libraries
library(zoo)
library(xts)
library(highfrequency)
library(TTR)
library(timeDate)
library(quantmod)

#Load data samples
data("sample_tdata");
data("sample_qdata");
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#Exploit the functionalities of the highfrequency package

#Match the trade and quote data
tqdata <- matchTradesQuotes(sample_tdata,sample_qdata);

#Get the inferred trade direction according to Lee-Ready classification
x <- getTradeDirection(tqdata);
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Compute the number of buy- and sell-initiated trades

tradeDirection <- matrix(x);

#Select buy-initiated trades
buy_side <- which( tradeDirection>0 );

#Compute the number of buy- and sell-trades 
num_buy_side <- length( matrix(buy_side) )
num_sell_side <- length( tradeDirection ) - length( matrix(buy_side) )

#Summarize the number of trades directions 
ntrades <- cbind(num_buy_side,num_sell_side)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Optimize the log-likelihood and compute PIN 

#Set the initial values for numerical optimization
initparams <- c(0.15, 0.05, 0.5, 0.5) 

#Compute the miniming parameter values for the negative of 
#the log-likelihood function
options(warn=-1)
param_optim <- optim(initparams, pin_likelihood, gr=NULL, ntrades)

}


\keyword{ pin_likelihood }