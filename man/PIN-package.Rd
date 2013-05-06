\name{PIN-package}
\alias{PIN-package}
\alias{PIN}
\docType{package}
\title{
Computes the probability of informed trading
}
\description{
This package uses information on trade initiation to compute the probability of informed trading of Easley and O'Hara (1987).
}
\details{
\tabular{ll}{
Package: \tab PIN\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-06\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
The package consists of the function pin_likelihood.R, which computes the log-likelihood for the PIN model. This takes as inputs the parameter values - in the order of epsilon, mu, alpha, delta - and the time series for the number of buy- and sell-initiated trades.
}
\author{
P. Zagaglia

Maintainer: P. Zagaglia <paolo.zagaglia@gmail.com>
}
\references{
D. Easley and M. O'Hara. Price, Trade Size, and Information in Securities Markets. Journal of Financial Economics, 19(1), 69-90, 1987.
}


\keyword{ package }


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

#Pick and rename the optimizing parameters
alph <- c( param_optim$par[3] )
miu  <- c( param_optim$par[2] )
delt <- c( param_optim$par[4] )
epsi <- c( param_optim$par[1] )

#Compute PIN 
pin <- (alph*miu)/(alph*miu+2*epsi) 
#----------------------------------------------------------------------

}
