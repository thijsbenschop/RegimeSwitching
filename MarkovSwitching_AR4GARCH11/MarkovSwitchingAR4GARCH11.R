## Clear history
rm(list = ls(all = TRUE))
graphics.off()

## Install and load packages
libraries = c("DEoptim")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
    install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# Set working directory
setwd("")

# Function to evaluate the loglikelihood, 2 states
condLogLikPar = function(theta, series){
	# calculates the conditional loglikelihood for a Markov Switching model 
	# with 2 states
	# Arguments
	# p11 transition probability to stay in state 1
	# p22 transition probability to stay in state 2
	# mu1, mu2, mean of distribution in resp. state 1 and 2
	# sigma1, sigma2, variance of distribution in resp. state 1 and 2
	# returns : value of conditional log-likelihood
	p11     = theta[1]
	p22     = theta[2]
	c1      = theta[3]
	c2      = theta[4]
	alpha01 = theta[5]
	alpha02 = theta[6]
	alpha11 = theta[7]
	alpha12 = theta[8]
	beta11  = theta[9]
	beta12  = theta[10]
	phi11   = theta[11]
	phi12   = theta[12]
	phi21   = theta[13]
	phi22   = theta[14]
	phi31   = theta[15]
	phi32   = theta[16]
	phi41   = theta[17]
	phi42   = theta[18]
	
	if ((alpha11 + beta11) >= 1 || (alpha12 + beta12) >= 1){logf = -999999999}
	# add restrictions on alpha 1 and beta 1
	else{
	n     = length(series)	# number of obs.
	p     = cbind(c(p11, 1-p22), c(1-p11, p22)) # transition matrix
	dzeta = cbind(rep(0, n), rep(0, n))
	f     = cbind(rep(0, n))
	eta   = cbind(rep(0,n), rep(0,n))
	
	dzetaInit = c((1-p[2,2])/(2-p[1,1]-p[2,2]), (1-p[1,1])/(2-p[2,2]-p[1,1]))
	# startvalue for iterations, assuming the Markov chain is ergodic
	# alternative  dzetaInit = c(0.5, 0.5)	
	
	# creating sigma and epsilon vectors
	sigma1    = cbind(rep(0,n))
	sigma2    = cbind(rep(0,n))
	sigma1[5] = sd(series)
	sigma2[5] = sd(series)
	epsilon21 = (series-c1)^2
	epsilon22 = (series-c2)^2
		#for (i in 2:n){
		#sigma1[i] = sqrt(alpha01 + (alpha11 * epsilon21[i-1]) + (beta11 * (sigma1[i-1])^2))
		#sigma2[i] = sqrt(alpha02 + (alpha12 * epsilon22[i-1]) + (beta12 * (sigma2[i-1])^2))
		#}
	#sigma1 = sqrt(sigma21)
	#sigma2 = sqrt(sigma22)
	
	for (i in 5:n){
		
		# Evaluate the densities under the two regimes
		############### create if else for different functional forms ##############
		mean1 = c1+phi11*series[i-1]+phi21*series[i-2] +phi31*series[i-3] +phi41*series[i-4] 
		mean2 = c2+phi12*series[i-1]+phi22*series[i-2] +phi32*series[i-3] +phi42*series[i-4]
		
		eta[i, 1] = dnorm(x=series[i], mean=mean1, sd=sigma1[i])
		eta[i, 2] = dnorm(x=series[i], mean=mean2, sd=sigma2[i])
		
		# Evaluate the conditional density of the ith observation
		if (i == 5){
		f[i] = t(p %*% c(eta[i,1], eta[i,2])) %*% dzetaInit
		}
		else{
		f[i] = t(p %*% c(eta[i,1], eta[i,2])) %*% c(dzeta[i-1, 1], dzeta[i-1, 2])
		}
		# Evaluate the state probabilities
		if(i==5){
		dzeta[i, 1] = dzetaInit[1]
		dzeta[i, 2]	= dzetaInit[2]
		}
		else{
		dzeta[i, 1] = (eta[i,1] * (p[,1] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		dzeta[i, 2] = (eta[i,2] * (p[,2] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		}

		# Calculating sigma2
		if(i == 5){
		sigma1[i+1] = sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		sigma2[i+1] = sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		}	
		else{
		sigma1[i+1] = sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		sigma2[i+1] = sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		}	
	}
	logf = sum(log(f[5:n]))
	}
	if(is.nan(logf)==TRUE){
		cat("Error : Returned not a number ", "\n")
		flush.console()
		logf = -999999999
		}
	#output = cbind(eta, f, dzeta)
	cat(logf, "\n")
	#cat(p11, " ", p22, " ", c1, " ", c2, "\n ", alpha01, " ", alpha02, " ", alpha11, " ", alpha12, "\n ", beta11, " ", beta12, "\n ", logf, "\n")
	flush.console()
	return(list(dzeta, sigma1, sigma2, epsilon21, epsilon22))	
}


# Function to evaluate the loglikelihood, 2 states, both AR(4)-GARCH(1,1)
condLogLik = function(theta, series){
	# calculates the conditional loglikelihood for a Markov Switching model 
	# with 2 states
	# Arguments
	# p11 transition probability to stay in state 1
	# p22 transition probability to stay in state 2
	# mu1, mu2, mean of distribution in resp. state 1 and 2
	# sigma1, sigma2, variance of distribution in resp. state 1 and 2
	# returns : value of conditional log-likelihood
	p11     = theta[1]
	p22     = theta[2]
	c1      = theta[3]
	c2      = theta[4]
	alpha01 = theta[5]
	alpha02 = theta[6]
	alpha11 = theta[7]
	alpha12 = theta[8]
	beta11  = theta[9]
	beta12  = theta[10]
	phi11   = theta[11]
	phi12   = theta[12]
	phi21   = theta[13]
	phi22   = theta[14]
	phi31   = theta[15]
	phi32   = theta[16]
	phi41   = theta[17]
	phi42   = theta[18]
	
	if ((alpha11 + beta11) >= 1 || (alpha12 + beta12) >= 1){logf = -999999999}
	# add restrictions on alpha 1 and beta 1
	else{
	n     = length(series)	# number of obs.
	p     = cbind(c(p11, 1-p22), c(1-p11, p22)) # transition matrix
	dzeta = cbind(rep(0, n), rep(0, n))
	f     = cbind(rep(0, n))
	eta   = cbind(rep(0,n), rep(0,n))
	
	dzetaInit = c((1-p[2,2])/(2-p[1,1]-p[2,2]), (1-p[1,1])/(2-p[2,2]-p[1,1]))
	# startvalue for iterations, assuming the Markov chain is ergodic
	# alternative  dzetaInit = c(0.5, 0.5)	
	
	# creating sigma and epsilon vectors
	sigma1    = cbind(rep(0,n))
	sigma2    = cbind(rep(0,n))
	sigma1[5] = sd(series)
	sigma2[5] = sd(series)
	epsilon21 = (series-c1)^2
	epsilon22 = (series-c2)^2
		#for (i in 2:n){
		#sigma1[i] = sqrt(alpha01 + (alpha11 * epsilon21[i-1]) + (beta11 * (sigma1[i-1])^2))
		#sigma2[i] = sqrt(alpha02 + (alpha12 * epsilon22[i-1]) + (beta12 * (sigma2[i-1])^2))
		#}
	#sigma1 = sqrt(sigma21)
	#sigma2 = sqrt(sigma22)
	
	for (i in 5:n){
		
		# Evaluate the densities under the two regimes
		############### create if else for different functional forms ##############
		mean1 = c1+phi11*series[i-1]+phi21*series[i-2] +phi31*series[i-3] +phi41*series[i-4] 
		mean2 = c2+phi12*series[i-1]+phi22*series[i-2] +phi32*series[i-3] +phi42*series[i-4]
		
		eta[i, 1] = dnorm(x=series[i], mean=mean1, sd=sigma1[i])
		eta[i, 2] = dnorm(x=series[i], mean=mean2, sd=sigma2[i])
		
		# Evaluate the conditional density of the ith observation
		if (i == 5){
		f[i] = t(p %*% c(eta[i,1], eta[i,2])) %*% dzetaInit
		}
		else{
		f[i] = t(p %*% c(eta[i,1], eta[i,2])) %*% c(dzeta[i-1, 1], dzeta[i-1, 2])
		}
		# Evaluate the state probabilities
		if(i==5){
		dzeta[i, 1] = dzetaInit[1]
		dzeta[i, 2]	= dzetaInit[2]
		}
		else{
		dzeta[i, 1] = (eta[i,1] * (p[,1] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		dzeta[i, 2] = (eta[i,2] * (p[,2] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		}

		# Calculating sigma2
		if(i == 5){
		sigma1[i+1] = sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		sigma2[i+1] = sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		}	
		else{
		sigma1[i+1] = sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		sigma2[i+1] = sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		}	
	}
	logf = sum(log(f[5:n]))
	}
	if(is.nan(logf)==TRUE){
		cat("Error : Returned not a number ", "\n")
		flush.console()
		logf = -999999999
		}
	#output = cbind(eta, f, dzeta)
	#cat(logf, "\n")
	#cat(p11, " ", p22, " ", c1, " ", c2, "\n ", alpha01, " ", alpha02, " ", alpha11, " ", alpha12, "\n ", beta11, " ", beta12, "\n ", logf, "\n")
	#flush.console()
	return(-logf)	
}


#Read-in data
data = read.table("dataCO2.txt", header=TRUE)

#Coefficient(s) of GARCH(1,1) without regime switching:
# constant, alpha0, alpha1, beta1
#         mu        omega       alpha1        beta1  
#-2.6460e-04   4.6171e-06   7.2568e-02   9.1988e-01  

# Test function with values of GARCH(1,1)
test = condLogLik(theta = c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988, 0.0988, 0.0988, -0.1391, -0.1391, 0.0795, 0.0795, 0.0609 , 0.0609), series = data[2:725, 5])
test
test = condLogLik(theta = parDE, series = data[2:725, 5])
test

#Plot state probs
dzetaTest = condLogLikDzeta(theta = c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988, 0.0988, 0.0988, -0.1391, -0.1391, 0.0795, 0.0795, 0.0609 , 0.0609), series = data[2:725, 5])
plot(dzetaTest[[1]][,1], type="l")

## Defining constaints
# Matrix with lienar combinations of parameters (p11. p22, c1, c2, alpha01, alpha02, alpa11, alpha12, beta11, beta12, phi11, phi12, phi21, phi22, phi31, phi32, phi41, phi42)
# (k x p) matrix, k constraints, p parameters 
#					
constraintMat = matrix(c(0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 1, 0, 0, 0, 0, 0, 0, 0,
						 0, 1, 0, 0, 0, 0, 0, 0,
						 0, 0, 1, 0, 0, 0, -1, 0,
						 0, 0, 0, 1, 0, 0, 0, -1,
						 0, 0, 0, 0, 1, 0, -1, 0,
						 0, 0, 0, 0, 0, 1, 0, -1,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0,
						 0, 0, 0, 0, 0, 0, 0, 0 ), 8, 18)

constraintVec = c(0, 0, 0, 0, 0, 0, -1, -1)

logRetFirstPeriod = data[2:725, 5]

## Other optimization alghorithms
# deoptim
library(DEoptim)
#		p11, p22, c1,   c2,   alf01, alf02, alf11, alf12, beta12, beta22, phi11,  phi12,  phi21,  phi22,  phi31,  phi32,  phi41,  phi42
lowParDE   = c(0,   0,   -0.3, -0.3, 0,     0,     0,     0,     0,    0, 	-1,	  -1,     -1,     -1, 	  -1,     -1,     -1,     -1)
upParDE    =  c(1,   1,    0.3,  0.3, 0.2,   0.2,   1,     1,     1,    1, 	 1,	   1,      1,      1, 	   1,      1,      1,     1)
controlDE  = list(NP = 240, itermax = 2500)
testDE     = DEoptim(fn      = condLogLik,
                     lower   = lowParDE,
                     upper   = upParDE,
                     series  = logRetFirstPeriod,
                     control = controlDE)
# Result
#Iteration: 1908 bestvalit: -1750.935434 bestmemit: 0.974023 0.881826 0.001113 -0.009027 0.000013 0.000223 0.007800 0.195207 0.864445 0.750989 -0.033932 0.301394 -0.063675 -0.210758 0.026095 0.196542 -0.031466 0.251212
#Iteration: 1638 bestvalit: -1750.834437 bestmemit: 0.974023 0.881826 0.001101 -0.009656 0.000013 0.000223 0.007800 0.204368 0.864445 0.750989 -0.022887 0.301394 -0.063675 -0.210758 0.028724 0.214878 -0.041642 0.251212
parDE = c(0.974023, 0.881826, 0.001101, -0.009656, 0.000013 , 0.000223, 0.007800, 0.204368, 0.864445, 0.750989, -0.022887, 0.301394, -0.063675, -0.210758, 0.028724, 0.214878, -0.041642, 0.251212)

# unconditional mean and sd
# state 1
parDE[3] / (1-parDE[11]-parDE[13]-parDE[15]-parDE[17])
sqrt(parDE[5] / (1-parDE[7]-parDE[9]))
# state 2
parDE[4] / (1-parDE[12]-parDE[14]-parDE[16]-parDE[18])
sqrt(parDE[6] / (1-parDE[8]-parDE[10]))

AIC = 2*(-1750.935434) + 2*18

# Constrained optimization
test7  = constrOptim(c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988, 0.0988, 0.0988, -0.1391, -0.1391, 0.0795, 0.0795, 0.0609 , 0.0609), ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test8  = constrOptim(test7$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test9  = constrOptim(test8$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test10 = constrOptim(test9$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test11 = constrOptim(test10$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test12 = constrOptim(test11$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test13 = constrOptim(test12$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test14 = constrOptim(test13$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test15 = constrOptim(test14$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test16 = constrOptim(test15$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test17 = constrOptim(test16$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test18 = constrOptim(test17$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)

cbind(test7$par, test8$par, test9$par, test10$par, test11$par, test12$par, test13$par, test14$par, test15$par, test16$par, test17$par, test18$par) 

####
#results of test 18	
parValues = test18$par

#$value
#[1] -2492.734

#$par
# [1]  1.604017e+00  3.607635e-01 -2.251704e-05 -1.613421e-01  4.475618e-10
# [6]  1.838259e+00  7.254282e-02  7.256048e-02  9.198838e-01  9.198772e-01
#[11]  1.058183e-02  4.204175e-02  7.197604e-03  1.524742e-02  6.525468e-02
#[16]  1.096865e-01 -1.724443e-03  1.949872e-01

final      = condLogLik(theta  = parValues,
                        series = data[2:725, 5])
dzetaFinal = condLogLikDzeta(theta  = parValues,
                             series = data[2:725, 5])

plot(dzetaFinal[,1],
     type = "l")

#different starting values
#test9 = constrOptim(c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.01, 0.01, 0.98, 0.98), ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)

#########################################################
# Plot dzeta, prob in high regime and plot of series
# in-sample period
dzeta  = condLogLikPar(theta  = parDE,
                       series = data[2:725, 5]) #calculate dzeta
dzeta1 = dzeta[[1]][,1]
dzeta2 = dzeta[[1]][,2]

par(mfrow=c(2,1))
jan08 = as.Date("01/01/08", "%d/%m/%y")
dec10 = as.Date("01/01/11", "%d/%m/%y")
DatePlot = as.Date(data$Date[6:725], , "%d/%m/%y")

par(mar = c(4,4,2,2))
plot(dzeta1[5:724]~DatePlot,
     type = "l",
     ylim = c(0.0,1.0),
     lwd  = 1.2,
     xlab = "Date",
     ylab = "Regime probabilities",
     xaxs = "i",
     yaxs = "i",
     xlim = c(jan08, dec10),
     yaxt = "n")

axis(side = 2,
     at   = c(0, 0.5, 1),
     las  = 1)

par(mar = c(4,4,0.5,2))
plot(data[6:725,5]~DatePlot,
     type = "l",
     lwd  = 1.2,
     xaxs = "i",
     yaxs = "i",
     xlim = c(jan08, dec10),
     xlab = "Date", 
     ylab = "Log returns",
     ylim = c(-0.15, 0.15),
     yaxt = "n")

axis(side = 2,
     at   = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15),
     las  = 1)

prob = cbind(dzeta1, dzeta2)
#write.table(prob, "probMSAR4GARCH.txt")

#compute the fitted values
series = data[2:725, 5]
fitted.y = rep(0,724)
fitted.y[1:4]=NA
for(i in 5:724){
		mean1 = parDE[3]+parDE[11]*series[i-1]+parDE[13]*series[i-2] +parDE[15]*series[i-3] +parDE[17]*series[i-4] 
		mean2 = parDE[4]+parDE[12]*series[i-1]+parDE[14]*series[i-2] +parDE[16]*series[i-3] +parDE[18]*series[i-4]		
		fitted.y[i] = dzeta1[i] * mean1 + dzeta2[i] * mean2
} 

par(new=TRUE)
plot(fitted.y,
     type = "l",
     ylim = c(-0.12, 0.12),
     col  = "red",
     xaxt = "n",
     yaxt = "n",
     ann  = FALSE)


#########################################################
####################### Forecasting #####################
#########################################################
### point forecast for period 2011-2012 (726-1183)
#2# reestimation, recursive
#matrix with estimate, error, estimated variance
recursivePar = parValues
prevPar = parValues
# reestimate the parameters for the longer time series
for (i in 726:1183){
	reestMSGARCH11 = constrOptim(prevPar,
	                             ui     = constraintMat,
	                             ci     = constraintVec,
	                             f      = condLogLik,
	                             series = data[2:i, 5],
	                             grad   = NULL)
	
	recursivePar = cbind(recursivePar, reestMSGARCH11$par) 
	prevPar      = reestMSGARCH11$par
	cat("iteration: ", i, "/n")
}

# load recursivePar
recursivePar = read.table("reestRecurMSAR4GARCH11.txt", header=TRUE)

# calculate point and density estimates, confidence intervals
pointForecastRec = rep(0, 458)
sigma2Forecast   = rep(0, 458)
p1Forecast       = rep(0, 458)
p2Forecast       = rep(0, 458)
c1Forecast       = rep(0, 458)
c2Forecast       = rep(0, 458)
dzeta1           = rep(0, 458)
dzeta2           = rep(0, 458)

for (i in 726:1183){	
	# Forecasting values for i
	# calculate dzeta[i], sigma[i]
	info = condLogLikPar(theta  = recursivePar[,(i-725)],
	                     series = data[2:(i-1), 5]) #calculate dzeta and sigma
	curPar = recursivePar[,(i-725)]	# estimated parameters
	dzeta1Cur = info[[1]][(i-2),1]
	dzeta2Cur = info[[1]][(i-2),2]
	p1 = dzeta1Cur * curPar[1] + dzeta2Cur * (1-curPar[2]) # prob of being in state 1 in forecast
	p2 = dzeta2Cur * curPar[2] + dzeta1Cur * (1-curPar[1]) # prob of being in state 2 in forecast
	p1Forecast[(i-725)] = p1
	meanCur1 = curPar[3] + curPar[11] * data[(i-1),5] + curPar[13] * data[(i-2),5]+ curPar[15] * data[(i-3),5]+ curPar[17] * data[(i-4),5]
	meanCur2 = curPar[4] + curPar[12] * data[(i-1),5] + curPar[14] * data[(i-2),5]+ curPar[16] * data[(i-3),5]+ curPar[18] * data[(i-4),5]
	c1Forecast[(i-725)] =meanCur1
	c2Forecast[(i-725)] =meanCur2
	#point estimate
	pointForecastRec[i-725] = p1 * curPar[3] + p2 * curPar[4]
	sigma2Cur = dzeta1Cur * (info[[2]][(i-2)])^2 + dzeta2Cur * (info[[3]][(i-2)])^2
	sigma2Forecast[i-725] = p1 * (curPar[5] + (curPar[7] * info[[4]][(i-2)]) + curPar[9]*sigma2Cur) + p2 * (curPar[6] + (curPar[8] * info[[5]][(i-2)]) + curPar[10]*sigma2Cur)
	}

estErrorRecur = pointForecastRec - data$logRet[726:1183]
resPlot       = cbind(estErrorRecur, sigma2Forecast)
#write.table(resPlot, "resMSAR4GARCH11.txt")

recursiveMSE = (1/length(estErrorRecur)) * sum(estErrorRecur^2)
recursiveMAE = (1/length(estErrorRecur)) * sum(abs(estErrorRecur))
recursiveMSE
recursiveMAE

# logreturns and predicted 95 percent confidence intervals
# calculate 2,5 and 97,5 quantiles
confInterval = cbind(rep(0, 458), rep(0, 458))
for (i in 1:458){
	confInterval[i,1] = qnorm(0.025, mean=pointForecastRec[i], sd=sqrt(sigma2Forecast[i]))
	confInterval[i,2] = qnorm(0.975, mean=pointForecastRec[i], sd=sqrt(sigma2Forecast[i]))
	}

conf = cbind(confInterval, pointForecastRec)
#write.table(conf, "confMSAR4GARCH.txt")

plot(data[726:1183,5],
     type = "l",
     ylim = c(-0.5, 0.5))

par(new = TRUE)
plot(confInterval[2:458,1],
     type = "l",
     ylim = c(-0.5, 0.5))

par(new = TRUE)
plot(confInterval[2:458,2],
     type = "l",
     ylim = c(-0.5, 0.5))

par(new = TRUE)
plot(pointForecastRec,
     type = "l",
     ylim = c(-0.3, 0.3),
     col  = "red",
     xaxt = "n",
     yaxt = "n",
     ann  = FALSE)

## Checking univformity of density forecasts
# Density transformation
u = c(rep(0, 458))
for (i in 1:length(u)){
	u[i] = pnorm(data$logRet[i+725],
	             mean = pointForecastRec[i],
	             sd   = sqrt(sigma2Forecast[i]))
}

hist(u, breaks = 20)

#write.table(u, "uMSAR4GARCH11.txt")

# test on uniformity
#Kolmogorov Smirnov
ks.test(u, "punif")