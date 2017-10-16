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
	
	if ((alpha11 + beta11) >= 1 || (alpha12 + beta12) >= 1){logf = -9999999}
		# add restrictions on alpha 1 and beta 1
	else{
	n       = length(series)	# number of obs.
	p       = cbind(c(p11, 1-p22), c(1-p11, p22)) # transition matrix
	dzeta   = cbind(rep(0, n), rep(0, n))
	f       = cbind(rep(0, n))
	eta     = cbind(rep(0,n), rep(0,n))
	
	dzetaInit = c((1-p[2,2])/(2-p[1,1]-p[2,2]), (1-p[1,1])/(2-p[2,2]-p[1,1]))
	# startvalue for iterations, assuming the Markov chain is ergodic
	# alternative  dzetaInit = c(0.5, 0.5)	
	
	# creating sigma and epsilon vectors
	sigma1 = cbind(rep(0,n))
	sigma2 = cbind(rep(0,n))
	sigma1[1] = sd(series)
	sigma2[1] = sd(series)
	epsilon21 = (series-c1)^2
	epsilon22 = (series-c2)^2
		#for (i in 2:n){
		#sigma1[i] = sqrt(alpha01 + (alpha11 * epsilon21[i-1]) + (beta11 * (sigma1[i-1])^2))
		#sigma2[i] = sqrt(alpha02 + (alpha12 * epsilon22[i-1]) + (beta12 * (sigma2[i-1])^2))
		#}
	#sigma1 = sqrt(sigma21)
	#sigma2 = sqrt(sigma22)
	
	for (i in 1:n){
		
		# Evaluate the densities under the two regimes
		############### create if else for different functional forms ##############
		eta[i, 1] = dnorm(x=series[i], mean=c1, sd=sigma1[i])
		eta[i, 2] = dnorm(x=series[i], mean=c2, sd=sigma2[i])
		
		# Evaluate the conditional density of the ith observation
		if (i == 1){
		f[i] = t(p %*% c(eta[i,1], eta[i,2])) %*% dzetaInit
		}
		else{
		f[i] = t(p %*% c(eta[i,1], eta[i,2])) %*% c(dzeta[i-1, 1], dzeta[i-1, 2])
		}
		# Evaluate the state probabilities
		if(i==1){
		dzeta[i, 1] = dzetaInit[1]
		dzeta[i, 2]	= dzetaInit[2]
		}
		else{
		dzeta[i, 1] = (eta[i,1] * (p[,1] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		dzeta[i, 2] = (eta[i,2] * (p[,2] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		}

		# Calculating sigma2
		if(i == 1){
		sigma1[i+1] = sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		sigma2[i+1] = sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		}	
		else{
		sigma1[i+1] = sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		sigma2[i+1] = sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		}	
	}
	logf = sum(log(f) * rep(1,n))
	}
	if(is.nan(logf)){
		logf = -9999999
		cat("Error not a number!!!", "\n")
		}
	#output = cbind(eta, f, dzeta)
	cat(logf, "\n")
	#cat(p11, " ", p22, " ", c1, " ", c2, "\n ", alpha01, " ", alpha02, " ", alpha11, " ", alpha12, "\n ", beta11, " ", beta12, "\n ", logf, "\n")
	flush.console()
	return(list(dzeta, sigma1, sigma2, epsilon21, epsilon22))	
}


# Function to evaluate the loglikelihood, 2 states, both GARCH(1,1)
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
	
	if ((alpha11 + beta11) >= 1 || (alpha12 + beta12) >= 1){logf = -9999999}
		# add restrictions on alpha 1 and beta 1
	else{
	n       = length(series)	# number of obs.
	p       = cbind(c(p11, 1-p22), c(1-p11, p22)) # transition matrix
	dzeta   = cbind(rep(0, n), rep(0, n))
	f       = cbind(rep(0, n))
	eta     = cbind(rep(0,n), rep(0,n))
	
	dzetaInit = c((1-p[2,2])/(2-p[1,1]-p[2,2]), (1-p[1,1])/(2-p[2,2]-p[1,1]))
	# startvalue for iterations, assuming the Markov chain is ergodic
	# alternative  dzetaInit = c(0.5, 0.5)	
	
	# creating sigma and epsilon vectors
	sigma1 = cbind(rep(0,n))
	sigma2 = cbind(rep(0,n))
	sigma1[1] = sd(series)
	sigma2[1] = sd(series)
	epsilon21 = (series-c1)^2
	epsilon22 = (series-c2)^2
		#for (i in 2:n){
		#sigma1[i] = sqrt(alpha01 + (alpha11 * epsilon21[i-1]) + (beta11 * (sigma1[i-1])^2))
		#sigma2[i] = sqrt(alpha02 + (alpha12 * epsilon22[i-1]) + (beta12 * (sigma2[i-1])^2))
		#}
	#sigma1 = sqrt(sigma21)
	#sigma2 = sqrt(sigma22)
	
	for (i in 1:n){
		
		# Evaluate the densities under the two regimes
		############### create if else for different functional forms ##############
		eta[i, 1] = dnorm(x=series[i], mean=c1, sd=sigma1[i])
		eta[i, 2] = dnorm(x=series[i], mean=c2, sd=sigma2[i])
		
		# Evaluate the conditional density of the ith observation
		if (i == 1){
		f[i] = t(p %*% c(eta[i,1], eta[i,2])) %*% dzetaInit
		}
		else{
		f[i] = t(p %*% c(eta[i,1], eta[i,2])) %*% c(dzeta[i-1, 1], dzeta[i-1, 2])
		}
		# Evaluate the state probabilities
		if(i==1){
		dzeta[i, 1] = dzetaInit[1]
		dzeta[i, 2]	= dzetaInit[2]
		}
		else{
		dzeta[i, 1] = (eta[i,1] * (p[,1] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		dzeta[i, 2] = (eta[i,2] * (p[,2] %*% c(dzeta[i-1, 1], dzeta[i-1, 2]))) /f[i]
		}

		# Calculating sigma2
		if(i == 1){
		sigma1[i+1] = sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		sigma2[i+1] = sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzetaInit[1]*(sigma1[i])^2 + dzetaInit[2]*(sigma2[i])^2)))
		}	
		else{
		sigma1[i+1] = sqrt(alpha01 + (alpha11 * epsilon21[i]) + (beta11 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		sigma2[i+1] = sqrt(alpha02 + (alpha12 * epsilon22[i]) + (beta12 * (dzeta[i,1]*(sigma1[i])^2 + dzeta[i,2]*(sigma2[i])^2)))
		}	
	}
	logf = sum(log(f) * rep(1,n))
	}
	if(is.nan(logf)){
		logf = -9999999
		cat("Error not a number!!!", "\n")
		}
	#output = cbind(eta, f, dzeta)
	#cat(logf, "\n")
	#cat(p11, " ", p22, " ", c1, " ", c2, "\n ", alpha01, " ", alpha02, " ", alpha11, " ", alpha12, "\n ", beta11, " ", beta12, "\n ", logf, "\n")
	#flush.console()
	return(-logf)	
}

# Read-in data
data = read.table("dataCO2.txt", header=TRUE)

#Coefficient(s) of GARCH(1,1) without regime switching:
# constant, alpha0, alpha1, beta1
#         mu        omega       alpha1        beta1  
#-2.6460e-04   4.6171e-06   7.2568e-02   9.1988e-01  

# Test function with values of GARCH(1,1)
test = condLogLik(theta = c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), series = data[2:725, 5])
test
#testa = condLogLik(theta = parDE, series = logRetFirstPeriod)


#Plot state probs
#dzetaTest = condLogLikDzeta(theta = parValues, series = data[2:725, 5])
#plot(dzetaTest[,1], type="l")

## Defining constaints
# Matrix with lienar combinations of parameters (p11. p22, c1, c2, alpha01, alpha02, alpa11, alpha12, beta11, beta12)
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
						  0, 0, 0, 0, 0, 1, 0, -1
						  ), 8, 10)
constraintMat
constraintVec =c(0, 0, 0, 0, 0, 0, -1, -1)

logRetFirstPeriod = data$logRet[2:725]

# p11, p22, c1, c2, alpha01, alpha02, alpha11, alpha12, beta11, beta12

## Other optimization alghorithms
# deoptim
library(DEoptim)
#			  p11, p22,  c1,   c2,  alf01, alf02, alf11, alf12, beta12, beta22
lowParDE = c(0,   0,   -0.3, -0.3, 0,     0,     0,     0,     0,      0)
upParDE =  c(1,   1,    0.3,  0.3, 0.2,   0.2,   1,     1,     1,      1)
controlDE = list(NP = 120, itermax = 500)
testDE = DEoptim(fn=condLogLik, lower=lowParDE, upper=upParDE, series=logRetFirstPeriod, control=controlDE)
# Result
#Iteration: 414 bestvalit: -1749.445659 bestmemit:    0.982090    0.992281   -0.004157    0.000865    0.000289    0.000052    0.103816    0.001251    0.723247    0.716569
parDE = c(0.982090,    0.992281,   -0.004157,    0.000865,    0.000289,    0.000052,    0.103816,    0.001251,    0.723247,    0.716569)
#Iteration: 500 bestvalit: -1749.494048 bestmemit:    0.987508    0.977368    0.000958   -0.004848    0.000055    0.000276    0.001098    0.090699    0.683173    0.751255
#Iteration: 699bestvalit -1749,399347
#parDE = c(0.979766, 0.987543, -0.004397, 0.000905, 0.000243, 0.000067, 0.091176, 0.000013, 0.765748, 0.611604)

llh = condLogLik(parDE, logRetFirstPeriod)

llhcor = llh/724*720
AIC = 2*llhcor + 2*10

# unconditonal standard deviation
# state 1
sqrt(parDE[5]/(1-parDE[7]-parDE[9])) #falsh herum?
# state 2
sqrt(parDE[6]/(1-parDE[8]-parDE[10]))

# unconditional probability
# state 1
(1-parDE[1])/(2-parDE[1]-parDE[2])
# state 2
(1-parDE[2])/(2-parDE[1]-parDE[2])


# nlminb
startPar = c(0.8, 0.2, -0.00026460, 0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988)
lowPar = c(0,0, -1, -1, 0, 0, 0, 0, 0, 0)
upPar = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
testNlminb = nlminb(start     = startPar,
                     objective = condLogLik,
                     series    = logRetFirstPeriod,
                     lower     = lowPar,
                     upper     = upPar)
testNlminb$par

# neldermead
library(neldermead)
fminsearch(fun=condLogLik, x0=startPar)
#?nlminb

#dfoptim
#?dfoptim

#Unconstrained optimization
test6 = optim(c(0.8, 0.2, -0.00026460, 0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988), condLogLik, series = logRetFirstPeriod)
test6b = optim(test6$par, condLogLik, series = logRetFirstPeriod)
test6c = optim(test6b$par, condLogLik, series = logRetFirstPeriod)
test6d = optim(test6c$par, condLogLik, series = logRetFirstPeriod)
test6e = optim(test6d$par, condLogLik, series = logRetFirstPeriod)
test6f = optim(test6e$par, condLogLik, series = logRetFirstPeriod)
test6g = optim(test6f$par, condLogLik, series = logRetFirstPeriod)
test6h = optim(test6g$par, condLogLik, series = logRetFirstPeriod)

initPar = c(0.8, 0.2, -0.00026460, 0.00026460, 0.0000046171, 0.0000046171, 0.072568, 0.072568, 0.91988, 0.91988)

# Constrained optimization
test7 = constrOptim(initPar, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test8 = constrOptim(test7$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
test9 = constrOptim(test8$par, ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)
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
#results of test 18	[1] -1749.773
parValues = test18$par
parValues = c(9.561845e-01,6.832274e-01,9.352358e-04,-1.301391e-02,5.633237e-07,3.251843e-05,7.237832e-02,7.220622e-02,9.201030e-01, 9.202217e-01)

final      = condLogLik(theta  = parValues,
                        series = data[2:725, 5])
dzetaFinal = condLogLikDzeta(theta  = parValues,
                             series = data[2:725, 5])

plot(dzetaFinal[,1], type="l")

#different starting values
#test9 = constrOptim(c(0.5, 0.5, -0.00026460, -0.00026460, 0.0000046171, 0.0000046171, 0.01, 0.01, 0.98, 0.98), ui=constraintMat, ci=constraintVec, condLogLik, series = logRetFirstPeriod, grad=NULL)

#########################################################
####################### Forecasting #####################
#########################################################
### point forecast for period 2011-2012 (726-1183)
#2# reestimation, recursive
#matrix with estimate, error, estimated variance
recursivePar = parDE 
prevPar = parDE
# reestimate the parameters for the longer time series
for (i in 726:1183){
	reestMSGARCH11 = constrOptim(prevPar,
	                              ui     = constraintMat,
	                              ci     = constraintVec,
	                              f      = condLogLik,
	                              series = data[2:i, 5],
	                              grad   = NULL)
	recursivePar = cbind(recursivePar, reestMSGARCH11$par) 
	prevPar = reestMSGARCH11$par
	cat("iteration: ", i, "/n")
	}


# load recursivePar
recursivePar = read.table("reestRecurMSGARCH11.txt", header=TRUE)

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
	info = condLogLikPar(theta = recursivePar[,(i-725)], series = data[2:(i-1), 5]) #calculate dzeta and sigma
	curPar = recursivePar[,(i-725)]	# estimated parameters
	dzeta1Cur = info[[1]][(i-2),1]
	dzeta2Cur = info[[1]][(i-2),2]
	dzeta1[i-725] = dzeta1Cur
	dzeta2[i-725] = dzeta2Cur
	
	p1 = dzeta1Cur * curPar[1] + dzeta2Cur * (1-curPar[2]) # prob of being in state 1 in forecast
	p2 = dzeta2Cur * curPar[2] + dzeta1Cur * (1-curPar[1]) # prob of being in state 2 in forecast
	
	#point estimate
	pointForecastRec[i-725] = p1 * curPar[3] + p2 * curPar[4]
	p1Forecast[i-725] = p1
	p2Forecast[i-725] = p2
	c1Forecast[i-725] = curPar[3]
	c2Forecast[i-725] = curPar[4]
	
	sigma2Cur = dzeta1Cur * (info[[2]][(i-2)])^2 + dzeta2Cur * (info[[3]][(i-2)])^2
	sigma2Forecast[i-725] = p1 * (curPar[5] + (curPar[7] * info[[4]][(i-2)]) + curPar[9]*sigma2Cur) + p2 * (curPar[6] + (curPar[8] * info[[5]][(i-2)]) + curPar[10]*sigma2Cur)
	cat("iteration: ", i, "\n")
	}

plot(c2Forecast,
     type = "l",
     ylim = c(-0.15, 0.2))

par(new=TRUE)
plot(c1Forecast,
     type = "l" ,
     ylim = c(-0.15, 0.2))

par(new=TRUE)
#mean(data$logRet[726:1183])
plot(data$logRet[726:1183],
     type = "l",
     ylim = c(-0.15, 0.2))

estErrorRecur = pointForecastRec - data$logRet[726:1183]
resPlot = cbind(estErrorRecur, sigma2Forecast)
#write.table(resPlot, "resMSGARCH11.txt")


recursiveMSE = (1/length(estErrorRecur)) * sum((estErrorRecur^2))
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
#write.table(conf, "confMSGARCH11.txt")

plot(pointForecastRec,
     type = "l",
     ylim = c(-0.005, 0.005))


# Plot dzeta, prob in high regime and plot of series
# in-sample period
dzeta = condLogLikPar(theta = parDE, series = data[2:725, 5]) #calculate dzeta
dzeta1 = dzeta[[1]][,1]
dzeta2 = dzeta[[1]][,2]
par(mfrow=c(2,1))
jan08 = as.Date("01/01/08", "%d/%m/%y")
dec10 = as.Date("01/01/11", "%d/%m/%y")
DatePlot = as.Date(data$Date[2:725], , "%d/%m/%y")

plot(dzeta1~DatePlot,
     type = "l",
     ylim = c(0.0,1.0),
     lwd  = 1.2,
     xlab = "Date",
     ylab = "Regime probabilities",
     xaxs = "i",
     yaxs = "i",
     xlim = c(jan08, dec10))

plot(data[2:725,5]~DatePlot,
     type = "l",
     lwd  = 1.2,
     xaxs = "i",
     yaxs = "i",
     xlim = c(jan08, dec10),
     xlab = "Date",
     ylab = "Log returns",
     ylim = c(-0.11, 0.11))

abline(h = parDE[3], col= "red", lwd = 1.2)
abline(h = parDE[4], col= "darkgreen", lwd = 1.2) 

prob = cbind(dzeta1, dzeta2)
#write.table(prob, "probMSGARCH11.txt")

# calculate the conditional estimated mean
meanCond = dzeta1 * parDE[3] + dzeta2 * parDE[4]

plot(meanCond,
     type = "l",
     col  = "red")

plot(data[2:725,5],
     type = "l")

plot(data[726:1183,5],
     type = "l",
     ylim = c(-0.2, 0.2))

par(new=TRUE)
plot(confInterval[2:458,1],
     type = "l",
     ylim = c(-0.2, 0.2))

par(new=TRUE)
plot(confInterval[2:458,2],
     type = "l",
     ylim = c(-0.2, 0.2))

par(new=TRUE)
plot(pointForecastRec,
     type = "l",
     ylim = c(-0.2, 0.2),
     col  = "red")


plot(data$logRet[726:1183],
     type = "l")

par(new=TRUE)
plot(estErrorRecur,
     type = "l")

# Density transformation
u = c(rep(0, 458))
for (i in 1:length(u)){
	u[i] = pnorm(data$logRet[i+725],
	              mean = pointForecastRec[i],
	              sd   = sqrt(sigma2Forecast[i]))
	}

hist(u, breaks=20)

# test on uniformity
#Kolmogorov Smirnov
ks.test(u, "punif")
