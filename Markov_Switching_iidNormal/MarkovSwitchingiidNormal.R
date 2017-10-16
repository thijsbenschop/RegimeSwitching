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

# Function to evaluate the loglikelihood, 2 states, both iid normal
condLogLikPar = function(theta, series){
	# calculates the conditional loglikelihood for a Markov Switching model 
	# with 2 states
	# Arguments
	# p11 transition probability to stay in state 1
	# p22 transition probability to stay in state 2
	# mu1, mu2, mean of distribution in resp. state 1 and 2
	# sigma1, sigma2, variance of distribution in resp. state 1 and 2
	# returns : value of conditional log-likelihood
	p11 = theta[1]
	p22 = theta[2]
	mu1 = theta[3]
	mu2 = theta[4] 
	sigma1 = theta[5]
	sigma2 = theta[6]
	
	n 		= length(series)	# number of obs.
	p 		= cbind(c(p11, 1-p22), c(1-p11, p22)) # transition matrix
	dzeta 	= cbind(rep(0, n), rep(0, n))
	f 		= cbind(rep(0, n))
	eta 	= cbind(rep(0,n), rep(0,n))
	
	dzetaInit = c((1-p[2,2])/(2-p[1,1]-p[2,2]), (1-p[1,1])/(2-p[2,2]-p[1,1]))
	# startvalue for iterations, assuming the Markov chain is ergodic
	# alternative  dzetaInit = c(0.5, 0.5)	
	
	for (i in 1:n){
		# Evaluate the densities under the two regimes
		############### create if else for different functional forms ##############
		eta[i, 1] = dnorm(x=series[i], mean=mu1, sd=sigma1)
		eta[i, 2] = dnorm(x=series[i], mean=mu2, sd=sigma2)
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
	}
	logf = sum(log(f[5:n]))
	cat(logf, " ")
	flush.console()
	output = cbind(eta, f, dzeta)
	return(list(dzeta))
}

# Function to evaluate the loglikelihood, 2 states, both iid normal
condLogLik = function(theta, series){
	# calculates the conditional loglikelihood for a Markov Switching model 
	# with 2 states
	# Arguments
	# p11 transition probability to stay in state 1
	# p22 transition probability to stay in state 2
	# mu1, mu2, mean of distribution in resp. state 1 and 2
	# sigma1, sigma2, variance of distribution in resp. state 1 and 2
	# returns : value of conditional log-likelihood
  p11 = theta[1]
  p22 = theta[2]
  mu1 = theta[3]
  mu2 = theta[4] 
  sigma1 = theta[5]
  sigma2 = theta[6]
	
  n 		= length(series)	# number of obs.
  p 		= cbind(c(p11, 1-p22), c(1-p11, p22)) # transition matrix
  dzeta 	= cbind(rep(0, n), rep(0, n))
  f 		= cbind(rep(0, n))
  eta 	    = cbind(rep(0,n), rep(0,n))
	
	dzetaInit = c((1-p[2,2])/(2-p[1,1]-p[2,2]), (1-p[1,1])/(2-p[2,2]-p[1,1]))
	# startvalue for iterations, assuming the Markov chain is ergodic
	# alternative  dzetaInit = c(0.5, 0.5)	
	
	for (i in 1:n){
		# Evaluate the densities under the two regimes
		############### create if else for different functional forms ##############
		eta[i, 1] = dnorm(x=series[i], mean=mu1, sd=sigma1)
		eta[i, 2] = dnorm(x=series[i], mean=mu2, sd=sigma2)
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
	}
	logf = sum(log(f[5:n]))
	if(is.nan(logf)==TRUE){
		logf = 99999999
		cat("Error: not a number", "/n")
		flush.console()
		}
	#cat(logf, " ")
	#flush.console()
	output = cbind(eta, f, dzeta)
	return(-logf)
}

# Read-in data
data = read.table("dataCO2.txt", header = TRUE)

logRetFirstPeriod = data[2:725,5]

## Other optimization alghorithms
#			  p11, p22,  mu1,  mu2,  sigma1, sigma2
lowParDE 	= c(0,   0,   -0.3, -0.3,   0,     0 )
upParDE 	= c(1,   1,    0.3,  0.3, 0.3,   0.3 )
controlDE 	= list(NP = 60, itermax = 2500)
testDE = DEoptim(fn 		= condLogLik,
				 lower 		= lowParDE,
				 upper 		= upParDE,
				 series 	= logRetFirstPeriod,
				 control 	= controlDE)
# Result
#Iteration: 2500 bestvalit: -1719.999284 bestmemit:    0.985531    0.976359    0.001405   -0.003738    0.016150    0.033617
#parDE = c(0.982090, 0.992281, -0.004157, 0.000865, 0.000289, 0.000052)
parDE = c(0.985531, 0.976359, 0.001405, -0.003738, 0.016150, 0.033617)


# Plot dzeta, prob in high regime and plot of series
# in-sample period
dzeta 		= condLogLikPar(theta 	= parDE,
							series 	= data[2:725, 5]) #calculate dzeta
dzeta1 		= dzeta[[1]][,1]
dzeta2 		= dzeta[[1]][,2]
par(mfrow 	= c(2,1))
jan08 		= as.Date("01/01/08", "%d/%m/%y")
dec10 		= as.Date("01/01/11", "%d/%m/%y")
DatePlot 	= as.Date(data$Date[6:725], ,"%d/%m/%y")

plot(dzeta1[5:724]~DatePlot,
	 type 	= "l",
	 ylim 	= c(0.0,1.0),
	 lwd 	= 2,
	 xlab 	= "Date",
	 ylab 	= "Regime probabilities",
	 xaxs 	= "i",
	 yaxs 	= "i",
	 xlim 	= c(jan08, dec10))

plot(data[6:725,5]~DatePlot,
	 type 	= "l",
	 lwd 	= 1.2,
	 xaxs 	= "i",
	 yaxs 	= "i",
	 xlim 	= c(jan08, dec10),
	 xlab 	= "Date",
	 ylab 	= "Log returns",
	 ylim 	= c(-0.11, 0.11))


llh = condLogLik(parDE, logRetFirstPeriod)

testnlm = nlm(f 		= condLogLik,
			  p 		= c(0.7, 0.8, mean(logRetFirstPeriod), mean(logRetFirstPeriod), sd(logRetFirstPeriod), sd(logRetFirstPeriod)),
			  series 	= logRetFirstPeriod,
			  iterlim 	= 100)

testOptim = optim(f 		= condLogLik,
				  p 		= c(0.7, 0.8, mean(logRetFirstPeriod), mean(logRetFirstPeriod), sd(logRetFirstPeriod), sd(logRetFirstPeriod)),
				  series 	= logRetFirstPeriod)

parValues = testnlm$estimate

# loop testing different starting values for p11 and p22
evalMat1 = matrix(rep(0, 121), 11, 11)
evalMat2 = matrix(rep(0, 121), 11, 11)
evalMat3 = matrix(rep(0, 121), 11, 11)
evalMat4 = matrix(rep(0, 121), 11, 11)
evalMat5 = matrix(rep(0, 121), 11, 11)
evalMat6 = matrix(rep(0, 121), 11, 11)
evalMat7 = matrix(rep(0, 121), 11, 11)
evalMat8 = matrix(rep(0, 121), 11, 11)

# loop testing different starting values for p11 and p22 

for (i in 0:10){
	for (j in 0:10){
		reestimate =nlm(f 		= condLogLik,
						p 		= c(0.7, 0.8, mean(logRetFirstPeriod), mean(logRetFirstPeriod), sd(logRetFirstPeriod), sd(logRetFirstPeriod)),
						series 	= logRetFirstPeriod,
						iterlim	= 100)

		evalMat1[i+1, j+1] = reestimate$estimate[1]
		evalMat2[i+1, j+1] = reestimate$estimate[2]
		evalMat3[i+1, j+1] = reestimate$estimate[3]
		evalMat4[i+1, j+1] = reestimate$estimate[4]
		evalMat5[i+1, j+1] = reestimate$estimate[5]
		evalMat6[i+1, j+1] = reestimate$estimate[6]
		evalMat7[i+1, j+1] = reestimate$minimum
		evalMat8[i+1, j+1] = reestimate$code	
		}
	}

AIC = 2*(test2$minimum) + 2*6

#########################################################
############## Forecasting ##############################
#########################################################

#2# reestimation, recursive
#matrix with estimate, error, estimated variance

recursivePar 	= parValues 
prevPar 		= parValues 
# reestimate the parameters for the longer time series
for (i in 726:1183){
	reestMSNormal 	= nlm(f 		= condLogLik,
						  p 		= prevPar,
						  series 	= data[(2:i),5],
						  iterlim 	= 100)

	recursivePar 	= as.data.frame(cbind(recursivePar, reestMSNormal$estimate)) 
	prevPar 		= reestMSNormal$estimate
	cat("iteration: ", i, "\n")
	}

# load recursivePar
recursivePar = read.table("reestRecurMSNormal.txt", header=TRUE)

llh = condLogLik(theta 	= c(recursivePar[,1]),
                series 	= data[2:725, 5])

llhcor 	= llh/724*720
AIC 	= 2*llhcor + 2*6

# calculate point and density estimates, confidence intervals
pointForecastRec 	= rep(0, 458)
sigma2Forecast 		= rep(0, 458)

p1Forecast = rep(0, 458)
p2Forecast = rep(0, 458)
c1Forecast = rep(0, 458)
c2Forecast = rep(0, 458)

for (i in 726:1183){	
	# Forecasting values for i
	# calculate dzeta[i], sigma[i]
	info 	= condLogLikPar(theta 	= recursivePar[,(i-725)],
							series 	= data[2:(i-1), 5]) #calculate dzeta and sigma
	curPar 	= recursivePar[,(i-725)]	# estimated parameters 
	# [1] p11,#[2] p22,#[3] mu1,#[4] mu2 ,#[5] sigma1,#[6] sigma2
	dzeta1Cur = info[[1]][(i-2),1]
	dzeta2Cur = info[[1]][(i-2),2]
	p1 = dzeta1Cur * curPar[1] + dzeta2Cur * (1-curPar[2]) # prob of being in state 1 in forecast
	p2 = dzeta2Cur * curPar[2] + dzeta1Cur * (1-curPar[1]) # prob of being in state 2 in forecast
	#point estimate
	mean1 = curPar[3]
	mean2 = curPar[4]
	
	pointForecastRec[i-725] = p1 * mean1 + p2 * mean2
	p1Forecast[i-725] 		= p1
	p2Forecast[i-725] 		= p2
	c1Forecast[i-725] 		= mean1
	c2Forecast[i-725] 		= mean2
	sigma2Forecast[i-725] 	= p1 * curPar[5]^2 + p2 * curPar[6]^2 
	}

##
plot(c2Forecast,
	 type = "l",
	 ylim = c(-0.1, 0.2))

par(new=TRUE)
plot(c1Forecast,
	 type = "l",
	 ylim = c(-0.1, 0.2),
     xaxt = "n",
     yaxt = "n",
     ann  = FALSE)

par(new=TRUE)
plot(data$logRet[726:1183],
	type = "l",
	ylim = c(-0.1, 0.2,),
	xaxt = "n",
	yaxt = "n",
	ann  = FALSE)

par(new=TRUE)
plot(pointForecastRec,
	 type 	= "l",
	 ylim 	= c(-0.1, 0.2),
	 col 	= "red",
     xaxt = "n",
     yaxt = "n",
     ann  = FALSE)

abline(h=mean(data$logRet[726:1183]))


estErrorRecur 	= pointForecastRec[1:458] - data$logRet[726:1183]
recursiveMSE 	= (1/length(estErrorRecur)) * sum(estErrorRecur^2)
recursiveMAE 	= (1/length(estErrorRecur)) * sum(abs(estErrorRecur))
recursiveMSE
recursiveMAE

resPlot = cbind(estErrorRecur, sigma2Forecast)

# logreturns and predicted 95 percent confidence intervals
# calculate 2,5 and 97,5 quantiles
confInterval = cbind(rep(0, 458), rep(0, 458)) 
for (i in 1:458){
	confInterval[i,1] = qnorm(0.025, mean=pointForecastRec[i], sd=sqrt(sigma2Forecast[i]))
	confInterval[i,2] = qnorm(0.975, mean=pointForecastRec[i], sd=sqrt(sigma2Forecast[i]))
	}

conf = cbind(confInterval, pointForecastRec)

# Plot series within estimated confidence intervals
plot(confInterval[,1],
	 type = "l",
	 ylim = c(-0.2, 0.2))

par(new=TRUE) 
plot(confInterval[,2],
	 type = "l",
	 ylim = c(-0.2, 0.2))

par(new=TRUE)
plot(data$logRet[726:1183],
	 type = "l",
	 ylim = c(-0.2, 0.2))

# Plot dzeta, prob in high regime and plot of series
par(mfrow=c(2,1))
plot(p1Forecast,
	 type = "l",
	 ylim = c(0,1))

plot(pointForecastRec,
	 type = "l")

##
plot(pointForecastRec,
	 type 	= "l",
	 ylim 	= c(-0.3, 0.3),
	 col 	= "red")

par(new=TRUE)
plot(data$logRet[726:1183],
	 type = "l",
	 ylim = c(-0.3, 0.3))

##
plot(pointForecastRec,
     type 	= "l",
     ylim 	= c(-0.3, 0.3),
     col 	= "red")

par(new=TRUE)
plot(estErrorRecur,
	 type = "l",
	 ylim = c(-0.3, 0.3))

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