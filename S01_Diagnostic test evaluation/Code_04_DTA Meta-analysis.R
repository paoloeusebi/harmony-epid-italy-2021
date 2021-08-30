library(tidyverse)
library(mada)
library(runjags)
library(rjags)
library(readxl)


d <- read_excel("S01_Diagnostic test evaluation/data/bastos_serological_covid_2020.xlsx",
                sheet = "elisa")
d

## DTA-MA: perfect reference test

# Forest plot of sensitivity
forest(madad(d),
       type = "sens",
       xlab = "Sensitivity",
       snames= d$Study)


# Forest plot of specificity
forest(madad(d),
       type = "spec",
       xlab = "Specificity",
       snames= d$Study)


# Data points with confidence ellipses on a ROC space
ROCellipse(d, pch = "")
points(fpr(d), sens(d))


#  Let's run the model with reitsma function (mada R package)
fit.reitsma <- reitsma(d)
summary(fit.reitsma)
round(summary(fit.reitsma)[[1]][,c(1,5,6)], 2)

plot(fit.reitsma, cex = 2,
     sroclwd = 2, plotsumm = T,predict = T,pch = 19,
     main = "")
points(fpr(d),
       sens(d), pch = 1)
legend("bottomright",
       c("data points", "summary estimate", "SROC", "95% conf. region", "95% pred.region"),
       pch = c(1, 19, NA, NA, NA),
       lwd = c(NA, 2, 2, 1, 1),
       lty = c(NA, NA, 1,1,3),
       bty = "n")



## Imperfect reference test

# Let's do it with rjags

dta_ma_imp <- "model {

	for(i in 1:l) {
		# Likelihood
		# se, sp are accuracy of CI
		# s2, c2 are accuracy of LU
		# pi is the prevalence

		cell[i,1:4] ~ dmulti(prob[i,1:4],n[i])

		prob[i,1] <- pi[i]*se[i]*s2+(1-pi[i])*(1-sp[i])*(1-c2)
		prob[i,2] <- pi[i]*se[i]*(1-s2)+(1-pi[i])*(1-sp[i])*c2
		prob[i,3] <- pi[i]*(1-se[i])*s2+(1-pi[i])*sp[i]*(1-c2)
		prob[i,4] <- pi[i]*(1-se[i])*(1-s2)+(1-pi[i])*sp[i]*c2

		# Expressing accuracy in terms of HSROC parameters

		b[i] <- exp((beta)/2)
		logit(se[i]) <- (theta[i] + 0.5*alpha[i])/b[i]
		logit(sp[i]) <- -(theta[i] - 0.5*alpha[i])*b[i]

		# Priors for CI accuracy

		theta[i] ~ dnorm(THETA,prec[1])
		alpha[i] ~ dnorm(LAMBDA,prec[2])


		# Priors for prevalence parameters

		pi[i] ~ dbeta(1,1)


	}

	# CI accuracy

	Se_overall <- 1/(1+exp((-THETA-0.5*LAMBDA)/exp(beta/2)))
	Sp_overall <- 1/(1+exp((THETA-0.5*LAMBDA)*exp(beta/2)))

	theta_new ~ dnorm(THETA,prec[1])
	alpha_new ~ dnorm(LAMBDA,prec[2])


	# Predicted values for CI in a new study

	Se_new <- 1/(1+exp(-(theta_new+0.5*alpha_new)/exp(beta/2)))
	Sp_new <- 1/(1+exp((theta_new-0.5*alpha_new)*exp(beta/2)))


	# Priors over the accuracy parameters of CI

	THETA ~ dunif(-2.6,2.6)
	LAMBDA ~ dunif(-5.2,5.2)
	beta ~ dunif(-1.3,1.3)

	for(j in 1:2) {

			prec[j] <- pow(sigma[j],-2)
			sigma[j] ~ dgamma(4,2)
	}

	# Priors over the accuracy parameters of the imperfect reference test (LU)

		###	Informative Priors	###
		s2 ~ dbeta(36.279, 13.761);
		c2 ~ dbeta(67.181, 9.597);

		###	Non-Informative 	###
		#s2 ~ dbeta(1,1);
		#c2 ~ dbeta(1,1);

	#data# l, n, cell
	#monitor# c2, s2, Se_overall ,Sp_overall


}"

# reordering the data
cell <- d %>% select(TP, FP, FN, TN) %>% as.data.frame() %>% as.matrix()
cell

l = length(unique(d$Study)) # number of studies
n = apply(cell, 1, sum) #sample size for each study

# initial values
inits1 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 100022)
inits2 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 300022)
inits3 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 500022)

results <- run.jags(dta_ma_imp,
                    n.chains = 3,
                    burnin = 5000,
                    sample = 10000)

results

plot(results,vars = "Se_overall",
     plot.type = c("trace", "density", "autocorr"))

plot(results,vars = "Sp_overall",
     plot.type = c("trace", "density", "autocorr"))

plot(results,vars = "s2",
     plot.type = c("trace", "density", "autocorr"))

plot(results,vars = "c2",
     plot.type = c("trace", "density", "autocorr"))


# Use another dataset

# add your code here!
