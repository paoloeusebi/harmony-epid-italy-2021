library(runjags)
library(rjags)
testjags()

# Data from Whitman et al. Test performance evaluation of SARS-CoV-2 serological assays.
# medRxiv. 2020;(https://doi.org/10.1101/2020.04.25.20074856). Accessed July 4, 2020

y <- c(17, 4, 14, 38)
t <- matrix(
  y,
  nrow = 2,
  byrow = F,
  dimnames = list(c("T+", "T-"),
                  c("D+", "D-")))
t

# Sample size
sum(t)

# Prevalence
p <- (t[1]+t[2])/(sum(t))
p

# Sensitivity
Se <- t[1]/(t[1]+t[2])
Se

# Specificity
Sp <- t[4]/(t[3]+t[4])
Sp

# Predictive values
PPV <- t[1]/(t[1]+t[3])
PPV

NPV <- t[4]/(t[2]+t[4])
NPV

# Likelihood ratios
PLR <- Se/(1-Sp)
PLR
NLR <- (1-Se)/Sp
NLR

# Diagnostic odds ratio
DOR <- PLR/NLR
DOR

# Overall accuracy
sum(diag(t))/sum(t)


# Bayesian model
model <- " model {

# likelihood

  y[1:4] ~ dmulti(prob[1:4], n)

  prob[1] <- p * Se
  prob[2] <- p * (1 - Se)
  prob[3] <- (1 - p) * (1 - Sp)
  prob[4] <- (1 - p) * Sp

# priors

  p ~ dbeta(1, 1)
  Se ~ dbeta(1,1)
  Sp ~ dbeta(1,1)


#data# n, y
#inits#
#monitor# Se, Sp, p

}

"

y # data
n <- sum(y) # total sample size
# initial values for the three chains
inits1 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 100022)
inits2 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 300022)
inits3 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 500022)

# Running
results <- run.jags(model,
                    n.chains = 3,
                    inits = list(inits1, inits2, inits3),
                    burnin = 1000,
                    sample = 5000)

print(results)

plot(results,
     vars = list("Se", "Sp", "p"),
     layout = c(3, 3),
     plot.type = c("trace", "histogram", "autocorr"))

# Bayesian model: calculating NPV according to specific prevalence

model <- "model {

# likelihood

  y[1:4] ~ dmulti(prob[1:4], n)

  prob[1] <- p * Se
  prob[2] <- p * (1 - Se)
  prob[3] <- (1 - p) * (1 - Sp)
  prob[4] <- (1 - p) * Sp

# priors

  p ~ dbeta(1, 1)
  Se ~ dbeta(1,1)
  Sp ~ dbeta(1,1)

# Computing NPV according to certain prevalence threshold

  NPV = (Sp*p)/((1-Se)*p + Sp*(1-p))
  NPV1 = (Sp*p1)/((1-Se)*p1 + Sp*(1-p1))

#data# n, y, p1
#inits#
#monitor# Se, Sp, p, NPV, NPV1

}

"
# The specific prevalence is provided as data for the model
p1 <- 0.005 # change from 0 to 1

# initial values
inits1 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 100022)
inits2 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 300022)
inits3 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 500022)

# Running
results <- run.jags(model,
                    n.chains = 3,
                    inits = list(inits1, inits2, inits3),
                    burnin = 1000,
                    sample = 5000)
print(results)
plot(results,
     vars = list("NPV1", "NPV"),
     layout = c(2, 1),
     plot.type = "histogram")

