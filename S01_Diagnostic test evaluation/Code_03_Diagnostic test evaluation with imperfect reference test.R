library(runjags)
library(rjags)
testjags()

# 0 = Test Negative; 1 = Test Positive
y <- c(48, 12, 4, 36); y
n <- sum(y); n

t <- matrix(
  y,
  nrow = 2,
  byrow = F,
  dimnames = list(0:1,
                  0:1))
t

# HW Model: 2 tests, 1 population
hw_2t_1p <- c("model{
  y ~ dmulti(prob, n)

  # Test1- Test2-
	prob[1] <- (p*((1-se[1])*(1-se[2]))) + ((1-p)*((sp[1])*(sp[2])))
  # Test1+ Test2-
	prob[2] <- (p*((se[1])*(1-se[2]))) + ((1-p)*((1-sp[1])*(sp[2])))
  # Test1- Test2+
	prob[3] <- (p*((1-se[1])*(se[2]))) + ((1-p)*((sp[1])*(1-sp[2])))
	# Test1+ Test2+
	prob[4] <- (p*((se[1])*(se[2]))) + ((1-p)*((1-sp[1])*(1-sp[2])))

	", "

  p ~ dbeta(1, 1)
  se[1] ~ dbeta(1, 1)
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)
  sp[2] ~ dbeta(1, 1)
  #data# y, n
  #monitor# p, se, sp
  #inits# p, se, sp
}
")

# choice of initial values (uncomment below if you want to try)
# p <- list(chain1=0.05, chain2=0.95)
# se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
# sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))

# initial values for 2 chains
inits1 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 100022)
inits2 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 300022)


# run
results <- run.jags(hw_2t_1p,
                    n.chains = 2,
                    burnin = 5000,
                    sample = 15000)

print(results)

plot(results,
     vars = list("se", "sp", "p"),
     layout = c(3, 3),
     plot.type = c("trace", "histogram", "autocorr"))

# Analysing simulated data is useful to check if
# we can recover parameter values.
se1 <- 0.9; sp1 <- 0.95
se2 <- 0.8; sp2 <- 0.99
prevalence <- 0.9; n <- 1000
truestatus <- rbinom(n, 1, prevalence)
Test1 <- rbinom(n, 1, (truestatus * se1) + ((1-truestatus) * (1-sp1)))
Test2 <- rbinom(n, 1, (truestatus * se2) + ((1-truestatus) * (1-sp2)))
t <- table(Test1, Test2); t
y <- as.numeric(t); y
n <- sum(y); n

results <- run.jags(hw_2t_1p,
                    n.chains =2,
                    inits = list(inits1, inits2),
                    burnin = 5000,
                    sample = 15000)
print(results)

# How well do we recover our parameters?


# Imposing constraints on the priors
hw_2t_1p_priors <- "model{
  y ~ dmulti(prob, n)

  # Test1- Test2-
	prob[1] <- (p * ((1-se[1])*(1-se[2]))) + ((1-p) * ((sp[1])*(sp[2])))
  # Test1+ Test2-
	prob[2] <- (p * ((se[1])*(1-se[2]))) + ((1-p) * ((1-sp[1])*(sp[2])))
  # Test1- Test2+
	prob[3] <- (p * ((1-se[1])*(se[2]))) + ((1-p) * ((sp[1])*(1-sp[2])))
  # Test1+ Test2+
	prob[4] <- (p * ((se[1])*(se[2]))) + ((1-p) * ((1-sp[1])*(1-sp[2])))

  p ~ dbeta(1, 1)
  se[1] ~ dbeta(HPSe[1,1], HPSe[1,2])T(1-sp[1], )
  sp[1] ~ dbeta(HPSp[1,1], HPSp[1,2])
  se[2] ~ dbeta(HPSe[2,1], HPSe[2,2])T(1-sp[2], )
  sp[2] ~ dbeta(HPSp[2,1], HPSp[2,2])
  #data# y, n, HPSe, HPSp
  #monitor# p, se, sp
  #inits# p, se, sp
}
"

# Note that we specify the prior hyperparameters as data so we can change these from R without havíng to edit the model file (this is optional!)

# Simulate again
se1 <- 0.9; sp1 <- 0.95
se2 <- 0.8; sp2 <- 0.99
prevalence <- 0.5
# Change N to be 10, 100 or 1000:
n <- 1000
truestatus <- rbinom(n, 1, prevalence)
Test1 <- rbinom(n, 1, (truestatus * se1) + ((1-truestatus) * (1-sp1)))
Test2 <- rbinom(n, 1, (truestatus * se2) + ((1-truestatus) * (1-sp2)))
t <- table(Test1, Test2); t
y <- as.numeric(t); y
n <- sum(y); n

HPSe <- matrix(c(1,1,1,1), nrow=2, ncol=2)
HPSp <- matrix(c(1,1,1,1), nrow=2, ncol=2)

# Initial values for two chains
p <- list(chain1=0.05, chain2=0.95)
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
results <- run.jags(hw_2t_1p_priors,
                    n.chains=2,
                    burnin = 5000,
                    sample = 10000)
results

# How well do we recover our parameters?

# H-W Model: 2 tests, m pops
hw_2t_mpops <- c("model{

                for (i in 1:m) {

                # likelihood

                y[i,1:4] ~ dmulti(prob[i,1:4],n[i])

                  # Test1- Test2-
                  prob[i,1] <- (pi[i] * ((1-se[1])*(1-se[2]))) + ((1-pi[i]) * ((sp[1])*(sp[2])))
                  # Test1+ Test2-
                  prob[i,2] <- (pi[i] * ((se[1])*(1-se[2]))) + ((1-pi[i]) * ((1-sp[1])*(sp[2])))
                  # Test1- Test2+
                  prob[i,3] <- (pi[i] * ((1-se[1])*(se[2]))) + ((1-pi[i]) * ((sp[1])*(1-sp[2])))
                  # Test1+ Test2+
                  prob[i,4] <- (pi[i] * ((se[1])*(se[2]))) + ((1-pi[i]) * ((1-sp[1])*(1-sp[2])))

                # priors for prevalence parameters
                pi[i] ~ dbeta(1,1)
                }


	", "
  se[1] ~ dbeta(1, 1)T(1-sp[1], )
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)T(1-sp[2], )
  sp[2] ~ dbeta(1, 1)

  #data# m, n, y
  #inits#
  #monitor# se, sp, pi
}
")

# Simulate two pops
se1 <- 0.9; sp1 <- 0.95 # Test 1
se2 <- 0.8; sp2 <- 0.99 # Test 2
pr1 <- 0.25; pr2 <- 0.5 # prevalences
n1 <- 10000; n2 <- 20000 # sample size

truestatus1 <- rbinom(n1, 1, pr1)
truestatus2 <- rbinom(n2, 1, pr2)
truestatus <- c(truestatus1, truestatus2)
pop <- c(rep(1,n1), rep(2,n2))

Test1 <- rbinom(n1+n2, 1, (truestatus * se1) + ((1-truestatus) * (1-sp1)))
Test2 <- rbinom(n1+n2, 1, (truestatus * se2) + ((1-truestatus) * (1-sp2)))

df <- data.frame(pop, Test1, Test2)

t <- table(df$Test1, df$Test2, df$pop); t
y_ <- as.numeric(t); y_

y <- matrix(
  y_,
  nrow = 2,
  byrow = T)
y

m = 2 # number of populations
n = apply(y, 1, sum); n

# initial values for two chains
inits1 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 100022)
inits2 = list(".RNG.name" ="base::Mersenne-Twister", ".RNG.seed" = 300022)

results <- run.jags(hw_2t_mpops,
                    n.chains=2,
                    inits = list(inits1, inits2),
                    burnin = 10000,
                    sample = 100000)
results
plot(results,
     vars = list("se", "sp", "pi"),
     layout = c(3, 6),
     plot.type = c("trace", "histogram", "autocorr"))

# How well do we recover our parameters?
