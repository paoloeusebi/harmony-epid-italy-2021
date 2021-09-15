library(runjags)
library(rjags)
testjags()


# Bayesian model
basicjags <- "model{
  # Likelihood part:
  Positives ~ dbinom(prevalence, TotalTests)

  # Prior part:
  prevalence ~ dbeta(2, 2)

  # Hooks for automatic integration with R:
  #data# Positives, TotalTests
  #monitor# prevalence
  #inits# prevalence
}
"

Positives <- 7
TotalTests <- 10

# initial values to be retrieved by runjags:
prevalence <- list(chain1=0.05, chain2=0.95)

# Running
results <- run.jags(basicjags,
                    n.chains = 2,
                    burnin = 1000,
                    sample = 5000)
print(results)
plot(results)

