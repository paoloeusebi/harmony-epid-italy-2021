#Load required libraries
library(tidyverse)
library(deSolve)

# Time points

time <- seq(from=1, to=100, by=1)

# SIR Model ---------------------------------------------------------------

## model inputs
initial_state_values=c(S=999999,I=1,R=0)
parameters=c(gamma=0.20,beta=0.5)

## SIR model function
sir_model <- function(time, state, parameters){
  with(as.list(c(state,parameters)),{
    N=S+I+R
    lambda=beta*(I/N)
    dS=-lambda*S
    dI=lambda*S-gamma*I
    dR=gamma*I

    return(list(c(dS,dI,dR)))
  }
  )
}

## solving the differential equations
output <- ode(
  y = initial_state_values,
  func = sir_model,
  parms = parameters,
  times = time
  ) %>%
  as.data.frame()

output

output2 <- output %>%
  gather(key = compartment, value = count, S:R)

output2

## plot the proportion of states (S, I, R) over time
ggplot(data = output2,
       aes(x = time, y = count/sum(initial_state_values), fill = compartment)) +
  geom_area(col="white") +
  scale_fill_manual(values = c("firebrick", "deepskyblue4", "yellowgreen")) +
  xlab("Time (days)") +
  ylab("Proportion of the population") +
  scale_color_discrete(name = "State")


# SEIR Model --------------------------------------------------------------
## model inputs
initial_state_values=c(S=500, E=499, I=1, R=0)
parameters=c(gamma=0.20, beta=0.5, R0=2, sigma=1)

## SEIR model function
seir_model <- function(time, state, parameters){
  with(as.list(c(state,parameters)),{
    dS <- -(R0*gamma) * S * I/(S + E + I + R)
    dE <- (R0*gamma) * S * I/(S + E + I + R) - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I

    return(list(c(dS, dE, dI, dR)))
  }
  )
}

## solving the differential equations
output <- ode(
  y = initial_state_values,
  func = seir_model,
  parms = parameters,
  times = time
  ) %>%
  as.data.frame()

output

output2 <- output %>%
  gather(key = compartment, value = count, S:R)

output2

## plot the proportion of states (S, E, I, R) over time
ggplot(data = output2,
       aes(x = time, y = count/sum(initial_state_values), fill = compartment)) +
  geom_area(col="white") +
  scale_fill_manual(values = c("violet","firebrick", "deepskyblue4", "yellowgreen")) +
  xlab("Time (days)") +
  ylab("Proportion of the population") +
  scale_color_discrete(name = "State")


# Exercise ----------------------------------------------------------------

# Run previous models changing the parameters
# Simulate a reasonable model (SIR and/or SEIR) for an epidemics of interest


