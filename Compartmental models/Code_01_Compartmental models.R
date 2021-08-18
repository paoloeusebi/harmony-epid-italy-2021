#Load required libraries
library(tidyverse)
library(deSolve)
library(reshape2)

# Model inputs

initial_state_values=c(S=999999,I=1,R=0)
parameters=c(gamma=0.20,beta=0.5)

# Time points

time=seq(from=1,to=100,by=1)

# SIR model function

sir_model <- function(time,state,parameters){
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


#Solving the differential equations
output <- ode(
  y = initial_state_values,
  func = sir_model,
  parms = parameters,
  times = time
  ) %>%
  as.data.frame()

output2 <- output %>%
  gather(key = compartment, value = count, S:R)

# To plot the proportion of susceptible, infected and recovered individuals over time
ggplot(data = output2,
       aes(x = time, y = count/1000000, fill = compartment)) +
  geom_area(col="white") +
  scale_fill_manual(values = c("firebrick", "deepskyblue4", "yellowgreen")) +
  xlab("Time (days)") +
  ylab("Proportion of the population") +
  scale_color_discrete(name = "State")
