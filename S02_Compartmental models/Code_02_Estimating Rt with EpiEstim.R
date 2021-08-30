library(EpiEstim)
library(tidyverse)

# from https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html

## load data
data(Flu2009)
## incidence:
head(Flu2009$incidence)
## serial interval (SI) distribution:
Flu2009$si_distr
## interval-ceonsored serial interval data:
## each line represents a transmission event,
## EL/ER show the lower/upper bound of the symptoms onset date in the infector
## SL/SR show the same for the secondary case
## type has entries 0 corresponding to doubly interval-censored data
## (see Reich et al. Statist. Med. 2009).
head(Flu2009$si_data)

# We can use the incidence R package to easily plot the daily incidence data:
library(incidence)
plot(as.incidence(Flu2009$incidence$I, dates = Flu2009$incidence$dates))

# Calculate Rt
# In this example, we only specify the mean and standard deviation of the serial interval. In that case an offset gamma distribution is used for the serial interval. In the following example, we use the mean and standard deviation of the serial interval for flu from Ferguson et al., Nature, 2005:
res_parametric_si <- estimate_R(Flu2009$incidence,
                                method = "parametric_si",
                                config = make_config(list(mean_si = 2.6,
                                                          std_si = 1.5)))
head(res_parametric_si$R)
plot(res_parametric_si, "R")
plot(res_parametric_si, legend = FALSE)



# Estimating R with a non parametric serial interval distribution
res_non_parametric_si <- estimate_R(Flu2009$incidence,
                                    method="non_parametric_si",
                                    config = make_config(list(
                                      si_distr = Flu2009$si_distr)))
plot(res_non_parametric_si, "R")
