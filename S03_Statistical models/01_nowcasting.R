#
# Nowcasting COVID-19 incidence indicators ----
#

source("RichardsGrowthGLM.R")
source("data.R")
load_install_package("ggplot2")

# New Positives ----

# Starting and ending date of analysis
start = "2020-02-25"; end = "2020-07-19"

dati = copy(dati_ita)
dati = dati[date >= start,]
dati = dati[date <= end,]
# response variable
dati[, y := nuovi_positivi]
# time variable
dati[, t := as.numeric(date - as.Date(start))+1]
# dummy variable for days after sundays and holidays
dati[, WeekEnd := (weekdays(date) %in% "Monday" | 
                   date %in% (Italian_holidays+1))]

# View(dati[, c("date", "t", "y", "WeekEnd")])

ggplot(dati, aes(x = date, y = y)) +
  geom_point(size = 1) +
  labs(x = NULL, y = "New positives") +
  theme_minimal()

# fit model
mod = growthGLM(di = dati$y,                  # response variable
                ti = dati$t,                  # time variable
                X = cbind(1, dati$WeekEnd),   # design matrix
                family = "Negative Binomial", # distribution
                tPred = max(dati$t)+14,       # prediction time horizon
                alpha = 0.05)                 # confidence level for predictions
str(mod)

mod$pars   # MLEs
mod$lik    # log-likelihood at MLEs
mod$R2diff # pseudo-Rsquared for first differences counts
mod$R2cum  # pseudo-Rsquared for cumulative counts

# table of parameters
tab = data.frame("Parameters" = mod$pars, "SE" = mod$se)
rownames(tab) = c("log(h)", "p", "s", "log(a)", "WeekEnd", "log(r)", "log(v)")
tab

# table of parameters (1-alpha)100% confidence intervals
alpha = 0.05
q = qnorm(1-alpha/2)
tab2 = data.frame("Parameters" = tab$Parameters,
                  "lower 95% CI" = tab$Parameters - q*tab$SE,
                  "upper 95% CI" = tab$Parameters + q*tab$SE,
                  check.names = FALSE)
# transform parameters to natural scale
tab2[c(1,4,6,7),] = exp(tab2[c(1,4,6,7),])
rownames(tab2) = c("h", "p", "s", "a", "WeekEnd", "r", "v")
tab2

peak <- floor(mod$pars[2] + log10(mod$pars[3])/exp(mod$pars[1]))
peak
dati$date[peak]



df = data.frame(date = c(dati$date, max(dati$date)+1:14),
                y = c(dati$y, rep(NA, 14)), 
                Y = c(cumsum(dati$y), rep(NA,14)),
                yfit = mod$linPredDiff,
                ylow = mod$lowdiff,
                yup = mod$updiff,
                Yfit = mod$linPredCum,
                Ylow = mod$lowcum,
                Yup = mod$upcum)

ggplot(df, aes(x = date, y = y)) +
  geom_point(size = 1) +
  geom_line(aes(y = yfit), col = "firebrick") +
  geom_ribbon(aes(ymin = ylow, ymax = yup), 
              fill = "firebrick", alpha = 0.3) +  
  geom_vline(xintercept = dati$date[peak,], lty = 2) +
  labs(x = NULL, y = "New positives") +
  theme_minimal()

ggplot(df, aes(x = date, y = Y)) +
  geom_point(size = 1) +
  geom_line(aes(y = Yfit), col = "firebrick") +
  geom_ribbon(aes(ymin = Ylow, ymax = Yup), 
              fill = "firebrick", alpha = 0.3) +  
  labs(x = NULL, y = "Cum Positives") +
  theme_minimal()


# New Deceased ----

# Starting and ending date of analysis
start = "2020-02-25"; end = "2020-07-19"
# start = "2020-08-16"; end = "2021-01-01"

dati = copy(dati_ita)
dati[, "nuovi_deceduti" := pmax(0, deceduti - shift(deceduti))]
dati = dati[date >= start,]
dati = dati[date <= end,]
# response variable
dati[, y := nuovi_deceduti]
# time variable
dati[, t := as.numeric(date - as.Date(start))+1]

# View(dati[, c("date", "t", "y")])

ggplot(dati, aes(x = date, y = y)) +
  geom_point(size = 1) +
  labs(x = NULL, y = "New deceased") +
  theme_minimal()

mod = growthGLM(di = dati$y,  # response variable
                ti = dati$t,  # time 
                X  = cbind(rep(1,length(dati$y))),  # design matrix
                family = "Negative Binomial",       # distribution
                tPred = max(dati$t)+14,  # prediction time horizon
                alpha = 0.05)            # confidence level for predictions
str(mod)

mod$pars   # MLEs
mod$lik    # log-likelihood at MLEs
mod$R2diff # pseudo-Rsquared for first differences counts
mod$R2cum  # pseudo-Rsquared for cumulative counts

# table of parameters
tab = data.frame("Parameters" = mod$pars, "SE" = mod$se)
rownames(tab) = c("log(h)", "p", "s", "log(a)", "log(r)", "log(v)")
tab

# table of parameters (1-alpha)100% confidence intervals
alpha = 0.05
q = qnorm(1-alpha/2)
tab2 = data.frame("Parameters" = tab$Parameters,
                  "lower 95% CI" = tab$Parameters - q*tab$SE,
                  "upper 95% CI" = tab$Parameters + q*tab$SE,
                  check.names = FALSE)
# transform parameters to natural scale
tab2[c(1,4,5,6),] = exp(tab2[c(1,4,5,6),])
rownames(tab2) = c("h", "p", "s", "a", "r", "v")
tab2

peak <- floor(mod$pars[2] + log10(mod$pars[3])/exp(mod$pars[1]))
peak
dati$date[peak]

df = data.frame(date = c(dati$date, max(dati$date)+1:14),
                y = c(dati$y, rep(NA, 14)), 
                Y = c(cumsum(dati$y), rep(NA,14)),
                yfit = mod$linPredDiff,
                ylow = mod$lowdiff,
                yup = mod$updiff,
                Yfit = mod$linPredCum,
                Ylow = mod$lowcum,
                Yup = mod$upcum)

ggplot(df, aes(x = date, y = y)) +
  geom_point(size = 1) +
  geom_line(size = 0.1) +
  geom_line(aes(y = yfit), col = "firebrick") +
  geom_ribbon(aes(ymin = ylow, ymax = yup),
              fill = "firebrick", alpha = 0.3) +
  geom_vline(xintercept = dati$date[peak], lty = 2) +
  labs(x = NULL, y = "New deceased") +
  theme_minimal()

ggplot(df, aes(x = date, y = Y)) +
  geom_point(size = 1) +
  geom_line(aes(y = Yfit), col = "firebrick") +
  geom_ribbon(aes(ymin = Ylow, ymax = Yup),
              fill = "firebrick", alpha = 0.3) +
  labs(x = NULL, y = "Cum deceased") +
  theme_minimal()



