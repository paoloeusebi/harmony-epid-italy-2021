#
# COVINDEX based on a GAM Beta regression model ----
#

source("CovindexGamBetaReg.R")
source("data.R")
load_install_package("ggplot2")
load_install_package("scales")

# Starting and ending date of analysis
start = "2020-03-01"; end = "2021-06-30"

dati = copy(dati_ita)
dati = dati[date >= start,]
dati = dati[date <= end,]
# create variables for COVINDEX analysis
dati[, TPR := positives/swabs]
# time variable
dati[, t := as.numeric(date - as.Date("2020-01-01"))]
# dati[, t := as.numeric(date - as.Date(start))+1]
# dummy variable for days after sundays and holidays 
dati[, WeekEnd := (weekdays(date) %in% c("Sunday", "Monday") | 
                   date %in% (Italian_holidays+1))]

# View(dati[, c("date", "t", "positives", "swabs", "TPR", "WeekEnd")])

ggplot(dati, aes(x = date, y = TPR, size = swabs)) +
  geom_point(alpha = 0.5) +
  scale_size(range = c(.01, 3), name="Swabs:  ",
             labels =  scales::label_number()) +
  labs(x = "", y = "Test positive rate") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels = monthyear_labels) +
  scale_y_continuous(lim = c(0,0.35), 
                     expand = c(0,0),
                     labels = scales::label_percent(accuracy = 1), 
                     sec.axis = dup_axis(name = ""),
                     minor_breaks = function(y) minor_breaks(y,2)) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", vjust = 1, size = 10),
        legend.pos = "bottom")

# fit model
mod = covindex_gam_betareg(
  y = dati$TPR, 
  t = dati$t, 
  x = data.frame(WeekEnd = dati$WeekEnd),
  wts = dati$swabs/mean(dati$swabs, na.rm=TRUE))
summary(mod)

# plot smoothing parameter selection
dt = data.table(k = mod$k, AIC = mod$AIC)
ggplot(dt, aes(x = k, y = AIC)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = mod$k_opt, lty = 2) +
  theme_minimal()

# check autocorrelation of residuals
acf(residuals(mod, type = "deviance"), lag.max = 21)
pacf(residuals(mod, type = "deviance"), lag.max = 21, ylim = c(-0.5, 0.5))

# predictions
pred = simpred_covindex_gam_betareg(
  mod,
  newdata = data.frame(t = dati$t, 
                       WeekEnd = FALSE) )

df = dati[, c("date", "TPR")]
CI_y = cbind(df, 
             list2DF(pred[c("mu", 
                            "lower_confint", "upper_confint", 
                            "lower_predint", "upper_predint")]))

ggplot(CI_y, aes(x = date, y = mu)) +
  # model fit
  geom_ribbon(aes(ymin = lower_confint, ymax = upper_confint),
              alpha = 1, fill = "grey50") +
  geom_ribbon(aes(ymin = lower_predint, ymax = upper_predint),
              alpha = 0.5, fill = "grey50") +
  geom_line(col = "dodgerblue2", lwd = 1) +
  # data
  geom_point(aes(y = TPR), size = 0.5) +
  #
  labs(y= "Test positive rate", x = "") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels = monthyear_labels) +
  scale_y_continuous(lim = c(0,0.34), 
                     breaks = seq(0,0.3,by=0.05),
                     expand = expansion(mult = 0, add = 0),
                     labels = scales::label_percent(accuracy = 1), 
                     sec.axis = dup_axis(name = ""),
                     minor_breaks = function(y) minor_breaks(y,2)) +
  theme_minimal()

# COVINDEX
covindex_sim = apply(pred$mu_sim, 2, covindex)
confint = apply(covindex_sim, 1, quantile, 
                prob = c(0.025, 0.975), na.rm = TRUE)
CI_covindex = cbind(df, 
                    covindex = covindex(pred$mu),
                    lower_confint = confint[1,], 
                    upper_confint = confint[2,])

ggplot(CI_covindex, aes(x = date, y = covindex)) +
  geom_hline(yintercept = 1, lty = 2) +
  # model fit
  geom_ribbon(aes(ymin = lower_confint, ymax = upper_confint),
              alpha = 1, fill = "grey") +
  geom_line(col = "dodgerblue2", lwd = 1) +
  #
  labs(y = "COVINDEX", x = "") + 
  scale_x_date(date_breaks = "1 month", 
               expand = c(0.01,0.01),
               date_labels =  monthyear_labels) +
  scale_y_continuous(limits = c(0.5,2), expand=c(0,0),
                     trans = "log10", 
                     breaks = seq(0.5, 2, by=0.1),
                     sec.axis = dup_axis(name=""),
                     minor_breaks = function(y) minor_breaks(y,1)) +
  theme_minimal()

