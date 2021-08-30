library(EpiEstim)
library(tidyverse)

# data munging ------------------------------------------------------------
wd <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

wd1 <- wd %>%
  as_tibble() %>%
  rename_with( ~ sub("X", "", .x), starts_with("X")) %>%
  gather(date, cases, 5:ncol(wd), factor_key = TRUE) %>%
  mutate(Country.Region = as.factor(Country.Region),
         date = as.Date(date, "%m.%d.%y"))

wd_long <- wd1 %>%
  mutate(date_num = as.numeric(date)) %>%
  group_by(Country.Region) %>% # maybe not needed
  arrange(date_num) %>%
  mutate(ncases = cases - lag(cases)) %>%
  ungroup()

ita_long <- wd_long %>%
  filter(Country.Region=="Italy") %>%
  select(date, cases, ncases) %>%
  mutate(ncases=ifelse(ncases<0, 0, ncases)) %>%
  filter(date>"2020-01-22")# to correct negative cases in 2020-06-19

## incidence data
head(ita_long)

# plot the daily incidence data
ggplot(data=ita_long, aes(x=date, y = ncases)) +
  geom_line()

# Calculate Rt
# In this example, we only specify the mean and standard deviation of the serial interval. In that case an offset gamma distribution is used for the serial interval. In the following example, we use the mean and standard deviation of the serial interval for flu from Ferguson et al., Nature, 2005:
res_parametric_si <- estimate_R(ita_long$ncases,
                                method = "parametric_si",
                                config = make_config(list(mean_si = 2.6,
                                                          std_si = 1.5)))
# check warnings!

# fit after a reasonable time
ita_long2 <- ita_long %>%
  filter(date>"2020-03-01")
# fit again
res_parametric_si <- estimate_R(ita_long2$ncases,
                                method = "parametric_si",
                                config = make_config(list(mean_si = 2.6,
                                                          std_si = 1.5)))

head(res_parametric_si$R)
plot(res_parametric_si, "R")
plot(res_parametric_si, legend = FALSE)
