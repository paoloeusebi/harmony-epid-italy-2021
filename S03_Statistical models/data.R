load_install_package("data.table")

## Data from Dipartimento della Protezione Civile
## COVID-19 Italia - Monitoraggio della situazione
## http://arcg.is/C1unv
## Source: https://github.com/pcm-dpc/COVID-19

# to download from DPC
# dati_ita = fread("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv")
# or read data downloaded previously
dati_ita = fread("dpc-covid19-ita-andamento-nazionale.csv")

## Prepare the data for the analyses ----

dati_ita$date = as.Date(dati_ita$data)
# consider only molecular swabs and corresponding positive cases
dati_ita[date >= "2021-01-15", "tamponi"] =
  dati_ita[date >= "2021-01-15", "tamponi_test_molecolare"]
dati_ita[date >= "2021-01-15", "totale_casi"] =
  dati_ita[date >= "2021-01-15", "totale_positivi_test_molecolare"]
# check for consistency
dati_ita[, positives := { x = c(NA, diff(totale_casi)); ifelse(x < 0, NA, x) } ]
dati_ita[, swabs  := { x = c(NA, diff(tamponi)); ifelse(x < 0, NA, x) } ]
# correct for some errors in the recorded data
dati_ita[date == "2020-12-17", "swabs"] <- 185320
dati_ita[date == "2020-06-19", "positives"] <- 251
dati_ita[date == "2021-03-22", "positives"] <- 13846
dati_ita[date == "2021-03-22", "positives"] <- 13846
# set italian holidays for years 2020 and 2021
Italian_holidays = as.Date(
  c("2020-01-01", "2020-01-06", "2020-04-12", "2020-04-13", "2020-04-25",
    "2020-05-01", "2020-06-02", "2020-08-15", "2020-11-01", "2020-12-08",
    "2020-12-25", "2020-12-26",
    "2021-01-01", "2021-01-06", "2021-04-04", "2021-04-05", "2021-04-25",
    "2021-05-01", "2021-06-02", "2021-08-15", "2021-11-01", "2021-12-08",
    "2021-12-25", "2021-12-26"))
