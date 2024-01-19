# PanelQGARCH

## Installation

```R
#install.packages("devtools")
library(devtools)
install_github("wyLI2020/PanelQGARCH")
```

## Usage

```R
fr_onerep_tau(tau, Y_N_T, if_group=c(TRUE, TRUE, TRUE), if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123_candid_lower=c(0,0,0), delta_123_candid_upper=c(2,2,2), testLength=NA)
fr_rolling_forecast(t0, tau, Y_N_T, if_group=c(TRUE, TRUE, TRUE), if_mean=FALSE, criterion="BIC3", tolerance=0.0001, num_candid_onesplit=6, delta_123_candid_lower=c(0,0,0), delta_123_candid_upper=c(2,2,2), testLength=NA, clsNum)
```

- **t0**: integer, moving window in rolling forecasting procedure
- **tau**: decimal in (0,1), quantile level
- **Y_N_T**: (N, T) matrix, panel data
- **clsNum**: integer, number of core in parallel

## Example

```R
library(PanelQGARCH)
data("Y_N100_T1761")

N <- nrow(Y_N_T); T <- ncol(Y_N_T)

# Employ the proposed three-stage estimation procedure to fit the entire data set.
est_tau5per_vec <- fr_onerep_tau(tau=0.05, Y_N_T)
est_tau10per_vec <- fr_onerep_tau(tau=0.1, Y_N_T)
est_tau25per_vec <- fr_onerep_tau(tau=0.25, Y_N_T)
est_tau75per_vec <- fr_onerep_tau(tau=0.75, Y_N_T)
est_tau90per_vec <- fr_onerep_tau(tau=0.9, Y_N_T)
est_tau95per_vec <- fr_onerep_tau(tau=0.95, Y_N_T)

# Compare the forecasting performance of the fitted panel quantile GARCH models before and after homogeneity pursuit by evaluating the VaRs of the 100 stock returns
forecast_tau5per <- fr_rolling_forecast(t0=1000, tau=0.05, Y_N_T, clsNum=10)
forecast_tau10per <- fr_rolling_forecast(t0=1000, tau=0.1, Y_N_T, clsNum=10)
forecast_tau25per <- fr_rolling_forecast(t0=1000, tau=0.25, Y_N_T, clsNum=10)
forecast_tau75per <- fr_rolling_forecast(t0=1000, tau=0.75, Y_N_T, clsNum=10)
forecast_tau90per <- fr_rolling_forecast(t0=1000, tau=0.9, Y_N_T, clsNum=10)
forecast_tau95per <- fr_rolling_forecast(t0=1000, tau=0.95, Y_N_T, clsNum=10)
fr_rolling_forecast_summary_average(forecast_tau5per, sig_level=0.1)
fr_rolling_forecast_summary_average(forecast_tau10per, sig_level=0.1)
fr_rolling_forecast_summary_average(forecast_tau25per, sig_level=0.1)
fr_rolling_forecast_summary_average(forecast_tau75per, sig_level=0.1)
fr_rolling_forecast_summary_average(forecast_tau90per, sig_level=0.1)
fr_rolling_forecast_summary_average(forecast_tau95per, sig_level=0.1)
rolling_forecast_0.1_mat <- rbind(cbind(fr_rolling_forecast_summary_average(forecast_tau5per, sig_level=0.1),
            fr_rolling_forecast_summary_average(forecast_tau10per, sig_level=0.1),
            fr_rolling_forecast_summary_average(forecast_tau25per, sig_level=0.1)),
      cbind(fr_rolling_forecast_summary_average(forecast_tau75per, sig_level=0.1),
            fr_rolling_forecast_summary_average(forecast_tau90per, sig_level=0.1),
            fr_rolling_forecast_summary_average(forecast_tau95per, sig_level=0.1)))
round(rolling_forecast_0.1_mat[,c(1,2,5,6,7,10,11,12,15)], digits = 2)
```

