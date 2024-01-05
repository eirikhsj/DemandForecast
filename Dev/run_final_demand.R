# Run final script

install.packages('DemandForecast', repos = NULL,type = "source")

library(DemandForecast)

prep = prep_demand_temp_data(include_na_volume = FALSE)
#
X_mat = prep$X_mat #Temperature field
date_demand = prep$date_demand
#dt_date = prep$dt_date
#dt_unified = merge(dt_demand, dt_date, by = c("date", "hour"))

#Specify time covariates
date_demand[,week:= week(date)]
date_demand[,week1 := cos(2*pi * week/52) ]
date_demand[,week2 := sin(2*pi * week/52) ]

date_demand[,month:= month(date)]
date_demand[,month1 := cos(2*pi * month/12) ]
date_demand[,month2 := sin(2*pi * month/12) ]

date_demand[,hour1 := cos(2*pi * hour/24) ]
date_demand[,hour2 := sin(2*pi * hour/24) ]

date_demand[,year:= year(date)]

date_demand[,w_day:= as.numeric(format(date, "%u"))]
date_demand[,w_day1 := cos(2*pi * w_day/7) ]
date_demand[,w_day2 := sin(2*pi * w_day/7) ]

date_demand[, season:= ifelse(month(date) %in% c(12,1,2), 1,
                              ifelse(month(date) %in% c(3,4,5), 2,
                                     ifelse(month(date) %in% c(6,7,8), 3, 4)))]

date_demand[,season1 := round(cos(2*pi * season/4), digits = 7)]
date_demand[,season2 := round(sin(2*pi * season/4), digits = 7)]

date_demand[,hourly_mean_grid := rowMeans(X_mat)]


#Then incorporate the NWP data
# PC1_quants = get_all_NWP_quantiles(path = "/mn/kadingir/datascience_000000/eirikhsj/sfe_nordic_temperature/",
#                                    pattern = 'sfe_nordic_temperature_',
#                                    start_month  = '01', start_year = 2012, forc_months = 133,
#                                    PC_ERA_path = "/mn/kadingir/datascience_000000/eirikhsj/PC_ERA_79_92.Rda",
#                                    pc_comp = 1,
#                                    reweight = FALSE,
#                                    rew_int = c(15,4))
# PC2_quants = get_all_NWP_quantiles(path = "/mn/kadingir/datascience_000000/eirikhsj/sfe_nordic_temperature/",
#                                    pattern = 'sfe_nordic_temperature_',
#                                    start_month  = '01', start_year = 2012, forc_months = 133,
#                                    PC_ERA_path = "/mn/kadingir/datascience_000000/eirikhsj/PC_ERA_79_92.Rda",
#                                    pc_comp = 2,
#                                    reweight = FALSE,
#                                    rew_int = c(15,4))
#setnames( PC2_quants$NWP, old = "NWP_PC1_mean", new = "NWP_PC2_mean")
load("PC1_quants.Rda")
load("PC2_quants.Rda")


load("New_Clima.Rda")

NWP_pred = merge(PC1_quants$NWP[date>as.Date("2012-12-31"),],
                 PC2_quants$NWP[date>as.Date("2012-12-31"),],
                 by = c("hour", "date", "init_date"))


#MODEL RUNS

S1_S2 = demand_forecast_final(X_mat, date_demand, NWP_pred, forc_start = '2016-01-01', forc_end = '2023-01-01', pred_win = 31, pred_lag= 0, train_y=3,
                                 reg_form = "volume ~ 1", p_comps=2, other_mods= NULL, comb = FALSE, custom = "s(PC1) + s(PC2)",
                                 incl_climatology = FALSE, no_pc = TRUE, cores = 44)
