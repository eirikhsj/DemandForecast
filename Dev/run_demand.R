#Run demand_forecast

library(DemandForecast)

library(data.table)
library(ncdf4)
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

forc_start = '2014-02-01'; forc_end = '2021-09-01'; pred_win = 30; pred_lag = 15; train_y = 5;
p_comps = 0;  reg_form = "volume ~ as.factor(hour) + as.factor(month) + year"


Best2 = demand_forecast(X_mat = X_mat, date_demand = date_demand, forc_start = '2015-06-01', forc_end = '2021-09-01', pred_win = 30, pred_lag = 15, train_y = 5,
                        p_comps = 6,  reg_form = "volume ~ as.factor(hour):as.factor(month) + as.factor(season) + as.factor(w_day):as.factor(season) + s(week1) + s(week2) + month + year")

save(Best2, file = 'Best2.RData')
