#Run demand_forecast

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

# for (i in 3:10){
#
#     mod = demand_forecast(X_mat, date_demand,forc_start = '2016-01-01', forc_end = '2023-01-01',
#                                       pred_win = 30, pred_lag = 15, train_y = 5, p_comps = i,  no_pc = TRUE,
# custom = paste0("as.factor(hour):as.factor(week)  + as.factor(w_day):as.factor(week) + s(week) + s(month) + season1 + season2 + year + s(PC1)+ s(PC",i,")") ,
# cores = 48, reg_form = "volume~1")
#     assign(paste0("Custom_comb",i, "_int_w_w"), mod)
#
#     save(list = c(paste0("Custom_comb",i, "_int_w_w")), file = paste0("Custom_comb",i, "_int_w_w",'.Rda'))
# }

for (i in 3:10){

    mod = demand_forecast(X_mat, date_demand,forc_start = '2016-01-01', forc_end = '2023-01-01',
                          pred_win = 30, pred_lag = 15, train_y = 5, p_comps = i,  no_pc = TRUE,
                          custom = paste0("as.factor(hour):as.factor(month)  + as.factor(w_day):as.factor(week) + s(week) + s(month) + season1 + season2 + year + s(PC1)+ s(PC",i,")") ,
                          cores = 48, reg_form = "volume~1")
    assign(paste0("Custom_comb",i, "_int_m_w"), mod)

    save(list = c(paste0("Custom_comb",i, "_int_m_w")), file = paste0("Custom_comb",i, "_int_m_w",'.Rda'))
}

for (i in 3:10){

    mod = demand_forecast(X_mat, date_demand,forc_start = '2016-01-01', forc_end = '2023-01-01',
                          pred_win = 30, pred_lag = 15, train_y = 5, p_comps = i,  no_pc = TRUE,
                          custom = paste0("as.factor(hour):as.factor(week)  + as.factor(w_day):as.factor(month) + s(week) + s(month) + season1 + season2 + year + s(PC1)+ s(PC",i,")") ,
                          cores = 48, reg_form = "volume~1")
    assign(paste0("Custom_comb",i, "_int_w_m"), mod)

    save(list = c(paste0("Custom_comb",i, "_int_w_m")), file = paste0("Custom_comb",i, "_int_w_m",'.Rda'))
}

for (i in 3:10){

    mod = demand_forecast(X_mat, date_demand,forc_start = '2016-01-01', forc_end = '2023-01-01',
                          pred_win = 30, pred_lag = 15, train_y = 5, p_comps = i,  no_pc = TRUE,
                          custom = paste0("as.factor(hour):as.factor(month)  + as.factor(w_day):as.factor(month) + s(week) + s(month) + season1 + season2 + year + s(PC1)+ s(PC",i,")") ,
                          cores = 48, reg_form = "volume~1")
    assign(paste0("Custom_comb",i, "_int_m_m"), mod)
    save(list = c(paste0("Custom_comb",i, "_int_m_m")), file = paste0("Custom_comb",i, "_int_m_m",'.Rda'))
}
