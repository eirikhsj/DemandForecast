#Run demand_forecast

library(DemandForecast)

prep = prep_demand_temp_data(include_na_volume = FALSE)
prep= prep_demand_temp_data(include_na_volume = TRUE,path = "/Users/Eirik/Desktop/Master2023/Data/")
X_mat = prep$X_mat                                      #Temperature field as columns
date_demand = prep$date_demand                          #Demand data with time covariates


Test1 = demand_forecast(X_mat, date_demand,forc_start = '2016-01-01', forc_end = '2023-01-01',
                pred_win = 30, pred_lag = 15, train_y = 5, p_comps = 0, no_pc = TRUE,
                custom = paste0("as.factor(hour)  + year + s(hourly_mean_grid)") ,
                cores = 4, reg_form = "volume~1")
