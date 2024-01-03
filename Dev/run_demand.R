#Run demand_forecast

library(DemandForecast)

#prep = prep_demand_temp_data(include_na_volume = FALSE)
prep= prep_demand_temp_data(include_na_volume = FALSE, path = "/Users/Eirik/Desktop/Master2023/Data/")
X_mat = prep$X_mat                                      #Temperature field as columns
date_demand = prep$date_demand                          #Demand data with time covariates


a = Sys.time()
Test1 = demand_forecast(X_mat, date_demand, forc_start = '2016-01-01', forc_end = '2023-01-01',
                pred_win = 30, pred_lag = 15, train_y = 5, p_comps =4, no_pc = TRUE,
                custom = FALSE ,
                cores = 22, reg_form = "volume~ hour", Setup = "Rolling")
print( Sys.time()-a)

a = Sys.time()
Test1 = demand_forecast(X_mat, date_demand, forc_start = '2016-01-01', forc_end = '2023-01-01',
                        pred_win = 30, pred_lag = 15, train_y = 5, p_comps = 1, reg_form = "volume~ hour", no_pc = TRUE,
                        custom = FALSE ,
                        cores = 24,  Setup = "Rolling")
print( Sys.time()-a)


a = Sys.time()
Good_model3 = demand_forecast(X_mat, date_demand, forc_start = '2014-01-01', forc_end = '2023-01-01',
                        pred_win = 30, pred_lag = 15, train_y = 5, p_comps = 2, no_pc = TRUE,
                        custom = FALSE ,incl_climatology =TRUE, comb = FALSE, int = 14,
                        cores = 20, reg_form = "volume~ s(PC1) +s(PC2) + as.factor(hour) + as.factor(week) + as.factor(w_day) + as.factor(month) + year + as.factor(season) + w_day:month",
                        Setup = "Rolling_test")
print( Sys.time()-a)



forc_start ='2014-01-01'
forc_end = '2023-01-01'
pred_win = 30
pred_lag= 15
train_y=5
#reg_form = "volume ~  s(PC1)"
p_comps = 1
other_mods= NULL
comb = FALSE
custom = FALSE
incl_climatology = FALSE
no_pc = TRUE
gam_lasso = FALSE
cores = 29
Setup = "Rolling_test"
