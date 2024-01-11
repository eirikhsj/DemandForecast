
#nr/project/stat/ClimateFutures/RenewableEnergy/SPMM/Volume/Data

library(DemandForecast)

##### 1. Get demand volume and ERA observed temperature data #################
PATH = "/nr/samba/user/esjavik/DemandData_link/"

MODPATH = "/nr/samba/user/esjavik/DemandModels/"

prep = prep_demand_temp_data(include_na_volume = FALSE, path = PATH) #Include NA to get pre-trained
#
X_mat = prep$X_mat #Temperature field
date_demand = prep$date_demand


##### 2. Get PC of ERA observed temperature data #################
if (file.exists(paste0(MODPATH,'PC_ERA_79_92.Rda'))){
    print('File exists - Loading file')
    load(paste0(MODPATH,'PC_ERA_79_92.Rda'))
} else {
    print('File does NOT exists - Creating file')
    PC_ERA = get_pre_trained_PC_ERA(file=PATH, X_mat,date_demand= date_demand,
        start_train = "1979-01-01", stop_train="1992-12-31", start_test = "1993-01-01", stop_test="2023-01-31",
        NWP=NA,run_again =TRUE)
    save(PC_ERA, file = paste0(MODPATH,'PC_ERA_79_92.Rda'))
}

##### 3. Get NWP quantiles #################

# Getting quantiles from jan 93 to may 2023 (the most current available)
if (file.exists(paste0(MODPATH,'NWP_quantiles_PC1_0193_0523.Rda'))){
    print('File exists - Loading file')
    load(paste0(MODPATH,'NWP_quantiles_PC1_0193_0523.Rda'))

} else {
    print('File does NOT exists - Creating file')
    NWP_quantiles_PC1_0193_0523 = get_all_NWP_quantiles(path = "/nr/samba/user/esjavik//bigdisk3_test/pro/sfe_daily_nordic_temperature/",
                                                        start_month  = '01', start_year = 1993, forc_months = 365,
                                                        PC_ERA_path = "/nr/samba/user/esjavik/DemandData/PC_ERA_79_92.Rda",
                                                        pc_comp = 1,reweight = FALSE)
    save(NWP_quantiles_PC1_0193_0523, file = paste0(PATH,'NWP_quantiles_PC1_0193_0523.Rda'))

}

##### 3. Prepare basic NWP run #################

ERA_NWP = get_ERA_NWP(ERA_path = '/nr/samba/user/esjavik/DemandModels/PC_ERA_79_92.Rda',
                       NWP_path = '/nr/samba/user/esjavik/DemandModels/NWP_quantiles_PC1_0193_0523.Rda',
                       quant = '90',
                       PC = 'PC1',
                       NWP_preds = 2,
                       reweight = FALSE,
                       tuning_k = 10
)


ERA_NWP2 = get_ERA_NWP(ERA_path = '/nr/samba/user/esjavik/DemandModels/PC_ERA_79_92.Rda',
                      NWP_path = '/nr/samba/user/esjavik/DemandModels/NWP_quantiles_PC1_0193_0523.Rda',
                      quant = '90',
                      PC = 'PC2',
                      NWP_preds = 2,
                      reweight = FALSE,
                      tuning_k = 10
)
ERA_NWP[,PC2 := ERA_NWP2[,PC2]]


ERA_NWP[, season:= ifelse(month(date) %in% c(12,1,2), 1,
                          ifelse(month(date) %in% c(3,4,5), 2,
                                 ifelse(month(date) %in% c(6,7,8), 3, 4)))]


##### 5. Get re-weighted NWP quantiles #################
NWP_quant_rew25_1993_2023_1 = get_all_NWP_quantiles(path = "/nr/samba/user/esjavik//bigdisk3_test/pro/sfe_daily_nordic_temperature/",
                                                    pattern = 'sfe_nordic_temperature_',
                                                    start_month  = '01', start_year = 1993, forc_months = 365,
                                                    PC_ERA_path = "/nr/samba/user/esjavik/DemandModels/PC_ERA_79_92.Rda",
                                                    pc_comp = 1,
                                                    reweight = TRUE, rew_int =c(15, 1)) # 1 day, last day 15th

NWP_quant_rew25_1993_2023_2 = get_all_NWP_quantiles(path = "/nr/samba/user/esjavik//bigdisk3_test/pro/sfe_daily_nordic_temperature/",
                                                    pattern = 'sfe_nordic_temperature_',
                                                    start_month  = '01', start_year = 1993, forc_months = 365,
                                                    PC_ERA_path = "/nr/samba/user/esjavik/DemandModels/PC_ERA_79_92.Rda",
                                                    pc_comp = 1,
                                                    reweight = TRUE, rew_int =c(15, 2))  # 2 days

NWP_quant_rew25_1993_2023_3 = get_all_NWP_quantiles(path = "/nr/samba/user/esjavik//bigdisk3_test/pro/sfe_daily_nordic_temperature/",
                                                    pattern = 'sfe_nordic_temperature_',
                                                    start_month  = '01', start_year = 1993, forc_months = 365,
                                                    PC_ERA_path = "/nr/samba/user/esjavik/DemandModels/PC_ERA_79_92.Rda",
                                                    pc_comp = 1,
                                                    reweight = TRUE, rew_int =c(15, 3)) # 3 days


# ERA_path = '/nr/samba/user/esjavik/DemandData/PC_ERA_79_92.Rda'
# NWP_path = '/nr/samba/user/esjavik/DemandData/NWP_quantiles_0193_0523.Rda'
# quant = '90'
# PC = 'PC1'
# NWP_preds = 2
# reweight = FALSE
# tuning_k = 10

t(Good_model[, lapply(.SD, function(x) {sqrt(mean( (x - volume)^2 ))}), .SDcols = c(paste0("clim_pred_", 1:14), paste0("pred_mod_clim1_", 1:14), paste0("pred_mod_clim2_", 7:14),paste0("pred_mod_clim3_", 7:14))])

t(Good_model[, lapply(.SD, function(x) {sqrt(mean( (x - volume)^2 ))}), .SDcols = c(paste0("clim_pred_", 1:4), paste0("pred_mod_clim1_", 1:4))])


t(Test1$Results[, lapply(.SD, function(x) {sqrt(mean( (x - volume)^2 ))}), .SDcols = c(paste0("clim_pred_", 1:14), paste0("pred_mod_clim1_", 1:14), paste0("pred_mod_clim2_", 7:14),paste0("pred_mod_clim3_", 7:14))])

a = Sys.time()
New_inter3_smooth2 = demand_forecast(X_mat, date_demand, forc_start = '2016-01-01', forc_end = '2023-01-01',
                                          pred_win = 30, pred_lag = 15, train_y = 3, p_comps = 2, no_pc = FALSE, incl_climatology = FALSE,
                                          reg_form = "volume~ as.factor(hour) + as.factor(w_day) + as.factor(month) + as.factor(season) + as.factor(week) + year + as.factor(w_day):as.factor(month) +as.factor(hour):as.factor(week) + as.factor(hour):as.factor(w_day)" ,
                                          cores = 27, custom = "res ~ s(PC1, by = month) + s(PC2) ", Setup = "Rolling")
print( Sys.time()-a)
dev.new()
plot(Test1$Results[month == 9, .(PC1, (volume-pred_1))], cex = 0.1, pch = 16, col = 'blue')
lines(seq(from = -400, to = 400), prz_9, col = 'red')


#Summarize testing done

#Approach1 : Using rolling_test to check climatology interaction
Good_model3 = demand_forecast(X_mat, date_demand, forc_start = '2016-01-01', forc_end = '2023-01-01', pred_win = 30, pred_lag = 15, train_y = 5, p_comps = 2, no_pc = TRUE,
                              custom = FALSE ,incl_climatology =TRUE, comb = FALSE, int = 7, TRUE,
                              cores = 27, reg_form = "volume~ s(PC1) +s(PC2) + as.factor(hour) + as.factor(week) + as.factor(w_day) + as.factor(month) + year + as.factor(season) + as.factor(w_day):as.factor(month)",
                              Setup = "Rolling_test")
t(Good_model3$Results[, lapply(.SD, function(x) {sqrt(mean( (x - volume)^2 ))}), .SDcols = c(paste0("clim_pred_", 7), paste0("pred_mod_clim1_", 7), paste0("pred_mod_clim2_", 7),paste0("pred_mod_clim3_", 7))])


#Approach2: anomaly

a = Sys.time()
New_inter3_smooth2 = demand_forecast(X_mat, date_demand, forc_start = '2016-01-01', forc_end = '2023-01-01',
                                     pred_win = 30, pred_lag = 15, train_y = 3, p_comps = 2, no_pc = FALSE, incl_climatology = FALSE,
                                     reg_form = "volume~ as.factor(hour) + as.factor(w_day) + as.factor(month) + as.factor(season) + as.factor(week) + year + as.factor(w_day):as.factor(month) +as.factor(hour):as.factor(week) + as.factor(hour):as.factor(w_day)" ,
                                     cores = 26, custom = "res ~ s(PC1_anomal, by = season) + s(PC2) ", Setup = "Rolling_anomaly")
print( Sys.time()-a)

t(New_inter3_smooth2$Results[, lapply(.SD, function(x) {sqrt(mean( (x - volume)^2 ))}), .SDcols = c( paste0("pred_",1:2))])


#Approach 3: extra interaction

a = Sys.time()
New_inter3 = demand_forecast(X_mat, date_demand, forc_start = '2016-01-01', forc_end = '2023-01-01',
                             pred_win = 30, pred_lag = 0, train_y = 3, p_comps = 2, no_pc = TRUE, incl_climatology = FALSE,
                             reg_form = "volume~ as.factor(hour) + as.factor(w_day) + as.factor(month) + as.factor(season) + as.factor(week) + year + as.factor(w_day):as.factor(month) +as.factor(hour):as.factor(week) + as.factor(hour):as.factor(w_day) " ,
                             cores = 22, custom = FALSE, Setup = "Rolling")
print( Sys.time()-a)
t(New_inter3$Results[, lapply(.SD, function(x) {sqrt(mean( (x - volume)^2 ))}), .SDcols = c( paste0("pred_",1:3))])
New_inter3$Results[month==12, lapply(.SD, function(x) {sqrt(mean( (x - volume)^2 ))}), .SDcols = c( paste0("pred_",3)), by = .(hour, month)][, .SD,keyby = pred_3]


New_inter3$Results[,PC1_clim := mean(PC1), by =c('year_day', 'hour')]
New_inter3$Results[,PC1_anomal := PC1 - PC1_clim]
New_inter3$Results[,res := volume-PC1,]




plot(ERA_NWP[year(date)>2007,sd(raw_loss), keyby = lead_time], type = 'l', ylim = c(0, 22), col = 'purple')
legend('bottomright', col = c('black','red', 'purple'), title = 'Model', legend = c('NWP1 (Weighted /direct)', 'NWP1 (Qreg All data)', "NWP1 (After 2007)"), pch = 16, lty = 1)

sam = sample(1:80000,5000)
pl = data.table(X_mat[sam,c(1,200,460)])
fig <- plot_ly(pl, x = ~V1, y = ~V2, z = ~V3, marker = list(size = 2.3, opacity = 0.3))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Temp K'),
                                   yaxis = list(title = 'Temp K'),
                                   zaxis = list(title = 'Tempk K')))

fig





library(raster)
imported_raster=raster('gpw_v4_population_count_rev11_2020_1_deg.tif')


nc = nc_open("gpw_v4_population_count_rev11_1_deg.nc")
names(nc[['var']])
pop = ncvar_get(nc, names(nc[['var']]))
pop = data.table(pop[,,1])

log_pop = pop[,names(pop):= lapply(.SD, function(x) log(replace(x, is.na(x), 1))), .SDcols = names(pop) ]

