
library(data.table)

#1) First load all  data
load("~/Desktop/Master2023/Data/PC_ERA/PC_ERA_79_23.Rda")
load("~/Desktop/Master2023/Data/u_component_of_wind_1hPa_pca_analysis.RData")
U1hPa_factor_loadings = dt_factor_loadings
U1hPa_factor_loadings[,hour := hour+1]
load("~/Desktop/Master2023/Data/u_component_of_wind_10hPa_pca_analysis.RData")
U10hPa_factor_loadings = dt_factor_loadings
U10hPa_factor_loadings[,hour := hour+1]
load("~/Desktop/Master2023/Data/mean_sea_level_pressure_nao_standardized_pca_analysis.RData")
NAO_factor_loadings = dt_factor_loadings
NAO_factor_loadings[,hour := hour+1]

#2) Combine into one data.table
dt_all = merge(PC_ERA_79_23$dt_train[, .(date, hour,
                                             "FL1_temp" = PC1,
                                             "FL2_temp" = PC2)],
               NAO_factor_loadings[, .(date,
                                             hour,
                                             "FL1_mslp" = V1,
                                             "FL2_mslp" = V2)],
                                    by = c("date", "hour"))
dt_all = merge(dt_all,
               U10hPa_factor_loadings[, .(date,
                                             hour,
                                             "FL1_u10" = V1,
                                             "FL2_u10" = V2)],
                                    by = c("date", "hour"))
# dt_all = merge(dt_all,
#                l_u_850$dt_factor_loadings[, .(date,
#                                               hour,
#                                               "FL1_u850" = V1,
#                                               "FL2_u850" = V2)],
#                                   by = c("date", "hour"))
dt_all = merge(dt_all,
               U1hPa_factor_loadings[, .(date,
                                            hour,
                                            "FL1_u1" = V1,
                                            "FL2_u1" = V2)],
                                    by = c("date", "hour"))

#dt_all = load_pca_framework(out_dir = '/Users/Eirik/Desktop/Master2023/Data') ## put out_dir in this with your path
dt_all[, year := year(date)]
dt_all[, month := month(date)]
dt_all[, day := mday(date)]


#3) Run models
mod_21_30_inter_nao = anomaly_lyocv(21, 30,form = "anomaly_future ~ FL_anomaly + FL_anomaly_30 + as.factor(month) * nao1_past")
mod_21_30_inter_nao_no_30anom = anomaly_lyocv(21, 30 ,form = "anomaly_future ~ FL_anomaly + as.factor(month) * nao1_past")
mod_21_30_inter_nao1_nao2 = anomaly_lyocv(21, 30 ,form = "anomaly_future ~ FL_anomaly + FL_anomaly_30 + as.factor(month) * nao1_past+ as.factor(month) * nao2_past")
mod_21_30_inter_u_1 = anomaly_lyocv(21, 30 ,form = "anomaly_future ~ FL_anomaly + FL_anomaly_30 + as.factor(month) * u_1hPa_past")
mod_21_30_inter_u_10 = anomaly_lyocv(21, 30 ,form = "anomaly_future ~ FL_anomaly + FL_anomaly_30 + as.factor(month) * u_1hPa_past10") #Not run
mod_21_30_just_anom = anomaly_lyocv(21, 30 ,form = "anomaly_future ~ FL_anomaly + FL_anomaly_30 + as.factor(month)")

#4) Save model output
save(mod_21_30_inter_nao, file = '/Users/Eirik/Desktop/Master2023/Output_data/output_loocv_nao_u_wind/mod_21_30_inter_nao.RData')
save(mod_21_30_inter_nao_no_30anom, file = '/Users/Eirik/Desktop/Master2023/Output_data/output_loocv_nao_u_wind/mod_21_30_inter_nao_no_30anom.RData')
save(mod_21_30_inter_nao1_nao2, file = '/Users/Eirik/Desktop/Master2023/Output_data/output_loocv_nao_u_wind/mod_21_30_inter_nao1_nao2.RData')
save(mod_21_30_inter_u_1, file = '/Users/Eirik/Desktop/Master2023/Output_data/output_loocv_nao_u_wind/mod_21_30_inter_u_1.RData')
save(mod_21_30_inter_u_10, file = '/Users/Eirik/Desktop/Master2023/Output_data/output_loocv_nao_u_wind/mod_21_30_inter_u_10.RData')
save(mod_21_30_just_anom, file = '/Users/Eirik/Desktop/Master2023/Output_data/output_loocv_nao_u_wind/mod_21_30_just_anom.RData')


#5) Plot performance

#Plotting Skill score
par(mfrow = c(1,1))
plot(mod_21_30_inter_nao$final_score[!is.na(score_mod),1 - mean(score_mod) / (mean(score_clim)), by = days_ahead], type = 'l', ylim = c(0, 0.15),
     main = 'Skill scores, different LYOCV models', xlab = 'Days ahead', ylab = 'Skill score')
lines(mod_21_30_inter_u_1$final_score[!is.na(score_mod),1 - mean(score_mod) / (mean(score_clim)), by = days_ahead], col ="red")
lines(mod_21_30_inter_u_10$final_score[!is.na(score_mod),1 - mean(score_mod) / (mean(score_clim)), by = days_ahead], col ="blue")

lines(mod_21_30_inter_nao1_nao2$final_score[!is.na(score_mod),1 - mean(score_mod) / (mean(score_clim)), by = days_ahead], col ="purple")
lines(mod_21_30_just_anom$final_score[!is.na(score_mod),1 - mean(score_mod) / (mean(score_clim)), by = days_ahead], col ="green")
lines(mod_21_30_inter_nao_no_30anom$final_score[!is.na(score_mod),1 - mean(score_mod) / (mean(score_clim)), by = days_ahead], col ="orange")
legend("topright", legend = c("NAO_1", "U_1",'U_10', "NAO_1 + NAO2", 'No Weather data', 'NAO_1 (no anom_30)'), lty = 1, lwd = 2,
       col = c('black', 'red', 'blue', 'purple', 'green', 'orange'), bty = "n")
abline(h = 0.1, lty = 2)

#mod_21_30_inter_nao[!is.na(score_mod),1 - mean(score_mod) / (mean(score_clim)), by = days_ahead]

#Plotting skill by season
par(mfrow = c(1,1))
plot(mod_21_30_inter_u_1$final_month_score[month%in%c(12,1,2),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead],
     type = 'l', ylim = c(0, 0.15), col = 'blue', main = 'Skill score by season (U_wind model)', xlab = 'Days ahead', ylab = 'Skill score')
lines(mod_21_30_inter_u_1$final_month_score[month%in%c(3,4,5),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead], col = 'green')
lines(mod_21_30_inter_u_1$final_month_score[month%in%c(6,7,8),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead], col = 'red')
lines(mod_21_30_inter_u_1$final_month_score[month%in%c(9,10,11),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead], col = 'brown')
legend("topright", legend = c("Winter", "Spring", "Summer", 'Autumn'), lty = 1, lwd = 2,
       col = c('blue', 'green', 'red', 'brown'), bty = "n")
abline(h = 0.1, lty = 2)

#Plotting skill for winter months (several models)
par(mfrow = c(1,1))
plot(mod_21_30_inter_nao$final_month_score[month%in%c(12,1,2),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead],
     type = 'l', ylim = c(0, 0.15), main = 'Skill score for winter months (monthly interaction models)', xlab = 'Days ahead', ylab = 'Skill score')
lines(mod_21_30_inter_u_1$final_month_score[month%in%c(12,1,2),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead], col = 'red')
lines(mod_21_30_inter_u_10$final_month_score[month%in%c(12,1,2),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead], col = 'blue')

lines(mod_21_30_inter_nao1_nao2$final_month_score[month%in%c(12,1,2),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead], col = 'purple')
lines(mod_21_30_just_anom$final_month_score[month%in%c(12,1,2),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead], col = 'green')
lines(mod_21_30_inter_nao_no_30anom$final_month_score[month%in%c(12,1,2),1 - mean(month_score_mod) / (mean(month_score_clim)), by = days_ahead], col = 'orange')
legend("topright", legend = c("NAO_1", "U_1",'U_10', "NAO_1 + NAO2", 'No Weather data', 'NAO_1 (no anom_30)'), lty = 1, lwd = 2,
       col = c('black', 'red', 'blue','purple', 'green', 'orange'), bty = "n")
abline(h = 0.1, lty = 2)


matplot(PC_ERA_79_23$dt_train[, .(PC1), by = yday(date)][,PC1], type = 'l')
lines(PC_ERA_79_23$dt_train[,.(date, "climatology" = mean(PC1)),.(month(date), mday(date), hour)][,climatology], col = 'red')




