
#' Leave Year Out CV model
#'
#' @param choose_days_ahead Integer. Number of days ahead we are predicting.
#' @param mslp_lag Integer. How many days the weather data are lagged
#' @param plotting Boolean. If TRUE it plots Skill plot and Pred vs actual for Anomaly and ERA-FL1
#' @param form String. Regression formula.
#' @param dt_all  data.table with input data.
#'
#' @return data.table
#'
#' @examples mod_21_30_inter_nao = anomaly_lyocv(21, 30,form = "anomaly_future ~ FL_anomaly + FL_anomaly_30 + as.factor(month) * nao1_past")
#' @export
anomaly_lyocv = function(choose_days_ahead = 21,mslp_lag = 30, plotting =FALSE, dt_all = dt_all,
                         form = "anomaly_future ~ FL_anomaly + FL_anomaly_30 + as.factor(month) * nao1_past"){

    score_results = list()
    month_score_results = list()
    k = 0
    for (days_ahead in 1:choose_days_ahead){
        dt_all[, nao1_past := c(rep(NA, mslp_lag * 24 - 1), zoo::rollmean(FL1_mslp, mslp_lag * 24))]
        dt_all[, nao2_past := c(rep(NA, mslp_lag * 24 - 1), zoo::rollmean(FL2_mslp, mslp_lag * 24))]
        dt_all[, u_1hPa_past := c(rep(NA, mslp_lag * 24 - 1), zoo::rollmean(FL1_u1, mslp_lag * 24))]
        dt_all[, u_2hPa_past := c(rep(NA, mslp_lag * 24 - 1), zoo::rollmean(FL2_u1, mslp_lag * 24))]
        dt_all[, u_1hPa_past10 := c(rep(NA, mslp_lag * 24 - 1), zoo::rollmean(FL2_u10, mslp_lag * 24))]
        #dt_all[, F_temp_past := c(rep(NA, mslp_lag * 24 - 1), zoo::rollmean(FL1_temp, mslp_lag * 24))]
        #dt_all[, F_temp_past2 := c(rep(NA, mslp_lag * 24 - 1), zoo::rollmean(FL2_temp, mslp_lag * 24))]

        years_all = dt_all[, unique(year(date))]
        dt_test_all = list()

        # LOOCV - leaving out 1 year at a time
        for(y in seq_along(years_all)){
            k = k+1
            year_y = years_all[y]
            dt_clim = dt_all[year != year_y,
                             .("FL_temp_climatology" = mean(FL1_temp)),
                             .(month, day, hour)]
            dt_y = merge(dt_all, dt_clim, by = c("month", "day","hour"))
            setkey(dt_y, "date")
            dt_y[, FL_anomaly := FL1_temp - FL_temp_climatology]
            dt_y[, FL_anomaly_30 := c(rep(NA, mslp_lag * 24 - 1), zoo::rollmean(FL_anomaly, mslp_lag *24))]
            dt_y[, anomaly_future := shift(FL_anomaly,
                                           n = days_ahead * 24,
                                           NA,
                                           "lead")]
            dt_y[, FL_clima_future := shift(FL_temp_climatology,
                                            n = days_ahead * 24,
                                            NA,
                                            "lead")]
            dt_y[, FL1_future := shift(FL1_temp,
                                       n = days_ahead * 24,
                                       NA,
                                       "lead")]
            dt_train = dt_y[ year != year_y]
            dt_test = dt_y[year == year_y]
            ##    f = "anomaly_future ~ FL_anomaly + FL_anomaly_30 + as.factor(month) *  nao1_past + as.factor(month) * nao2_past"
            f = form
            #f = "anomaly_future ~ FL_anomaly + FL_anomaly_30 + as.factor(month) * nao2_past + as.factor(month)* u_1hPa_past"
            mod = lm(as.formula(f), data = dt_train)
            dt_test[, y_hat := predict(mod, newdata = dt_test)]
            dt_test[, days_ahead := days_ahead]
            dt_test = na.omit(dt_test)

            score_results[[k]] = data.table(year = year_y,
                                            days_ahead = days_ahead,
                                            score_clim = dt_test[, sqrt(mean(anomaly_future^2))],
                                            score_mod = dt_test[, sqrt(mean((anomaly_future - y_hat)^2))])
            month_score_results[[k]] = data.table(year = year_y,
                                                  days_ahead = days_ahead,
                                                  dt_test[, .(month_score_clim = sqrt(mean(anomaly_future^2)),
                                                              month_score_mod = sqrt(mean((anomaly_future - y_hat)^2))), by = month])
            dt_test_all[[y]] = dt_test
            print(score_results[[k]])
        }
        dt_results = rbindlist(dt_test_all)
        dt_results[, error := anomaly_future - y_hat]
        dt_results = dt_results[!is.na(error)]

        if (plotting == TRUE){
            par(mfrow = c(2,2))
            #Prediction plot Anomaly scale
            plot(dt_results[year %in% 2000& hour == 6, anomaly_future], type = 'l', ylim = c(-170, 170),
                 main = "Predicted vs Actual of Anomaly_FL1  (2000, hour 6)", ylab = 'Anomaly FL1', xlab = 'Day in year 2000')
            lines(dt_results[year %in% 2000& hour == 6, y_hat], col = 'red')
            abline(h = 0, lty = 2)
            legend("topright", legend = c("Anomaly", "Predicted Anomaly"), lty = 1, lwd = 2,
                   col = c('black', 'red'), bty = "n")

            #Prediction plot FL1 scale
            plot(dt_results[year %in% 2000& hour == 6, FL1_future], type = 'l', ylim = c(-200, 330),
                 main = "Predicted vs Actual of FL1 (2000, hour 6)",ylab = 'FL1', xlab = 'Day in year 2000')
            lines(dt_results[year %in% 2000& hour == 6, FL_clima_future], col = 'green')
            lines(dt_results[year %in% 2000& hour == 6, FL_clima_future + y_hat], col = 'red')
            legend("topright", legend = c("FL1", "Climatology","Climatology+Predicted Anomaly"), lty = 1, lwd = 2,
                   col = c('black', 'green','red'), bty = "n")

            #Skillplot by month Anomaly scale
            plot(dt_results[, 1 - sqrt(mean((anomaly_future - y_hat)^2))/sqrt(mean((anomaly_future)^2)), by = month], ylim = c(0, 0.57), type = 'l',
                 main = paste0('Skill score for ', days_ahead, ' days ahead, by month'), ylab = 'Skill score', xlab ='Month')

        }

    }
    final_score = rbindlist(score_results)
    final_month_score = rbindlist(month_score_results)
    out = c()
    out$final_score = final_score
    out$final_month_score = final_month_score
    return(out)
}
