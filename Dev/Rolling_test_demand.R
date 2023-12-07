Rolling_test = function(i,X_mat, date_demand, init_days,pred_win, pred_lag, train_y,
                   reg_form, p_comps, other_mods, comb, custom,incl_climatology, no_pc,num_thread,gam_lasso){
    init_day = init_days[i]
    target_days = seq(init_day+ pred_lag, length.out = pred_win,  by = '1 days')
    print(paste('Forecast made on:', init_day))
    train_cutoff = seq(init_day,  length.out = 2, by = paste0('-',train_y, ' year'))[2]
    I_train = date_demand[(date > train_cutoff) & (date < init_day), I]
    I_test = date_demand[date %in% target_days, I]
    print(length(I_train))
    if (p_comps > 0){
        pc_data  = get_pca(X_mat, I_train, I_test, p_comps) # **** SVD ****
        rm(X_mat)
        gc()
        dt_train = pc_data$dt_train
        dt_test  = pc_data$dt_test
    } else{
        dt_train = date_demand[I_train, ]
        dt_test  = date_demand[I_test, ]
    }
    dt_train= dt_train[!is.na(volume)]
    print(paste('We are training on:', dim(dt_train)[1], 'observations ', 'on iteration ', i))
    max_year_train = dt_train[, max(year)]
    dt_test = dt_test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error
    n = nrow(dt_test)
    results = dt_test
    results[,'init_date' := rep(init_day, length(dt_test$date))]
    results[,'lead_time':= ((as.integer(difftime(date, init_day, units = 'days')))*24)+ hour]
    results[,paste0('N_train') := rep(dim(dt_train)[1], dim(dt_test)[1])]
    if(incl_climatology == TRUE){

        results = Run_climatology(results = results, dt_train = dt_train, dt_test = dt_test, reg_form = reg_form,
                                  intervals = 14, basic = TRUE)

        }
    return(results)
}


Run_climatology = function(results, dt_train, dt_test, reg_form,
                           intervals = 8, basic = TRUE, daily_cycle = FALSE, weekday_cycle = FALSE){
    for (j in 1:intervals){

        ## **** Step 1:Fix leap year issue ****
        if(results[,max(.SD) , .SDcols = c(paste0('yr_by_',j))] >
           dt_train[,max(.SD) , .SDcols = c(paste0('yr_by_',j))]){
            print('Problem period')
            results[get(paste0('yr_by_',j))==max(get(paste0('yr_by_',j))), paste0('yr_by_',j) := max(get(paste0('yr_by_',j)))-1]
        }

        ## **** Step 2:Find Daily cycle  ****
        dt_train[, ave_vol := mean(volume), by = c(paste0('yr_by_',j), "hour")]
        climatology = dt_train[, .(ave_vol = mean(volume)), by = c(paste0('yr_by_',j), "hour")]
        results = merge(climatology, results, by= c(paste0('yr_by_',j), 'hour') )
        setnames(dt_train, "ave_vol", paste0('clim_pred_',j))
        setnames(results, "ave_vol", paste0('clim_pred_',j))

        ## **** Step 5:Run Models  ****
        mod_clim1 = mgcv::gam(as.formula(paste0(reg_form," + clim_pred_",j) ), data = dt_train)


        ## **** Step 6:Predict and log results  ****
        pred_demand_test_clim1 = predict(mod_clim1, newdata = results)
        results[,paste0('pred_','mod_clim1_',j) := pred_demand_test_clim1]

        ## **** Step 3:Find Weekday cycle  ****
        if (j>=7){
            dt_train[, ave_vol_wd := mean(volume), by = c(paste0('yr_by_',j), "w_day")]
            climatology2 = dt_train[, .(ave_vol_wd = mean(volume)), by = c(paste0('yr_by_',j), "w_day")]
            results = merge(climatology2, results, by= c(paste0('yr_by_',j), 'w_day') )
            setnames(dt_train, "ave_vol_wd", paste0('ave_vol_wd_',j))
            setnames(results, "ave_vol_wd", paste0('ave_vol_wd_',j))
            mod_clim2 = mgcv::gam(as.formula(paste0(reg_form, "+ ave_vol_wd_", j ) ) , data = dt_train)
            mod_clim3 = mgcv::gam(as.formula(paste0(reg_form," + clim_pred_",j," + ave_vol_wd_", j ) ) , data = dt_train)
            pred_demand_test_clim2 = predict(mod_clim2, newdata = results)
            pred_demand_test_clim3 = predict(mod_clim3, newdata = results)
            results[,paste0('pred_','mod_clim2_',j) := pred_demand_test_clim2]
            results[,paste0('pred_','mod_clim3_',j) := pred_demand_test_clim3]
        }

        if(j ==3){
            #print(paste0("PC1 + CLIM: ",  round(results[,sqrt(sum( (volume- get(paste0('pred_','mod_clim_',j)) )^2)/.N)],2)))
           # print(paste0("Climatology: ", round(results[,sqrt(sum( (volume- get(paste0('clim_pred_',j       )) )^2)/.N)],2)))
        }
    }
    return(results)

}

## **** Step 5:WRONG  ****
dt_test_new = dt_all[I %in% I_test,  ]
cols = c("date", "I", "hour", "w_day", "week", "month", "year", "PC1", "volume", paste0('yr_by_',j))
dt_all = rbind(dt_test[,.SD, .SDcols = cols], dt_train[,.SD, .SDcols= cols])
dt_all[, ave_vol := mean(volume), by = c(paste0('yr_by_',j), "hour")]
dt_all[, ave_vol_wd := mean(volume), by = c(paste0('yr_by_',j), "w_day")]
setnames(dt_all, "ave_vol", paste0('clim_pred_',j))
setnames(dt_all, "ave_vol_wd", paste0('ave_vol_wd_',j))
