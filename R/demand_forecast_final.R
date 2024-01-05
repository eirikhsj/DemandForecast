
#' Rolling CV forecast model for energy demand.
#'
#' @param forc_start String. Date for first initalization date.
#' @param forc_end String. Date for last initialization date.
#' @param pred_win Integer. Span of days in prediction window, default is 30.
#' @param pred_lag Integer. Days after initialization day the prediction window begins, default is 15.
#' @param train_y Integer. Years of training data.
#' @param reg_form String. Regression formula before PC covariates have been added.
#' @param p_comps Integer. Number of pc-comps in model.
#' @param other_mods List of other models to run or default is NULL.
#' @param comb Boolean. FALSE allows specific combination of PCs to be run without running all combinations.
#' @param custom Boolean. Specifies a custom pc based model.
#' @param incl_climatology Boolean. Includes climatology model.
#'
#' @return A data.table with actual and predicted demand values.
#'
#' @examples  Mod1 = demand_forecast_final(forc_start ='2016-01-01',forc_end = '2021-03-01',pred_win = 30, train_y = 3, pred_lag = 0,
#' reg_form = "volume ~ as.factor(hour) + as.factor(month) + year", p_comps = 3, other_mods= NULL,comb = TRUE)
#' @export
demand_forecast_final = function(X_mat, date_demand, NWP_pred, forc_start, forc_end, pred_win = 31,
                                 pred_lag= 0, train_y=3,reg_form, p_comps, other_mods= NULL,
                                 comb = TRUE, custom = FALSE,incl_climatology = FALSE, no_pc = TRUE, cores = 20, setup = "Rolling_final_ensamble"){

    last_poss_pred = range(date_demand$date)[2] - pred_win - pred_lag

    #1) Fix dates
    start = as.Date(forc_start)
    end = as.Date(forc_end)
    all_days = seq(start, end,  by = '1 days')
    init_days = all_days[mday(all_days)==1]

    ## **** Run parallel cores ****
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1) #Set number of threads

    detailed_results = parallel::mclapply(seq_along(init_days_all),
                                setup,
                                X_mat = X_mat, date_demand = date_demand,init_days= init_days, pred_win = pred_win, pred_lag= pred_lag, train_y=train_y,
                                reg_form= reg_form, p_comps= p_comps, other_mods= other_mods, comb = comb, custom = custom,
                                incl_climatology = incl_climatology, no_pc = no_pc,gam_lasso =gam_lasso, NWP_pred = NWP_pred,
                                mc.cores = cores)

    ## **** Return results ****
    Results = rbindlist(detailed_results)
    out = list()
    out$Results = Results
    print('Demand forecast is completed')
    return(out)
}

#' @export
Rolling_final = function(i, X_mat, date_demand, init_days,pred_win, pred_lag, train_y,
                   reg_form, p_comps, other_mods, comb, custom,incl_climatology, no_pc,num_thread,gam_lasso,NWP_pred = NWP_pred){

    ## ***** Step 1: Form the training and test datasets ****
    init_day = init_days[i]
    days_in_m = days_in_month(init_day)
    target_days = seq(init_day+ pred_lag, length.out = as.integer(days_in_m),  by = '1 days')
    print(paste('Forecast made on:', init_day))

    train_cutoff = seq(init_day,  length.out = 2, by = paste0('-',train_y, ' year'))[2]

    I_train = date_demand[(date > train_cutoff) & (date < init_day), I]
    I_test = date_demand[date %in% target_days, I]

    if (p_comps > 0){
        pc_data  = get_pca(X_mat, I_train, I_test, p_comps) # **** SVD ****
        dt_train = pc_data$dt_train
        dt_test  = pc_data$dt_test
    } else{
        dt_train = date_demand[I_train, ]
        dt_test  = date_demand[I_test, ]
    }
    print(paste('We are training on:', dim(dt_train)[1], 'observations'))
    max_year_train = dt_train[, max(year)]
    dt_test = dt_test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error
    n = nrow(dt_test)

    ## ***** Step 2: Train models ****

    if (p_comps > 0){                              #PC GAM MODELS
            mod = mgcv::gam(as.formula(sprintf("%s + %s", reg_form, custom)), data = dt_train)
    }

    ## **** Step 3: Check fit ****
    dt_test = merge(dt_test, NWP_pred, by = c("hour", "date"))
    results = dt_test
    results[,'lead_time':= (((as.integer(difftime(date, init_date, units = 'days')))*24)+ hour)/6]
    dt_PC = results[,.(date, hour, lead_time, PC1, PC2)]

    pred_demand_PC1_actual = predict(mod, newdata = dt_test)
    results[,'pred_demand_PC1_actual' := pred_demand_PC1_actual]
    dt_test[, PC1:= NULL]
    dt_test[, PC2:= NULL]
    pred_names1 = c("mean","q10","q20","q30","q40",
                    "q50","q60","q70","q80","q90")
    setnames(dt_test, old = "NWP_PC2_mean", new = "PC2")
    for (name in pred_names1){ #Fold all quantiles into the prediction
        PC1_name = paste0("NWP_PC1_",name)
        #PC2_name = paste0("NWP_PC2_",name)
        setnames(dt_test, old = PC1_name, new = "PC1")
        #setnames(dt_test, old = PC2_name, new = "PC2")
        pred_demand_test = predict(mod, newdata = dt_test)
        results[,paste0('pred_',name) := pred_demand_test]
        dt_test[, PC1:= NULL]
        #dt_test[, PC2:= NULL]
    }
    for (q in seq(10,90,10)){ #apply model error on each prediction
        pred_name= paste0("pred_q",q)
        col = data.table(t(results[, lapply(get(pred_name), function(x) quantile(x + rnorm(100, 0, mod$sig2^(0.5)),q/100))]))
        results[,paste0("final_pred_q_",q) := col]
    }

    final_results = merge(results, dt_PC, by = c("date", "hour", "lead_time"))

    if (is.data.table(final_results)==FALSE){
        print(str(final_results))
        print('Trying to convert')
        final_results = as.data.table(final_results)
        if(dim(final_results)[1] < 100){
            print('This is strange')
            print(final_results)
        }
    }
    return(final_results)
}



#' @export
Rolling_final_ensamble = function(i, X_mat, date_demand, init_days,pred_win, pred_lag, train_y,
                         reg_form, p_comps, other_mods, comb, custom,incl_climatology, no_pc,num_thread,gam_lasso,NWP_pred = NWP_pred){

    ## ***** Step 1: Form the training and test datasets ****
    init_day = init_days[i]
    days_in_m = lubridate::days_in_month(init_day)
    target_days = seq(init_day+ pred_lag, length.out = as.integer(days_in_m),  by = '1 days')
    print(paste('Forecast made on:', init_day))

    train_cutoff = seq(init_day,  length.out = 2, by = paste0('-',train_y, ' year'))[2]

    I_train = date_demand[(date > train_cutoff) & (date < init_day), I]
    I_test = date_demand[date %in% target_days, I]

    if (p_comps > 0){
        pc_data  = get_pca(X_mat = X_mat, date_demand = date_demand, I_train = I_train, I_test = I_test, p_comps = p_comps) # **** SVD ****
        dt_train = pc_data$dt_train
        dt_test  = pc_data$dt_test
    } else{
        dt_train = date_demand[I_train, ]
        dt_test  = date_demand[I_test, ]
    }
    print(paste('We are training on:', dim(dt_train)[1], 'observations'))
    max_year_train = dt_train[, max(year)]
    dt_test = dt_test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error
    n = nrow(dt_test)

    ## ***** Step 2: Train models ****
    if (p_comps > 0){                        #PC GAM MODELS
        t1 = Sys.time()
        mod = mgcv::gam(as.formula(sprintf("%s + %s", reg_form, custom)), data = dt_train)
        print(paste0('Time = ', Sys.time()- t1))
    }

    ## **** Step 3: Select Temperature Input ****


    ## **** Step 4: NWP ****
    y = year(init_day)
    m = month(init_day)
    date_fetch = as.Date(paste0(y, '-', m,'-01'))
    path = "~/bigdisk3/pro/sfe_daily_nordic_temperature/"
    pattern = 'sfe_nordic_temperature_'
    dt_file = paste0(y,'_',m)

    files_i = paste0(path, pattern, dt_file, '.nc4')
    nc = nc_open(files_i)
    forec = ncvar_get(nc, attributes(nc$var)$names[1])
    year_nwp = year(date_fetch) #year
    month_nwp = month(date_fetch) #month
    forc_name = paste0('forc_', year_nwp, '_', month_nwp)
    assign(forc_name, forec)
    print(c(year_nwp,month_nwp))
    nc_close(nc)

    #Get forecast and run pca- delete dt when done and store pca-results
    forecast = get(forc_name)
    print(dim(forecast))
    pc_data = get_pca(X_mat, date_demand= NULL, I_train, I_test, p_comps =p_comps,
                      NWP = forecast[,,,1:dim(forecast)[4]],
                      U = PC_ERA$U,
                      mu = PC_ERA$mu)
    rm(forecast)
    gc()

    if (any(is.na(pc_data$NWP_PC_mat[,pc_comp,]))==TRUE){
        pc_nwp = t(na.omit(t(pc_data$NWP_PC_mat[,pc_comp,])))
    } else{
        pc_nwp = pc_data$NWP_PC_mat[,pc_comp,]
    }
    ## **** Step 3: Check fit ****
    dt_test = dt_test[hour %in% c(6,12,18,24)]
    dt_results = cbind(dt_test, pc_nwp[1:dim(dt_test)[1],])

    #dt_results[,lapply(predict(mod, newdata = dt_test))]
    col_name = paste0('V',1:dim(pc_nwp)[2])
    col_out = paste0('pred_', col_name, sep = "")

    dt_results[, c(col_out) := lapply(.SD, multi_col_pred), .SDcols = col_name]
    dt_results[, pred_obs_temp := lapply(.SD, multi_col_pred), .SDcols = 'PC1']

    ## **** Step 4: Find Quantile ****
    dt_quant = data.table(t(dt_results[,.SD, .SDcols = col_out]))
    dt_quant_expand = replicate(1000, dt_quant, simplify = FALSE)
    dt_quant = rbindlist(dt_quant_expand)
    dt_quant_err = dt_quant[,.SD + matrix(rnorm(.N*dim(.SD)[2], 0, mod$sig2^(0.5)), nrow = .N, ncol = dim(.SD)[2])]

    dt_quant = data.table(t(dt_quant[,lapply(.SD ,  function(x) quantile(x, seq(0.1, 0.9, 0.1)))]))
    setnames(dt_quant, paste0('pred_q_',seq(0.1, 0.9, 0.1)))
    final_results = cbind(dt_results, dt_quant)
    return(final_results)
}


#' @export
multi_col_pred = function(col_name){
    res = data.table(hour = dt_results$hour,
                     w_day = dt_results$w_day,
                     week = dt_results$week,
                     month = dt_results$month,
                     season = dt_results$season,
                     year = dt_results$year,
                     PC1 = col_name)

    return(predict(mod, newdata = res))
}

