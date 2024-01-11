
#' Rolling CV forecast model for energy demand.
#'
#' @param forc_start Date for first initalization date.
#' @param forc_end Date for last initialization date.
#' @param pred_win Integer. Span of days in prediction window, default is 30.
#' @param pred_lag Integer. Days after initialization day the prediction window begins, default is 15.
#' @param train_y Integer. Years of training data.
#' @param reg_form String. Regression formula before PC covariates have been added.
#' @param p_comps Integer. Number of pc-comps in model.
#' @param other_mods List of other models to run or default is NULL.
#' @param comb Boolean. FALSE allows specific combination of PCs to be run without running all combinations.
#' @param custom Boolean. Specifies a custom pc based model.
#' @param incl_climatology Boolean. TRUE includes climatology model.
#' @param no_pc Boolean. TRUE includes basic (no pc) model.
#' @param gam_lasso Boolean. TRUE sets lasso penalty.
#' @param cores Integer. Sets number of parallelized cores, default is 4.
#'
#' @return A data.table with actual and predicted demand values.
#' @export
#' @examples  Mod1 = demand_forecast(forc_start ='2014-01-01',forc_end = '2023-01-01',pred_win = 30, train_y = 5,
#' pred_lag = 15,reg_form = "volume ~ as.factor(hour) + as.factor(month) + year", p_comps = 3, other_mods= NULL,comb = TRUE)
#'
#'
demand_forecast = function(X_mat, date_demand, forc_start, forc_end, pred_win = 30, pred_lag= 15, train_y=5,
                           reg_form, p_comps, other_mods= NULL, comb = FALSE, custom = FALSE,
                           incl_climatology = FALSE, no_pc = TRUE, gam_lasso = FALSE, cores = 4, int = 1,
                           Setup = "Rolling"){

    last_poss_pred = range(date_demand$date)[2] - pred_win - pred_lag

    # **** Fix dates ****
    start = as.Date(forc_start)
    end = as.Date(forc_end)
    all_days = seq(start, end,  by = '1 days')
    init_days = all_days[mday(all_days)==1]

    # **** Run parallel cores ****
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1) #Set number of threads

    detailed_results = parallel::mclapply(seq_along(init_days),
                      Setup,
                      X_mat = X_mat, date_demand = date_demand,init_days= init_days, pred_win = pred_win, pred_lag= pred_lag, train_y=train_y,
                      reg_form= reg_form, p_comps= p_comps, other_mods= other_mods, comb = comb, custom = custom,
                      incl_climatology = incl_climatology, no_pc = no_pc, gam_lasso =gam_lasso, int = int,
                      mc.cores = cores, mc.silent= FALSE, mc.preschedule = FALSE)

    # **** Return results ****
    Results = tryCatch({rbindlist(detailed_results, fill = TRUE)},
                       error = function(cond){print('Could not use rbind, there probably was an error in a core')
                           print(cond)
                           Results = detailed_results
                       })
    out = list()
    out$Results = Results
    print('Demand forecast is completed')
    return(out)

}

#' Basic Model function for Demand prediction
#' Function which computes demand forecast.
#' @param i Integer. Iteration
#' @param pred_win Integer. Span of days in prediction window, default is 30.
#' @param pred_lag Integer. Days after initialization day the prediction window begins, default is 15.
#' @param train_y Integer. Years of training data.
#' @param reg_form String. Regression formula before PC covariates have been added.
#' @param p_comps Integer. Number of pc-comps in model.
#' @param other_mods List of other models to run or default is NULL.
#' @param comb Boolean. FALSE allows specific combination of PCs to be run without running all combinations.
#' @param custom String. Specifies a custom pc based model.
#' @param incl_climatology Boolean. TRUE includes climatology model.
#' @param no_pc Boolean. TRUE includes basic (no pc) model.
#' @param gam_lasso Boolean. TRUE sets lasso penalty.
#' @param int Integer. Selects number of climatology models to test.
#' @return data.table
#' @export
Rolling = function(i,X_mat, date_demand, init_days, pred_win, pred_lag, train_y,
                   reg_form, p_comps, other_mods, comb, custom, incl_climatology, no_pc,gam_lasso, int){
    ## **** Step 1: Form the training and test datasets ****
    init_day = init_days[i]
    days_in_m = lubridate::days_in_month(init_day)
    target_days = seq(init_day+ pred_lag, length.out = as.integer(days_in_m),  by = '1 days')
    print(paste('Forecast made on:', init_day))

    train_cutoff = seq(init_day,  length.out = 2, by = paste0('-',train_y, ' year'))[2]

    I_train = date_demand[(date > train_cutoff) & (date < init_day), I]
    I_test = date_demand[date %in% target_days, I]

    if (p_comps > 0){
        pc_data  = get_pca(X_mat, date_demand, I_train, I_test, p_comps) # **** SVD ****

        dt_train = pc_data$dt_train
        dt_test  = pc_data$dt_test

        rm(pc_data)
        gc()
    } else{
        dt_train = date_demand[I_train, ]
        dt_test  = date_demand[I_test, ]
    }

    if (is.null(other_mods) == TRUE){  #Clean-up unless we need X_mat later.
        rm(X_mat)
        gc()
    }

    dt_train= dt_train[!is.na(volume)]
    print(paste('We are training on:', dim(dt_train)[1], 'observations'))
    max_year_train = dt_train[, max(year)]
    max_week_train = dt_train[, max(week)]
    dt_test = dt_test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error
    dt_test = dt_test[week > max_week_train, week := max_week_train] #amounts to pretending that week 53 = week 52  to avoid factor error

    dt_train = dt_train[!is.na(volume)]
    dt_test = dt_test[!is.na(volume)]

    n = nrow(dt_test)

    ## **** Step 2: Train models ****
    mods = list()
    if (no_pc == TRUE){
        print(reg_form)
        mods[[1]] = mgcv::gam(as.formula(reg_form), data = dt_train) #Just Time covariate models
    } else{
        mods[[1]] = mgcv::gam(as.formula("volume ~ 1"), data = dt_train)
    }

    if (p_comps > 0){                              #PC GAM MODELS
        if(custom == FALSE){
            for (j in 2:(p_comps+1)){
                P = paste0("s(PC", 1:(j-1))
                PC = paste(P, collapse=') + ')
                print(sprintf("%s + %s)", reg_form, PC))
                mods[[j]] = mgcv::gam(as.formula(sprintf("%s + %s)", reg_form, PC)), select = gam_lasso, data = dt_train)
            }
        } else{ #Custom order of PC components
            print(sprintf("%s + %s)", reg_form, custom))
            mods[[2]] = mgcv::gam(as.formula(sprintf("%s + %s", reg_form, custom)), data = dt_train)
            }
    }

    ## **** Step 3: Check fit ****
    results = dt_test
    results[,'init_date' := rep(init_day, length(dt_test$date))]
    results[,'lead_time':= ((as.integer(difftime(date, init_day, units = 'days')))*24)+ hour]
    results[,paste0('N_train') := rep(dim(dt_train)[1], dim(dt_test)[1])]
    for (mod in 1:length(mods)){
        pred_demand_test = predict(mods[[mod]], newdata = dt_test)
        pred_demand_train = predict(mods[[mod]], newdata = dt_train)
        RMSE_train = sqrt(mean((dt_train$volume - pred_demand_train)^2))
        RMSE = sqrt(sum((dt_test$volume - pred_demand_test)^2)/n)
        MAE = sum(abs(dt_test$volume - pred_demand_test))/n
        print(sprintf("RMSE mod_%i for initalization %s is %f", mod, format(init_day,"%Y-%m-%d"), RMSE))
        results[,paste0('pred_',mod) := pred_demand_test]
        results[,paste0('RMSE_train_',mod) := rep(RMSE_train, dim(dt_test)[1])]

    }
    ## **** Step 4: Climatology ****
    if(incl_climatology == TRUE){
        climatology = dt_train[, .(clima_pred = mean(volume)), by = .(month_day = format(date, format ="%m-%d"), hour)]

        if ("02-29" %in% format(dt_test$date, format ="%m-%d") & !("02-29" %in%climatology$month_day)){ #Leap year issue
            leap = climatology[month_day=="02-28",]
            leap[,month_day:=  rep("02-29", 24)]
            climatology = rbind(climatology, leap)
        }

        results[,month_day := format(date, format ="%m-%d")]
        results = climatology[results, on = c('month_day', 'hour')]

        print(paste0("Climatology: ", results[,sqrt(sum((volume - clima_pred)^2)/n)] ))

        # for (j in 1:14){
        #     print(j)
        #     climatology = dt_train[, .(ave_vol = mean(volume)), by = c(paste0('yr_by_',j), "hour")]
        #     results = merge(climatology, results, by= c(paste0('yr_by_',j), 'hour') )
        #     setnames(results, "ave_vol", paste0('clim_pred_',j))
        #     print(paste0("Climatology: ", results[,sqrt(mean( (volume- get(paste0('clim_pred_',j)) )^2))]  ))
        # }
    }

    ## **** Step 5: Other Models ****
    if (is.null(other_mods) == FALSE){ #Call other models
        nr = length(mods)
        for (k in 1:length(other_mods)){                                #OTHER MODELS
            mods[[nr+k]] = get(other_mods[[k]])(dt_train, date_demand, I_test, X_mat[I_train,], X_mat[I_test,], m) #tbats needs m, lasso needs X_mat

            if(mods[[nr+k]]$mod$call[1] == "glmnet()"){
                print("glmnet-lasso")
                pred_names = rownames(coef(mods[[nr+k]]))[2:length(rownames(coef(mods[[nr+k]])))]

                if (other_mods[k] == "lasso_temp_and_time2"){
                    dt_test = dt_test[, c("volume", "hour", "month", "year", "season", "week", "w_day") :=.(volume, factor(hour), factor(month), year, season, factor(week), factor(w_day))]
                    dt_test = dt_test[, .(volume, hour, month, year, season, week, w_day)]
                    dt_mod_test = model.matrix(~ hour*week + w_day*month + week + month + season + year + volume, data = dt_test)

                    dat_test = data.table(dt_mod_test, X_mat[I_test,])
                    test = as.matrix(dat_test[, ..pred_names])
                    print(dim(test))
                    pred_demand_test = data.table(predict(mods[[nr+k]], newx = test))
                    pred_demand_test[, c("volume", "date", "hour") := date_demand[I_test, c("volume", "date", "hour")]]
                } else{
                    test_dat = data.table(dt_test, X_mat[I_test,])
                    test = as.matrix(test_dat[,..pred_names])
                    pred_demand_test = data.table(predict(mods[[nr+k]], newx = test))
                    pred_demand_test[, c("volume", "date", "hour") := dt_test[, c("volume", "date", "hour")]]
                }
                results = merge(results, pred_demand_test, by = c("volume","date" ,"hour"))
            } else if(mods[[nr+k]]$mod$call[1] == "xgb.train()"){
                print("xgb - results")
                pred_demand_test = mods[[nr+k]]$pred_demand_test
                results = merge(results, pred_demand_test, by = c("volume","date" ,"hour"))

                #importance_matrix = xgb.importance(model = mods[[mod]])
                #xgb_features[k,] = importance_matrix[1:10,1][[1]]
                #print(importance_matrix[1:10,1][[1]])
                #xgb_gain[k,] = importance_matrix[1:10,2][[1]]
            }
        }
    }

    if (is.data.table(results)==FALSE){
        print(str(results))
        print('Trying to convert')
        results = as.data.table(results)
        if(dim(results)[1] < 100){
            print('This is strange')
            print(results)
        }
    }
    return(results)
}




#' Basic Model function for Demand prediction
#' Function which computes demand forecast.
#' @param i Integer. Iteration
#' @param pred_win Integer. Span of days in prediction window, default is 30.
#' @param pred_lag Integer. Days after initialization day the prediction window begins, default is 15.
#' @param train_y Integer. Years of training data.
#' @param reg_form String. Regression formula before PC covariates have been added.
#' @param p_comps Integer. Number of pc-comps in model.
#' @param other_mods List of other models to run or default is NULL.
#' @param comb Boolean. FALSE allows specific combination of PCs to be run without running all combinations.
#' @param custom String. Specifies a custom pc based model.
#' @param incl_climatology Boolean. TRUE includes climatology model.
#' @param no_pc Boolean. TRUE includes basic (no pc) model.
#' @param gam_lasso Boolean. TRUE sets lasso penalty.
#' @param int Integer. Selects number of climatology models to test.
#' @return data.table
#' @export
Rolling_time = function(i,X_mat, date_demand, init_days, pred_win, pred_lag, train_y,
                   reg_form, p_comps, other_mods, comb, custom, incl_climatology, no_pc,gam_lasso, int){
    ## **** Step 1: Form the training and test datasets ****
    init_day = init_days[i]
    target_days = seq(init_day+ pred_lag, length.out = pred_win,  by = '1 days')
    print(paste('Forecast made on:', init_day))

    train_cutoff = seq(init_day,  length.out = 2, by = paste0('-',train_y, ' year'))[2]

    I_train = date_demand[(date > train_cutoff) & (date < init_day), I]
    I_test = date_demand[date %in% target_days, I]

    if (p_comps > 0){
        pc_data  = get_pca(X_mat, date_demand, I_train, I_test, p_comps) # **** SVD ****

        dt_train = pc_data$dt_train
        dt_test  = pc_data$dt_test

        rm(pc_data)
        gc()
    } else{
        dt_train = date_demand[I_train, ]
        dt_test  = date_demand[I_test, ]
    }

    if (is.null(other_mods) == TRUE){  #Clean-up unless we need X_mat later.
        rm(X_mat)
        gc()
    }

    dt_train= dt_train[!is.na(volume)]
    print(paste('We are training on:', dim(dt_train)[1], 'observations', 'on iteration ', i, 'starting on', format(init_day,"%Y-%m-%d")))
    max_year_train = dt_train[, max(year)]
    max_week_train = dt_train[, max(week)]
    dt_test = dt_test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error
    dt_test = dt_test[week > max_week_train, week := max_week_train] #amounts to pretending that week 53 = week 52  to avoid factor error
    dt_train = dt_train[!is.na(volume)]
    dt_test = dt_test[!is.na(volume)]

    results = dt_test
    results[,'init_date' := rep(init_day, length(dt_test$date))]
    results[,'lead_time':= ((as.integer(difftime(date, init_day, units = 'days')))*24)+ hour]
    results[, paste0('N_train') := rep(dim(dt_train)[1], dim(dt_test)[1])]

    ## **** Step 2: Train models ****
    mods = list()
    print(reg_form)
    mods[[1]] = mgcv::gam(as.formula(reg_form), data = dt_train) #Just Time covariate models
    dt_train[, pred_demand_train := predict(mods[[1]], newdata=.SD),]
    dt_train[, res := volume - pred_demand_train]

    results[, paste0('pred_',1) := predict(mods[[1]], newdata=.SD),]
    results[, res := volume - get(paste0('pred_',1))]

    MAE = results[,mean(abs(volume - get(paste0('pred_',1))))]
    MAPE = results[,mean( abs((volume - get(paste0('pred_',1))) /volume))]
    RMSE = results[,sqrt(mean((volume - get(paste0('pred_',1)))^2))]

    print(sprintf("For mod_%i on initalization %s RMSE is %.2f, MAE is %.2f, MAPE is %.3f", 1, format(init_day,"%Y-%m-%d"), RMSE,MAE,MAPE))

    #Custom order of PC components
    print(custom)
    mods[[2]] = mgcv::gam(as.formula(custom), data = dt_train)

    results[,res_pred := predict(mods[[2]], newdata =.SD),]
    results[,paste0('pred_',2) := res_pred + pred_1]

    MAE = results[,mean(abs(volume - get(paste0('pred_',2))))]
    MAPE = results[,mean( abs((volume - get(paste0('pred_',2))) /volume))]
    RMSE = results[,sqrt(mean((volume - get(paste0('pred_',2)))^2))]
    print(sprintf("For mod_%i on initalization %s RMSE is %.2f, MAE is %.2f, MAPE is %.3f", 2, format(init_day,"%Y-%m-%d"), RMSE,MAE,MAPE))


    if (is.data.table(results)==FALSE){
        print(str(results))
        print('Trying to convert')
        results = as.data.table(results)
        if(dim(results)[1] < 100){
            print('This is strange')
            print(results)
        }
    }
    return(results)
}



#' Basic Model function for Demand prediction
#' Function which computes demand forecast based on anomalies
#' @param i Integer. Iteration
#' @param pred_win Integer. Span of days in prediction window, default is 30.
#' @param pred_lag Integer. Days after initialization day the prediction window begins, default is 15.
#' @param train_y Integer. Years of training data.
#' @param reg_form String. Regression formula before PC covariates have been added.
#' @param p_comps Integer. Number of pc-comps in model.
#' @param other_mods List of other models to run or default is NULL.
#' @param comb Boolean. FALSE allows specific combination of PCs to be run without running all combinations.
#' @param custom String. Specifies a custom pc based model.
#' @param incl_climatology Boolean. TRUE includes climatology model.
#' @param no_pc Boolean. TRUE includes basic (no pc) model.
#' @param gam_lasso Boolean. TRUE sets lasso penalty.
#' @param int Integer. Selects number of climatology models to test.
#' @return data.table
#' @export
Rolling_anomaly = function(i,X_mat, date_demand, init_days,pred_win, pred_lag, train_y,
                        reg_form, p_comps, other_mods, comb, custom, incl_climatology, no_pc,gam_lasso, int){
    ## **** Step 1: Form the training and test datasets ****
    init_day = init_days[i]
    target_days = seq(init_day+ pred_lag, length.out = pred_win,  by = '1 days')
    print(paste('Forecast made on:', init_day))

    train_cutoff = seq(init_day,  length.out = 2, by = paste0('-',train_y, ' year'))[2]

    I_train = date_demand[(date > train_cutoff) & (date < init_day), I]
    I_test = date_demand[date %in% target_days, I]

    if (p_comps > 0){
        pc_data  = get_pca(X_mat, date_demand, I_train, I_test, p_comps) # **** SVD ****

        dt_train = pc_data$dt_train
        dt_test  = pc_data$dt_test

        rm(pc_data)
        gc()
    } else{
        dt_train = date_demand[I_train, ]
        dt_test  = date_demand[I_test, ]
    }

    if (is.null(other_mods) == TRUE){  #Clean-up unless we need X_mat later.
        rm(X_mat)
        gc()
    }

    dt_train= dt_train[!is.na(volume)]
    print(paste('We are training on:', dim(dt_train)[1], 'observations', 'on iteration ', i, 'starting on', format(init_day,"%Y-%m-%d")))
    max_year_train = dt_train[, max(year)]
    max_week_train = dt_train[, max(week)]
    dt_test = dt_test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error
    dt_test = dt_test[week > max_week_train, week := max_week_train] #amounts to pretending that week 53 = week 52  to avoid factor error
    dt_train = dt_train[!is.na(volume)]
    dt_test = dt_test[!is.na(volume)]

    results = dt_test
    results[,'init_date' := rep(init_day, length(dt_test$date))]
    results[,'lead_time':= ((as.integer(difftime(date, init_day, units = 'days')))*24)+ hour]
    results[, paste0('N_train') := rep(dim(dt_train)[1], dim(dt_test)[1])]

    ## **** Step 2: Train models ****
    mods = list()
    print(reg_form)
    mods[[1]] = mgcv::gam(as.formula(reg_form), data = dt_train) #Just Time covariate models
    dt_train[, pred_demand_train := predict(mods[[1]], newdata=.SD),]
    dt_train[, res := volume - pred_demand_train]

    results[, paste0('pred_',1) := predict(mods[[1]], newdata=.SD),]
    results[, res := volume - get(paste0('pred_',1))]

    MAE = results[,mean(abs(volume - get(paste0('pred_',1))))]
    MAPE = results[,mean( abs((volume - get(paste0('pred_',1))) /volume))]
    RMSE = results[,sqrt(mean((volume - get(paste0('pred_',1)))^2))]

    print(sprintf("For mod_%i on initalization %s RMSE is %.2f, MAE is %.2f, MAPE is %.3f", 1, format(init_day,"%Y-%m-%d"), RMSE,MAE,MAPE))

    #Custom order of PC components
    print(custom)

    dt_train[,PC1_clim := mean(PC1), by =c('year_day', 'hour')]
    dt_train[,PC1_anomal := PC1 - PC1_clim]
    dt_train[,PC1_anom_stand:=(PC1_anomal-mean(PC1_anomal))/sd(PC1_anomal)]

    mods[[2]] = mgcv::gam(as.formula(custom), data = dt_train)


    climatology = dt_train[,.(PC1_clim = mean(PC1)), by =c('year_day', 'hour')]
    results = merge(climatology, results, by =c('year_day', 'hour') )
    results[,PC1_anomal := PC1 - PC1_clim]
    meanPC1_anomal = dt_train[,mean(PC1_anomal)]
    sdPC1_anomal = dt_train[,sd(PC1_anomal)]
    results[,PC1_anom_stand:=(PC1_anomal-meanPC1_anomal)/sdPC1_anomal]



    results[,res_pred := predict(mods[[2]], newdata =.SD),]
    results[,paste0('pred_',2) := res_pred + pred_1]

    MAE = results[,mean(abs(volume - get(paste0('pred_',2))))]
    MAPE = results[,mean( abs((volume - get(paste0('pred_',2))) /volume))]
    RMSE = results[,sqrt(mean((volume - get(paste0('pred_',2)))^2))]
    print(sprintf("For mod_%i on initalization %s RMSE is %.2f, MAE is %.2f, MAPE is %.3f", 2, format(init_day,"%Y-%m-%d"), RMSE,MAE,MAPE))


    if (is.data.table(results)==FALSE){
        print(str(results))
        print('Trying to convert')
        results = as.data.table(results)
        if(dim(results)[1] < 100){
            print('This is strange')
            print(results)
        }
    }
    return(results)
}


#' Test function for Climatology
#' Function which computes demand forecast.
#' @param i Integer. Iteration
#' @param pred_win Integer. Span of days in prediction window, default is 30.
#' @param pred_lag Integer. Days after initialization day the prediction window begins, default is 15.
#' @param train_y Integer. Years of training data.
#' @param reg_form String. Regression formula before PC covariates have been added.
#' @param p_comps Integer. Number of pc-comps in model.
#' @param other_mods List of other models to run or default is NULL.
#' @param comb Boolean. FALSE allows specific combination of PCs to be run without running all combinations.
#' @param custom Boolean. Specifies a custom pc based model.
#' @param incl_climatology Boolean. TRUE includes climatology model.
#' @param no_pc Boolean. TRUE includes basic (no pc) model.
#' @param gam_lasso Boolean. TRUE sets lasso penalty.
#' @param int Integer. Selects number of climatology models to test.
#' @return data.table
#' @export
Rolling_test = function(i, X_mat, date_demand, init_days,pred_win, pred_lag, train_y,
                        reg_form, p_comps, other_mods, comb, custom,incl_climatology, no_pc,gam_lasso, int){
    init_day = init_days[i]
    target_days = seq(init_day+ pred_lag, length.out = pred_win,  by = '1 days')
    print(paste('Forecast made on:', init_day))
    train_cutoff = seq(init_day,  length.out = 2, by = paste0('-',train_y, ' year'))[2]
    I_train = date_demand[(date > train_cutoff) & (date < init_day), I]
    I_test = date_demand[date %in% target_days, I]
    print(length(I_train))
    if (p_comps > 0){
        pc_data  = get_pca(X_mat, date_demand, I_train, I_test, p_comps) # **** SVD ****
        rm(X_mat)
        gc()
        dt_train = pc_data$dt_train
        dt_test  = pc_data$dt_test
    } else{
        dt_train = date_demand[I_train, ]
        dt_test  = date_demand[I_test, ]
    }
    dt_train = dt_train[!is.na(volume)]
    dt_test = dt_test[!is.na(volume)]
    print(paste('We are training on:', dim(dt_train)[1], 'observations ', 'on iteration ', i))
    max_year_train = dt_train[, max(year)]
    dt_test = dt_test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error
    n = nrow(dt_test)
    results = dt_test
    results[,'init_date' := rep(init_day, length(dt_test$date))]
    results[,'lead_time':= ((as.integer(difftime(date, init_day, units = 'days')))*24)+ hour]
    results[,paste0('N_train') := rep(dim(dt_train)[1], dim(dt_test)[1])]
    if(incl_climatology){

        results = Rolling_climatology(results, dt_train, dt_test, reg_form,
                                  intervals = int, basic = TRUE)
        print(paste0('Completed climatology ', i, ' ', Sys.time()))
    }
    return(results)
}

#' Test function for Climatology2
#' Function which computes demand forecast.
#' @param results Integer. Iteration
#' @param dt_train Integer. Span of days in prediction window, default is 30.
#' @param dt_test Integer. Days after initialization day the prediction window begins, default is 15.
#' @param reg_form String. Regression formula before PC covariates have been added.
#' @param intervals Integer. Number of model intervals to test.
#' @param basic Boolean
#' @param daily_cycle Boolean. FALSE allows specific combination of PCs to be run without running all combinations.
#' @param weekly_cycle Boolean.
#' @return data.table
#' @export
Rolling_climatology = function(results, dt_train, dt_test, reg_form,
                           intervals = 7, basic = TRUE, daily_cycle = FALSE, weekday_cycle = FALSE){
    for (j in 7:intervals){

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

        ## **** Step 3:Run Models  ****
        mod_clim1 = mgcv::gam(as.formula(paste0(reg_form," + clim_pred_",j) ), data = dt_train)


        ## **** Step 4:Predict and log results  ****
        pred_demand_test_clim1 = predict(mod_clim1, newdata = results)
        results[,paste0('pred_','mod_clim1_',j) := pred_demand_test_clim1]

        ## **** Step 5:Find Weekday cycle  ****
        if (j>=7){
            dt_train[, ave_vol_wd := mean(volume), by = c(paste0('yr_by_',j), "w_day")]
            climatology = dt_train[, .(ave_vol_wd = mean(volume)), by = c(paste0('yr_by_',j), "w_day")]
            results = merge(climatology, results, by= c(paste0('yr_by_',j), 'w_day') )
            setnames(dt_train, "ave_vol_wd", paste0('ave_vol_wd_',j))
            setnames(results, "ave_vol_wd", paste0('ave_vol_wd_',j))
            mod_clim2 = mgcv::gam(as.formula(paste0(reg_form, "+ ave_vol_wd_", j ) ) , data = dt_train)
            mod_clim3 = mgcv::gam(as.formula(paste0(reg_form," + clim_pred_",j," + ave_vol_wd_", j ) ) , data = dt_train)
            pred_demand_test_clim2 = predict(mod_clim2, newdata = results)
            pred_demand_test_clim3 = predict(mod_clim3, newdata = results)
            results[,paste0('pred_','mod_clim2_',j) := pred_demand_test_clim2]
            results[,paste0('pred_','mod_clim3_',j) := pred_demand_test_clim3]
        }


        #print(paste0("PC1 + CLIM: ",  round(results[,sqrt(mean( (volume- get(paste0('pred_','mod_clim_',j)) )^2))],2)))
        #print(paste0("Climatology: ", round(results[,sqrt(mean( (volume- get(paste0('clim_pred_',j       )) )^2))],2)))

    }
    return(results)

}
