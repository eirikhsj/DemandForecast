
#' Rolling CV forecast model for energy demand.
#'
#' @param forc_start Date for first initalization date.
#' @param forc_end Date for last initialization date.
#' @param pred_win Span of days in prediction window, default is 30.
#' @param pred_lag Days after initialization day the prediction window begins, default is 15.
#' @param train_y Years of training data.
#' @param reg_form Regression formula before PC covariates have been added.
#' @param p_comps Integer. Number of pc-comps in model.
#' @param other_mods List of other models to run or default is NULL.
#' @param comb Boolean. FALSE allows specific combination of PCs to be run without running all combinations.
#' @param custom Boolean.Specifies a custom pc based model.
#' @param incl_climatology Boolean.Includes climatology model.
#' @param no_pc Boolean.Includes basic (no pc) model.
#'
#' @return A data.table with actual and predicted demand values.
#' @export
#' @import mgcv
#' @import forecast
#' @import xgboost
#' @import glmnet
#' @import parallel
#' @import data.table
#' @examples  Mod1 = demand_forecast(forc_start ='2014-01-01',forc_end = '2021-03-01',pred_win = 30, train_y = 5, pred_lag = 15,reg_form = "volume ~ as.factor(hour) + as.factor(month) + year", p_comps = 3, other_mods= NULL,comb = TRUE)
#'
#'
demand_forecast = function(X_mat, date_demand, forc_start, forc_end, pred_win = 30, pred_lag= 15, train_y=5,
                           reg_form, p_comps, other_mods= NULL, comb = TRUE, custom = FALSE,
                           incl_climatology = FALSE, no_pc = TRUE, cores = 4){

    last_poss_pred = range(date_demand$date)[2] - pred_win - pred_lag

    #1) Fix dates
    start = as.Date(forc_start)
    end = as.Date(forc_end)
    all_days = seq(start, end,  by = '1 days')
    init_days_all = all_days[mday(all_days)==1]

    ## **** Run parallel cores ****
    blas_set_num_threads(1)
    omp_set_num_threads(1) #Set number of threads

    detailed_results = mclapply(seq_along(init_days_all),
                      "Rolling",
                      X_mat = X_mat, date_demand = date_demand,init_days= init_days_all, pred_win = pred_win, pred_lag= pred_lag, train_y=train_y,
                      reg_form= reg_form, p_comps= p_comps, other_mods= other_mods, comb = comb, custom = custom,
                      incl_climatology = incl_climatology, no_pc = no_pc,
                      mc.cores = cores)

    ## **** Return results ****
    Results = rbindlist(detailed_results)
    out = list()
    out$Results = Results
    print('Demand forecast has completed')
    return(out)
}

#' @export
Rolling = function(i,X_mat, date_demand, init_days,pred_win, pred_lag, train_y,
                   reg_form, p_comps, other_mods, comb, custom,incl_climatology, no_pc,num_thread){
    ## ***** Step 1: Form the training and test datasets ****
    init_day = init_days[i]
    target_days = seq(init_day+ pred_lag, length.out = pred_win,  by = '1 days')
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
                mods[[j]] = mgcv::gam(as.formula(sprintf("%s + %s)", reg_form, PC)), data = dt_train)
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
        climatology = dt_train[, .(ave_vol = mean(volume)), by = .(month_day = format(date, format ="%m-%d"), hour)]

        if ("02-29" %in% format(dt_test$date, format ="%m-%d") & !("02-29" %in%climatology$month_day)){ #Leap year issue
            leap = climatology[month_day=="02-28",]
            leap[,month_day:=  rep("02-29", 24)]
            climatology = rbind(climatology, leap)
        }
        clima_pred = climatology[month_day %in% format(dt_test$date, format ="%m-%d"), ave_vol]
        results[,'clima_pred' := clima_pred]
        print(paste0("Climatology: ", sqrt(sum((dt_test$volume - clima_pred)^2)/n)))
    }

    ## **** Step 5: Other Models ****
    if (is.null(other_mods) == FALSE){ #Call other models
        nr = length(mods)
        for (k in 1:length(other_mods)){                                #OTHER MODELS
            mods[[nr+k]] = get(other_mods[[k]])(dt_train, X_mat[I_train,], m) #tbats needs m, lasso needs X_mat

            if(mods[[nr+k]]$call[1] == "glmnet()"){
                print("glmnet-lasso")
                pred_names = rownames(coef(mods[[nr+k]]))[2:length(rownames(coef(mods[[nr+k]])))]

                if (other_mods[k] == "lasso_temp_and_time"){
                    dt_test = date_demand[,.(volume, hour, month, year, season, week, w_day)]
                    dt_test$hour = factor(dt_test$hour)
                    dt_test$w_day = factor(dt_test$w_day)
                    dt_test$week = factor(dt_test$week)
                    dt_test$month = factor(dt_test$month)
                    dt_test$season = factor(dt_test$season)
                    dt_test = dt_test[I_test,]
                    dt_mod_test = model.matrix(~ . - 1, data = dt_test)

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

#' @export
lasso_just_temp = function(dt_train, X_mat, m){
    #ls = sort(round(exp((1:1000)/2000)^14 -1, 2), decreasing = TRUE)
    #ls = sort(round(exp((1:1000)/1000)^3 +1, 2), decreasing = TRUE)
    ls = sort(seq(3,5,length.out = 1000), decreasing = TRUE)
    #ls = sort(round(seq(0,20, length.out = 1000), 2), decreasing = TRUE)

    dat = data.table(dt_train[,.(volume)], X_mat)
    Y_inp = as.matrix(dat[,volume])
    X_inp = as.matrix(dat[, .SD, .SDcols = -'volume'])
    set.seed(100)
    l1 = glmnet(x = X_inp, y =Y_inp, alpha = 1, lambda = ls)

    return(l1)
}

#' @export
lasso_temp_and_time = function(dt_train, X_mat, m){
    #ls = sort(round(exp((1:1000)/2000)^14 -1, 2), decreasing = TRUE)
    ls = sort(round(seq(0,20, length.out = 1000), 2), decreasing = TRUE)
    #ls = sort(seq(0,8,length.out = 1000), decreasing = TRUE)
    dt_train = dt_train[,.(volume, hour, month, year, season, week, w_day)]
    dt_train$hour = factor(dt_train$hour)
    dt_train$w_day = factor(dt_train$w_day)
    dt_train$week = factor(dt_train$week)
    dt_train$month = factor(dt_train$month)
    dt_train$season = factor(dt_train$season)
    dt_mod_train = model.matrix(~ . - 1, data = dt_train)
    dat_train = data.table(dt_mod_train, X_mat)

    Y_inp = as.matrix(dat_train[,volume])
    X_inp = as.matrix(dat_train[, .SD, .SDcols = -'volume'])
    set.seed(100)
    l1 = glmnet(x = X_inp, y =Y_inp, alpha = 1, lambda = ls)

    return(l1)
}


# res = pred_demand_test[, lapply(.SD, function(x) sqrt(mean((volume - x)^2)))]
# pred_demand_test = pred_demand_test[,1]
