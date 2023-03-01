
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
#' @param custom Boolean.Specifies a custom pc bases model.
#' @param incl_climatology Boolean.Includes climatology model.
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
demand_forecast = function(X_mat1, date_demand1, forc_start1, forc_end1, pred_win1 = 30, pred_lag1= 15, train_y1=5,
                           reg_form1, p_comps1, other_mods1= NULL, comb1 = TRUE, custom1 = FALSE,
                           incl_climatology1 = TRUE, cores = 4){

    last_poss_pred = range(date_demand$date)[2] - pred_win - pred_lag

    #1) Fix dates
    start = as.Date(forc_start)
    end = as.Date(forc_end)
    all_days = seq(start, end,  by = '1 days')
    init_days1 = all_days[mday(all_days)==1]


    ## **** Run parallel cores ****
    detailed_results = mclapply(seq_along(init_days1),
                      "Rolling",
                      X_mat = X_mat1, date_demand = date_demand1,init_days= init_days1, pred_win = pred_win1, pred_lag= pred_lag1, train_y=train_y1,
                      reg_form= reg_form1, p_comps= p_comps1, other_mods= other_mods1, comb = comb1, custom = custom1,
                      incl_climatology = incl_climatology1,
                      mc.cores = cores)

    # detailed_results = lapply(seq_along(init_days),
    #                             "Rolling",
    #                             X_mat = X_mat, date_demand = date_demand,init_days= init_days, pred_win = pred_win, pred_lag= pred_lag, train_y=train_y,
    #                             reg_form= reg_form, p_comps= p_comps, other_mods= other_mods, comb = comb, custom = custom,
    #                             incl_climatology = incl_climatology)

    Results = rbindlist(detailed_results)

    out = list()
    out$Results = Results
    #out$mods = mods
    print('Demand forecast has completed')
    return(out)
}

#' @export
Rolling = function(i,X_mat, date_demand, init_days,pred_win, pred_lag, train_y,
                   reg_form, p_comps, other_mods, comb, custom,incl_climatology
){

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
    #print(paste('We are training on:', dim(dt_train)[1], 'observations'))
    max_year_train = dt_train[, max(year)]
    dt_test = dt_test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error
    n = nrow(dt_test)

    ## ***** Step 2: Train models ****
    mods = list()
    mods[[1]] = gam(as.formula(reg_form), data = dt_train) #Just Time covariate models

    j = 1
    if (p_comps > 0){                              #PC GAM MODELS
        if(custom == FALSE){
            for (j in 2:(p_comps+1)){
                P = paste0("s(PC", 1:(j-1))
                PC = paste(P, collapse=') + ')
                mods[[j]] = gam(as.formula(sprintf("%s + %s)", reg_form, PC)), data = dt_train)
            }
        } else{ mods[[2]] = gam(as.formula(sprintf("%s + %s", reg_form, custom)), data = dt_train) }
    }

    ## **** Step 3: Check fit ****
    #print(paste0("Test : ", dt_test[1,date], " - ", dt_test[length(dt_test$date),date], " "))
    results = dt_test
    results[,'init_date' := rep(init_day, length(dt_test$date))]
    results[,'lead_time':= ((as.integer(difftime(date, init_day, units = 'days')))*24)+ hour]
    for (mod in 1:length(mods)){
        pred_demand_test = predict(mods[[mod]], newdata = dt_test)
        RMSE = sqrt(sum((dt_test$volume - pred_demand_test)^2)/n)
        MAE = sum(abs(dt_test$volume - pred_demand_test))/n
        #print(sprintf("RMSE mod_%i for initalization %s is %f", mod, format(init_day,"%Y-%m-%d"), RMSE))
        results[,paste0('pred_',mod) := pred_demand_test]
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
        #print(paste0("Climatology: ", sqrt(sum((dt_test$volume - clima_pred)^2)/n)))
    }


    #detailed_results[[i]] = results
    #diff = round(difftime(Sys.time(),start_time, units="secs"),4)
    #print(paste0('This round took ', diff, ' seconds. Only ', round((length(init_days)-i)*diff/60,3), ' minutes left'))
    return(results)
}