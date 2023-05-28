
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
#'
#' @return A data.table with actual and predicted demand values.
#'
#' @examples  Mod1 = demand_forecast_final(forc_start ='2014-01-01',forc_end = '2021-03-01',pred_win = 30, train_y = 5, pred_lag = 15,reg_form = "volume ~ as.factor(hour) + as.factor(month) + year", p_comps = 3, other_mods= NULL,comb = TRUE)
#' @export
demand_forecast_final = function(X_mat, date_demand, NWP_pred, forc_start, forc_end, pred_win = 45, pred_lag= 0, train_y=3,
                           reg_form, p_comps, other_mods= NULL, comb = TRUE, custom = FALSE,
                           incl_climatology = FALSE, no_pc = TRUE, cores = 44){

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
                                "Rolling_final",
                                X_mat = X_mat, date_demand = date_demand,init_days= init_days_all, pred_win = pred_win, pred_lag= pred_lag, train_y=train_y,
                                reg_form= reg_form, p_comps= p_comps, other_mods= other_mods, comb = comb, custom = custom,
                                incl_climatology = incl_climatology, no_pc = no_pc,gam_lasso =gam_lasso,NWP_pred = NWP_pred,
                                mc.cores = cores)

    ## **** Return results ****
    Results = rbindlist(detailed_results)
    out = list()
    out$Results = Results
    print('Demand forecast has completed')
    return(out)
}

#' @export
Rolling_final = function(i,X_mat, date_demand, init_days,pred_win, pred_lag, train_y,
                   reg_form, p_comps, other_mods, comb, custom,incl_climatology, no_pc,num_thread,gam_lasso,NWP_pred = NWP_pred){
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

    if (p_comps > 0){                              #PC GAM MODELS
            mod = mgcv::gam(as.formula(sprintf("%s + %s", reg_form, custom)), data = dt_train)
    }

    ## **** Step 3: Check fit ****
    dt_test = merge(dt_test, NWP_pred, by = c("hour", "date"))
    dt_PC = dt_test[,.(date, hour, PC1, PC2)]
    results = dt_test
    results[,'lead_time':= (((as.integer(difftime(date, init_date, units = 'days')))*24)+ hour)/6]

    pred_demand_PC1_actual = predict(mod, newdata = dt_test)
    results[,'pred_demand_PC1_actual' := pred_demand_PC1_actual]
    dt_test[, PC1:= NULL]
    dt_test[, PC2:= NULL]
    pred_names1 = c("mean","q10","q20","q30","q40",
                    "q50","q60","q70","q80","q90")
    for (name in pred_names1){ #Fold all quantiles into the prediction
        PC1_name = paste0("NWP_PC1_",name)
        PC2_name = paste0("NWP_PC2_",name)
        setnames(dt_test, old = PC1_name, new = "PC1")
        setnames(dt_test, old = PC2_name, new = "PC2")
        pred_demand_test = predict(mod, newdata = dt_test)
        results[,paste0('pred_',name) := pred_demand_test]
        dt_test[, PC1:= NULL]
        dt_test[, PC2:= NULL]
    }
    final_results = merge(results, dt_PC, by = c("hour", "date"))

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

