#' @title Rolling Temperature PC model using NWP forecasts
#'
#' @param forc_start Date. Quantile of interest.
#' @param forc_end Date. Quantile of interest.
#' @param q Float. Quantile of interest.
#' @param ERA_NWP Data.table. Data of ERA observations and quantiles of NWP forecasts.
#' @param predictors Integer. Number of NWP predictors.
#' @param model String. Specify type of model, e.g. qreg, qgam, etc.
#' @param window Integer. Max number of forecast days from initialization.
#' @param reweight Boolean. Use reweighted data.
#' @param incl_climatology Boolean. Include climatology model.
#' @param incl_others Boolean. Include climatology model.
#' @param cores Integer. Number of cores to be used.
#' @param formula String. Regression formula.
#' @param df_spline Integer. Number of basis splines.
#'
#' @return Model output with predictions, loss and model coefficients.
#' @import quantreg
#' @import data.table
#' @import splines
#' @import stats
#'
#' @examples mod = rolling_mod_NWP('2007-01-01', '2022-05-01', 0.9, ERA_NWP, 1, model = 'qreg', window = 60, reweight=FALSE)
#' @name rolling_mod_NWP

# forc_start=as.Date('2012-01-01'); forc_end=as.Date('2023-01-01'); q=0.9; ERA_NWP; predictors=1; model='spline'; window = 125; reweight = FALSE;
# incl_climatology =TRUE; formula = 'PC1 ~ 1'; incl_other = FALSE; skill_interval = 0; cores = 4;df_spline = 8

#' @export
rolling_mod_NWP = function(forc_start=as.Date('2007-01-01'), forc_end=as.Date('2023-01-01'), q=0.9, ERA_NWP,
                           formula = 'PC1 ~ 1', predictors=1, model='qreg', window = 60, reweight = FALSE,
                           incl_climatology =FALSE, incl_other = FALSE, skill_interval = 0, cores = 4, df_spline = 8){

    #1) Fix dates
    start =as.Date(forc_start)
    end = as.Date(forc_end)
    all_days = seq(start, end,  by = '1 days')

    if (reweight==TRUE){
        init_days = all_days[mday(all_days)==16]  #reweighting starts on the 15th
        re = 're'
    }else{
        init_days = all_days[mday(all_days)==1]
        re = ''
    }

    #2) Build regression formula and fetch variables
    incl_vars = c('date','PC1', 'PC2','hour', 'week', 'month', 'season', 'year', 'init_date', 'lead_time')
    pred_vars = c()

    if (predictors > 0){
        for (l in 1:predictors){
            pred_vars = c(pred_vars, paste0('NWP', l, '_', re, q*100))
            formula = paste0(formula, ' + ', paste0('NWP', l, '_', re, q*100))}
    }
    if (incl_other[1] == FALSE){
        incl_vars = c(incl_vars, pred_vars)
    } else{
        pred_vars = incl_other[2]
        incl_vars = c(incl_vars, pred_vars, incl_other[1])
        print(incl_vars)
    }

    ERA_NWP_vars = ERA_NWP[lead_time <= window*4,.SD, .SDcols =incl_vars]

    if (skill_interval>0){
        ERA_NWP_vars[, c(paste0("PC_", skill_interval, "days_roll"), paste0("NWP1_roll", skill_interval)):=
                                            .(rollapply(PC1, width = skill_interval*4, partial = TRUE, FUN = mean),
                                              rollapply(NWP1_90, width = skill_interval*4, partial = TRUE, FUN = mean)), by = .(init_date)]
    }


    ## **** Run parallel cores ****
    blas_set_num_threads(1)
    omp_set_num_threads(1) #Set number of threads

    #3) Forecast iteration
    detailed_results = mclapply(seq_along(init_days),
                                "Rolling_nwp",
                                ERA_NWP_vars = ERA_NWP_vars,q= q, init_days= init_days, window = window, reweight = reweight, model = model,
                                predictors = predictors, incl_climatology = incl_climatology, formula = formula, pred_vars = pred_vars, df_spline= df_spline,
                                mc.cores = cores)

    #4) Store and return
    Results = rbindlist(detailed_results,use.names=FALSE)
    out = c()
    out$Results = Results
    print('Temperature forecast has completed')
    print('温度预报完成')
    return(out)
}

#' @export
Rolling_nwp = function(i, ERA_NWP_vars, q, init_days, window, reweight, model, predictors,
                       incl_climatology, formula, pred_vars, df_spline){
    ## 3a) Time keeping
    init_day = init_days[i]
    target_days = seq(init_day, length.out = window,  by = '1 days')
    print(paste('Forecast made on:', init_day))
    ERA_NWP_time = ERA_NWP_vars[date <= target_days[length(target_days)],]
    ERA_NWP_final= na.omit(ERA_NWP_time)

    ## 3b) Split train-test
    train = ERA_NWP_final[date<init_day, .SD, keyby = .(date,hour)]
    if (reweight ==TRUE){
        test = ERA_NWP_final[date %in%target_days, .SD, keyby = .(date,hour)] #Here we only use 1 predictor, cant use init_day but no confusion in target_days
    } else{
        test = ERA_NWP_final[init_date == init_day, .SD, keyby = .(date,hour)] #Use of several preds means target_days selects to many dates
    }
    results = test
    # if (max(train$year) > max(test$year)){
    #     test = test[year > max(train$year), year := max(train$year)]
    # }

    ## 3c) Run qr-reg
    if (model == 'qreg'){
        print(formula)
        qreg = rq(formula, data = train, tau = c(q))
        #qreg = gam(as.formula(formula), data = train)
        train_l = pinball_loss(q, predict(qreg), train$PC1)
        test_l = pinball_loss(q, predict(qreg, newdata = test), test$PC1)
        results[,'pred' := predict(qreg, newdata = test)]

    }  else if(model == 'spline'){
        print(formula)
        spline_var = paste0("NWP1_",q*100)
        spline_form = paste0("test$PC1 ~ bs(test$",spline_var,", df=",df_spline,")")
        X_test = model.matrix(as.formula(spline_form))
        qreg = quantreg::rq(PC1 ~ splines::bs(get(spline_var), df=df_spline), data=train, tau=c(q))
        test_l = pinball_loss(q,  predict(qreg, newdata = test), test$PC1)
        results[, 'pred':= predict(qreg, newdata = test)]
    }

    ## 3d) Register loss
    results[,'test_loss' := test_l]
    print(paste0('Ave pinball loss for model issued on ', init_day, ' is = ', round(mean(test_l),digits = 2)))

    ## 3e) Register Beta coefficients

    if(predictors >0 & model == 'qreg'){
        betas = data.table(t(coef(qreg)))[,..pred_vars]
        colnames(betas) = paste0('coef_',pred_vars)
        betas = betas[rep(1,dim(test)[1]),] #unnecessary storage use, expand later
        results = cbind(results, betas)
    }

    ## 3f) Include Climatology
    if(incl_climatology == TRUE){
        if (skill_interval>0){
            pc_int = paste0("PC_", skill_interval, "days_roll")
            climatology = train[, .(quant =quantile(get(pc_int),probs = q)), by = .(month_day = format(date, format ="%m-%d"), hour)]
        } else{
            climatology = train[, .(quant =quantile(PC1,probs = q)), by = .(month_day = format(date, format ="%m-%d"), hour)]
        }

        if ("02-29" %in% format(test$date, format ="%m-%d") & !("02-29" %in%climatology$month_day)){ #Leap year issue
            print('Correcting leap year')
            leap = climatology[month_day=="02-28",]
            leap[,month_day:=  rep("02-29", 4)]
            climatology = rbind(climatology, leap)
        }
        pred_clima = climatology[month_day %in% format(test$date, format ="%m-%d"), .(month_day, hour,  clima_pred= quant)]
        results[,month_day:=  format(date, format ="%m-%d")]
        results = merge(results,pred_clima, by = c("month_day", "hour"))

    }

    return(results)
}

