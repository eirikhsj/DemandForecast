#' @title Rolling model
#'
#' @param q Quantile of interest.
#' @param ERA_NWP Data of ERA observations and quantiles of NWP forecasts.
#' @param predictors Number of NWP predictors.
#' @param model Specify type of model, e.g. qreg, qgam, etc.
#' @param window Max number of forecast days from initialization.
#' @param hour_v Boolean. Include hour covariate.
#' @param week_v Boolean. Include week covariate.
#' @param month_v Boolean. Include month covariate.
#' @param year_v Boolean. Include year covariate.
#'
#' @return Model output with predictions, loss and model coefficients.
#' @export
#' @import quantreg
#' @import data.table
#' @import stats
#'
#' @examples mod = rolling_mod('2007-01-01', '2022-05-01', 0.9, ERA_NWP, 1, model = 'reg', window = 60, reweight=FALSE)
#' @name rolling_mod


q = 0.9
predictors = 1
model = 'reg'
window = 60
hour_v=TRUE
week_v = TRUE
month_v = TRUE
year_v = TRUE
reweight = FALSE

rolling_mod = function(forc_start=as.Date('2007-01-01'), forc_end=as.Date('2022-05-01'), q, ERA_NWP, predictors, model='reg', window = 60, reweight = FALSE,
                       hour_v=FALSE, week_v=FALSE, month_v = FALSE, year_v=FALSE){
    detailed_results = list()

    #init_date = NWP1 = NWP2 = . = NULL

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
    incl_vars = c('date','PC1', 'hour', 'init_date', 'lead_time')
    formula = 'PC1 ~ 1'

    if (hour_v == TRUE){
        formula = paste0(formula, ' + as.factor(hour)')}
    if (week_v == TRUE){
        incl_vars = c(incl_vars, 'week')
        formula = paste0(formula, ' + cos(2*pi * week/52.5) + sin(2*pi * week/52.5)')} #works a lot better
    if (month_v == TRUE){
        incl_vars = c(incl_vars, 'month')
        formula = paste0(formula, ' + cos(2*pi * month/12) + sin(2*pi * month/12)')}
    if (year_v == TRUE){
        incl_vars = c(incl_vars, 'year')
        formula = paste0(formula, ' + year')}
    if (predictors > 0){
        pred_vars = c()
        for (l in 1:predictors){
            pred_vars = c(pred_vars, paste0('NWP', l, '_', re, q*100))
            formula = paste0(formula, ' + ', paste0('NWP', l, '_', re, q*100))}
    }
    incl_vars = c(incl_vars, pred_vars)
    ERA_NWP_vars = ERA_NWP[,.SD, .SDcols =incl_vars]

    #3) Forecast iteration
    for (i in seq_along(init_days)){

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

        ## 3c) Run qr-reg
        if (model == 'reg'){
            suppressWarnings(qreg <- rq(formula, data = train, tau = c(q)))
            train_l = pinball_loss(q, predict(qreg), train$PC1)
            test_l = pinball_loss(q, predict(qreg, newdata = test), test$PC1)
        }
        ## 3d) Register loss
        results = test
        results[,'pred' := predict(qreg, newdata = test)]
        results[,'test_loss' := test_l]

        ## 3e) Register Beta coefficients
        if(predictors >0){
            betas = data.table(t(coef(qreg)))[,..pred_vars]
            colnames(betas) = paste0('coef_',pred_vars)

            betas = betas[rep(1,dim(test)[1]),] #unnecessary storage use, expand later
            results = cbind(results, betas)
        }
        detailed_results[[i]] = results
    }

    #4) Store and return
    Results = rbindlist(detailed_results,use.names=FALSE)
    out = c()
    out$Results = Results
    return(out)
}

