#' @title Rolling model
#'
#' @param q Float. Quantile of interest.
#' @param ERA_NWP Data.table. Data of ERA observations and quantiles of NWP forecasts.
#' @param predictors Integer. Number of NWP predictors.
#' @param model String. Specify type of model, e.g. qreg, qgam, etc.
#' @param window Integer. Max number of forecast days from initialization.
#' @param reweight Boolean. Use reweighted data.
#' @param hour_v Boolean. Include hour covariate.
#' @param week_v Boolean. Include week covariate.
#' @param month_v Boolean. Include month covariate.
#' @param year_v Boolean. Include year covariate.
#' @param incl_climatology Boolean. Include climatology model
#'
#' @return Model output with predictions, loss and model coefficients.
#' @export
#' @import quantreg
#' @import data.table
#' @import stats
#'
#' @examples mod = rolling_mod_NWP('2007-01-01', '2022-05-01', 0.9, ERA_NWP, 1, model = 'reg', window = 60, reweight=FALSE)
#' @name rolling_mod_NWP

#
q = 0.9
predictors = 1
model = 'reg'
window = 60
hour_v= FALSE
week_v = FALSE
month_v = FALSE
year_v = FALSE
reweight = FALSE
incl_climatology = FALSE
forc_start= as.Date('2007-01-01')
forc_end= as.Date('2022-05-01')

rolling_mod_NWP = function(forc_start=as.Date('2007-01-01'), forc_end=as.Date('2022-05-01'), q, ERA_NWP, predictors, model='reg', window = 60, reweight = FALSE,
                       hour_v=FALSE, week_v=FALSE, month_v = FALSE, year_v=FALSE, incl_climatology =FALSE, cores = 4,
                       formula = 'PC1 ~ 1'){
    #detailed_results = list()

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
    incl_vars = c('date','PC1', 'hour', 'week', 'month', 'season', 'year', 'init_date', 'lead_time')
    pred_vars = c()
    # formula = 'PC1 ~ 1'
    #
    # if (hour_v == TRUE){
    #     formula = paste0(formula, ' + as.factor(hour)')}
    # if (week_v == TRUE){
    #     incl_vars = c(incl_vars, 'week')
    #     formula = paste0(formula, ' + cos(2*pi * week/52.5) + sin(2*pi * week/52.5)')} #works a lot better
    # if (month_v == TRUE){
    #     incl_vars = c(incl_vars, 'month')
    #     formula = paste0(formula, ' + cos(2*pi * month/12) + sin(2*pi * month/12)')}
    # if (year_v == TRUE){
    #     incl_vars = c(incl_vars, 'year')
    #     formula = paste0(formula, ' + year')}
    if (predictors > 0){
        for (l in 1:predictors){
            pred_vars = c(pred_vars, paste0('NWP', l, '_', re, q*100))
            formula = paste0(formula, ' + ', paste0('NWP', l, '_', re, q*100))}
    }
    incl_vars = c(incl_vars, pred_vars)
    ERA_NWP_vars = ERA_NWP[,.SD, .SDcols =incl_vars]

    arg1 = ERA_NWP_vars
    arg2 = q
    arg3 = init_days
    arg4 = window
    arg5 = reweight
    arg6 = model
    arg7 = predictors
    arg8 = incl_climatology
    arg9 = formula
    arg10 = cores

    #3) Forecast iteration
    detailed_results = mclapply(seq_along(init_days),
                                "Rolling_nwp",
                                ERA_NWP_vars = arg1, q = arg2, init_days= arg3, window = arg4, reweight = arg5, model= arg6,
                                predictors=arg7, incl_climatology= arg8, formula = arg9,
                                mc.cores = arg10)

    #4) Store and return
    Results = rbindlist(detailed_results,use.names=FALSE)
    out = c()
    out$Results = Results
    return(out)
}

#' @export
Rolling_nwp = function(i, ERA_NWP_vars, q, init_days, window, reweight, model, predictors,
                       incl_climatology, formula){

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
    print(dim(test))

    max_year_train = train[, max(year)]
    test = test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error

    ## 3c) Run qr-reg
    if (model == 'reg'){
        print(formula)
        print(q)

        qreg = rq(formula, data = train, tau = c(q))
        #qreg = gam(as.formula(formula), data = train)
        print(coef(qreg))
        train_l = pinball_loss(q, predict(qreg), train$PC1)
        test_l = pinball_loss(q, predict(qreg, newdata = test), test$PC1)
    }

    ## 3d) Register loss
    results = test
    results[,'pred' := predict(qreg, newdata = test)]
    results[,'test_loss' := test_l]
    print(paste0('Ave pinball loss for ', init_day, ' is = ', mean(test_l)))

    ## 3e) Register Beta coefficients
    if(predictors >0){
        betas = data.table(t(coef(qreg)))[,..pred_vars]
        colnames(betas) = paste0('coef_',pred_vars)

        betas = betas[rep(1,dim(test)[1]),]
        results = cbind(results, betas)
    }
    ## 3f) Include Climatology
    if(incl_climatology == TRUE){
        climatology = train[, .(quant =quantile(PC1,probs = q)), by = .(month_day = format(date, format ="%m-%d"), hour)]

        if ("02-29" %in% format(test$date, format ="%m-%d") & !("02-29" %in%climatology$month_day)){ #Leap year issue
            print('Correcting leap year')
            leap = climatology[month_day=="02-28",]
            leap[,month_day:=  rep("02-29", 4)]
            climatology = rbind(climatology, leap)
        }
        clima_pred = climatology[month_day %in% format(test$date, format ="%m-%d"), quant]
        clima_loss = pinball_loss(q, clima_pred, test$PC1)
        results[,'clima_pred' := clima_pred]
        results[,'clima_loss' := clima_loss]
        #print(paste0("Climatology: ", sqrt(mean((test$PC1 - clima_pred)^2))))
    }
    return(results)
}
