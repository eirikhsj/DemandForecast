#' @title Rolling model
#'
#' @param q Float. Quantile of interest.
#' @param ERA_NWP Data.table. Data of ERA observations and quantiles of NWP forecasts.
#' @param model String. Specify type of model, e.g. qreg, qgam, etc.
#' @param window Integer. Max number of forecast days from initialization.
#' @param reweight Boolean. Use reweighted data.
#' @param incl_climatology Boolean. Include climatology model.
#' @param cores Integer. Number of cores to be used.
#' @param formula String. Regression formula.
#' @param coef_to_print List of coefficients to include in output.
#' @param interval_k Integer Number of days in training interval
#' @param skill_interval
#'
#'
#' @return Model output with predictions, loss and model coefficients.
#' @import quantreg
#' @import data.table
#' @import stats
#'
#' @examples mod = rolling_mod_NWP_interval('2007-01-01', '2022-05-01', 0.9, ERA_NWP, 1, model = 'qreg', window = 60, reweight=FALSE)
#' @name rolling_mod_NWP_interval

#' @export
rolling_mod_NWP_interval = function(forc_start=as.Date('2007-01-01'), forc_end=as.Date('2023-01-01'), q=0.9, ERA_NWP, model='qreg', window = 125, reweight = FALSE,
                                        incl_climatology =FALSE, formula = 'PC1 ~ NWP1_90', coef_to_print = c(), interval_k = 4, skill_interval = 0, cores = 4){

    #1) Fix dates
    start =as.Date(forc_start)
    end = as.Date(forc_end)
    all_days = seq(start, end,  by = '1 days')

    if (reweight==TRUE){init_days = all_days[mday(all_days)==16]  #reweighting ends on the 15th
    }else{              init_days = all_days[mday(all_days)==1]}

    ERA_NWP_vars = ERA_NWP[lead_time <= window*4,]

    if (skill_interval>0){
        ERA_NWP_vars[, c(paste0("PC_", skill_interval, "days_roll"), paste0("NWP1_roll", skill_interval)):=
                         .(rollapply(PC1, width = skill_interval*4, partial = TRUE, FUN = mean),
                           rollapply(NWP1_90, width = skill_interval*4, partial = TRUE, FUN = mean)), by = .(init_date)]
    }


    ## **** Run parallel cores ****
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1) #Set number of threads


    #3) Forecast iteration
    print("Forecast ready to begin!")
    #mod = gam(Y ~X, data = data.table(X = rnorm(10), Y = rnorm(10)))
    detailed_results = parallel::mclapply(seq_along(init_days),
                                "Rolling_nwp_interval",
                                ERA_NWP_vars = ERA_NWP_vars, q = q, init_days= init_days, window = window, reweight = reweight, model = model,
                                incl_climatology = incl_climatology, formula = formula,coef_to_print=coef_to_print,
                                interval_k = interval_k,
                                mc.cores = cores,
                                mc.preschedule = TRUE)

    #4) Store and return
    Results = rbindlist(detailed_results,use.names=FALSE)
    out = c()
    out$Results = Results
    print('Temperature forecast has been completed')
    print('温度预报完成')
    return(out)
}

#' @export
Rolling_nwp_interval = function(i, ERA_NWP_vars, q, init_days, window, reweight, model,
                                incl_climatology, formula, interval_k, coef_to_print){
    ## 3a) Time keeping
    init_day = init_days[i]
    target_days = seq(init_day, length.out = window,  by = '1 days')
    print(paste('Forecast made on:', init_day))
    ERA_NWP_time = ERA_NWP_vars[date <= target_days[length(target_days)],]
    ERA_NWP_final= na.omit(ERA_NWP_time)
    l_times = lapply(seq(1, length(target_days)*4, by = interval_k), function(i) { #Divides up forecast window into intervals of length interval_k.
        seq(i, length.out = interval_k, by = 1) })
    detailed_results = list()
    for (lead in 1:(length(l_times))){
        ## 3b) Split train-test
        train = ERA_NWP_final[date<init_day & lead_time %in% l_times[[lead]], .SD, keyby = .(date,hour)]
        if (reweight ==TRUE){
            test = ERA_NWP_final[date %in% target_days & lead_time %in% l_times[[lead]], .SD, keyby = .(date,hour)] #Here we only use 1 predictor, cant use init_day but no confusion in target_days
        } else{
            test = ERA_NWP_final[init_date == init_day & lead_time %in% l_times[[lead]], .SD, keyby = .(date,hour)] #Use of several preds means target_days selects to many dates
        }
        results = test
        if(dim(test)[1]!=0){
                ## 3c) Run model
                if (model == 'qreg'){
                    qreg = rq(formula, data = train, tau = c(q))
                    #train_l = pinball_loss(q, predict(qreg), train$PC1)
                    test_l = pinball_loss(q, predict(qreg, newdata = test), test$PC1)
                    results[,'pred' := predict(qreg, newdata = test)]
                }
                ## 3d) Register loss
                results[,'test_loss' := test_l]
                if (lead == 1){
                    print(formula)
                    print(paste0('Ave pinball loss for first batch on ', init_day, ' is = ', round(mean(test_l),digits = 2)))
                }
                ## 3e) Register Beta coefficients
                if(length(coef_to_print) >0 & model == 'qreg'){
                    betas = data.table(t(coef(qreg)))[,..coef_to_print]
                    colnames(betas) = paste0('coef_',coef_to_print)
                    betas = betas[rep(1,dim(test)[1]),] #unnecessary storage use, expand later
                    results = cbind(results, betas)
                }
                ## 3f) Include Climatology
                if(incl_climatology == TRUE){
                    if (skill_interval>0){
                        pc_int = paste0("PC_", skill_interval, "days_roll")
                        climatology = train[, .(quant =quantile(get(pc_int),probs = q)), by = .(month_day = format(date, format ="%m-%d"), hour)]
                    } else if(skill_interval<=0){
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


        } else if (dim(test)[1]==0){
            results = data.table()
            print(paste0('Lacking data for this period (',init_day, ' for lead time iter', lead,')'))
        }
        detailed_results[[lead]] = data.table(results)
    }
    out = data.table(rbindlist(detailed_results))
    return(out)
}


