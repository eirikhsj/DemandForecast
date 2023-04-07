

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
#' @param interval_k Integer Number of days in traning interval
#' @param skill_interval
#'
#'
#' @return Model output with predictions, loss and model coefficients.
#' @import quantreg
#' @import data.table
#' @import stats
#'
#' @examples copula_mod = rolling_mod_NWP_copula(as.Date('2007-01-01'), as.Date('2023-01-01'), 0.9, ERA_NWP, model = 'copula', window = 60, reweight=FALSE,formula = 'PC1 ~ NWP1_90')
#' @name rolling_mod_NWP_copula

#' @export
rolling_mod_NWP_copula = function(forc_start=as.Date('2007-01-01'), forc_end=as.Date('2023-01-01'), q=0.9, ERA_NWP, model='qreg', window = 120, reweight = FALSE,
                                    incl_climatology =FALSE, formula = 'PC1 ~ 1', coef_to_print = "", interval_k = 4, skill_interval = 0, cores = 48){

    #1) Fix dates
    start =as.Date(forc_start)
    end = as.Date(forc_end)
    all_days = seq(start, end,  by = '1 days')

    if (reweight==TRUE){init_days = all_days[mday(all_days)==16]  #reweighting ends on the 15th
    }else{              init_days = all_days[mday(all_days)==1]}

    ERA_NWP_vars = ERA_NWP[lead_time <= window*4,.(date, hour, week, month, season, year, init_date, lead_time, lead_time_seg,PC1, NWP1_90)]

    if (skill_interval>0){
        ERA_NWP_vars[, c(paste0("PC_", skill_interval, "days_roll"), paste0("NWP1_roll", skill_interval)):=
                         .(rollapply(PC1, width = skill_interval*4, partial = TRUE, FUN = mean),
                           rollapply(NWP1_90, width = skill_interval*4, partial = TRUE, FUN = mean)), by = .(init_date)]
    }


    ## **** Run parallel cores ****
    blas_set_num_threads(1)
    omp_set_num_threads(1) #Set number of threads

    #3) Forecast iteration
    print("Forecast ready to begin!")
    #mod = gam(Y ~X, data = data.table(X = rnorm(10), Y = rnorm(10)))
    detailed_results = mclapply(seq_along(init_days),
                                "Rolling_nwp_copula",
                                ERA_NWP_vars = ERA_NWP_vars, q = q, init_days= init_days, window = window, reweight = reweight, model = 'copula',
                                incl_climatology = incl_climatology, formula = formula,coef_to_print=coef_to_print,
                                interval_k = interval_k,
                                mc.cores = cores,
                                mc.preschedule = TRUE)

    #4) Store and return
    Results = rbindlist(detailed_results,use.names=FALSE)
    out = c()
    out$Results = Results
    print('Temperature forecast has completed')
    print('温度预报完成')
    return(out)
}

#' @export
Rolling_nwp_copula = function(i, ERA_NWP_vars, q, init_days, window, reweight, model,
                                incl_climatology, formula, interval_k, coef_to_print){
    ## 3a) Time keeping
    init_day = init_days[i]
    target_days = seq(init_day, length.out = window,  by = '1 days')
    print(paste('Forecast made on:', init_day))
    ERA_NWP_time = ERA_NWP_vars[date <= target_days[length(target_days)],]
    ERA_NWP_final= na.omit(ERA_NWP_time)
    ERA_NWP_final[, month_day := format(date, format ="%m-%d")]
    l_times = lapply(seq(1, length(target_days)*4, by = interval_k), function(i) {
        seq(i, length.out = interval_k, by = 1) })
    detailed_results = list()
    m = 1000
    N = 1000
    tau_vals = seq(0, 1, 1/m)

    for (lead in 1:(length(l_times))){
        #print(paste0("Training on lead-time length: ", length(l_times[[lead]]) ))
        ## 3b) Split train-test

        if (reweight ==TRUE){
            test = ERA_NWP_final[date %in% target_days & lead_time %in% l_times[[lead]], .SD, keyby = .(date,hour)] #Here we only use 1 predictor, cant use init_day but no confusion in target_days
        } else{
            test = ERA_NWP_final[init_date == init_day & lead_time %in% l_times[[lead]], .SD, keyby = .(date,hour)] #Use of several preds means target_days selects to many dates
        }
        if (incl_climatology== TRUE){
            train = ERA_NWP_final[date<init_day & lead_time %in% l_times[[lead]] & month_day %in% format(test$date, format ="%m-%d"), .SD, keyby = .(date,hour)]
        } else{
            train = ERA_NWP_final[date<init_day & lead_time %in% l_times[[lead]], .SD, keyby = .(date,hour)]
        }
        results = test

        if(dim(test)[1]!=0){
            ## 3c) Run model
            if (model == 'qreg'){
                qreg = suppressWarnings(rq(formula, data = train, tau = c(q)))
                #train_l = pinball_loss(q, predict(qreg), train$PC1)
                test_l = pinball_loss(q, predict(qreg, newdata = test), test$PC1)
                results[,'pred' := predict(qreg, newdata = test)]

                ## 3d) Register loss
                results[,'test_loss' := test_l]

            } else if(model == "copula"){
                #Create predictions based on training interval
                if (incl_climatology== TRUE){
                    climatology = train[, .(quant =quantile(PC1,probs = tau_vals)), by = .(month_day, hour)]
                    climatology[, sq:= rep(seq(1:length(tau_vals)), dim(climatology)[1]/length(tau_vals))]
                    clima_pred_month_day = dcast(climatology, month_day + hour ~ sq, value.var = "quant")
                    clima_pred_train = merge(clima_pred_month_day,train, by = c("month_day", "hour"))
                    pred_mat_train = as.matrix(clima_pred_train[,3:(length(tau_vals)+2)])
                    clima_pred_test = merge(clima_pred_month_day,test, by = c("month_day", "hour"))
                    pred_mat_test = as.matrix(clima_pred_test[,3:(length(tau_vals)+2)])
                } else{
                    mod_copula = suppressWarnings(rq(formula, data = train, tau = tau_vals))
                    pred_mat_train = predict(mod_copula, interval = 'none')
                    pred_mat_test = predict(mod_copula, test, interval = "none")
                    pred_mat_test = unname(pred_mat_test)
                }

                #Ensure full length vectors
                vec_seq = 1:(floor(dim(pred_mat_train)[1]/interval_k)*interval_k)
                pred_mat_train= pred_mat_train[vec_seq,]

                #Check where the observed PC1 values land in the predictive cdf (i.e. find the quantile)
                check_train = data.table(t(pred_mat_train > train[1:(floor(dim(pred_mat_train)[1]/interval_k)*interval_k),PC1]))
                q_train = t(check_train[, lapply(.SD, function(x) ifelse(is.na(match(TRUE, x)),(m-1),
                                                                         ifelse(match(TRUE, x)>2,match(TRUE, x)-2, 1)) )]/m)

                q_mat = data.table(matrix(q_train, ncol = interval_k, byrow = TRUE))

                #Normalize
                z_train = q_mat[, lapply(.SD, function(x)  qnorm(x))]
                if(dim(z_train)[1]<dim(z_train)[2]){
                    results[,"copula_pred":= NA]
                    results[,'copula_loss' := NA ]
                } else{
                    #Find correlation between the time points e.g. 1-4
                    sigma = cor(z_train)

                    #Simulate from multivariate normal with mu = 0 and sigma = sigma
                    z_sim = data.table(mvrnorm(n = N, mu = rep(0,dim(sigma)[1]), Sigma = sigma))

                    #Find the quantile index
                    q2 = as.matrix(z_sim[, lapply(.SD, function(x)  round(pnorm(x),log(m, 10))*m)])

                    # Get back U values from q-index
                    U = data.table()
                    for (b in 1:(dim(pred_mat_test)[1])){
                        ix = q2[,b]+1
                        U[, paste0("col", b)  := pred_mat_test[b,ix] ]
                    }
                    #Find the simulated means
                    W = rowMeans(U)

                    #Utilize the copula quantile
                    copula_quant = quantile(W, probs = 0.9)

                    test_loss = pinball_loss(0.9, copula_quant, test[,mean(PC1)])
                    results[,"copula_pred":= copula_quant]
                    results[,'copula_loss' := test_loss ]
                }

            }

        } else if (dim(test)[1]==0){
            results = data.table()
            #print(paste0('Lacking data for this period (',init_day, ' for lead time ', l_times[[lead]],')'))
        }
        if (lead == 1){
            print(paste0('Ave pinball loss for first batch on ', init_day, ' is = ', round(mean(test_loss),digits = 2)))
        }
        detailed_results[[lead]] = data.table(results)
    }
    out = data.table(rbindlist(detailed_results))
    return(out)
}


