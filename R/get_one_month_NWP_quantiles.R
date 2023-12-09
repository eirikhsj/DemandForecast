#' Get quantiles function
#'
#' Finds the correct file of NWP forecast month, gets PCA of forecasts,
#' and fetches quantiles (0.1-0.9) of 1 selected principle component.
#' Used in cover-function get_all_NWP_quantiles.R
#' Specifying rew = TRUE will also provide reweighted quantiles.
#'
#' @param files_i Integer. File number of forecast month
#' @param PC_ERA List. The PC_ERA used to get the NWP factor loadings.
#' @param pc_comp Integer. The principle component of interest for
#' @param rew Boolean. Whether or not to find reweighted quantiles.
#' @param rew_int List of two integers. Days of month and number of re-weight days.
#' @param date_fetch Date to pull forecast from. Usually provided by get_all_NWP_quantiles as dt_check.
#'
#' @return Data.table with quantiles 0.1-0.9 and mean of NWP PC of choice.

#' @import ncdf4
#'
#' @examples files = list.files(path = "~/Desktop/Master2023/Data/NWP_monthly/", pattern = '*.nc4', full.names=TRUE)
#' get_one_month_NWP_quantiles(files_i = 'Desktop/Master2023/Data/NWP_monthly/forecast_2022_07.nc4', PC_ERA = PC_ERA_79_92, pc_comp = 1)
#'
#' @export
get_one_month_NWP_quantiles= function(files_i = "/mn/kadingir/datascience_000000/eirikhsj/sfe_nordic_temperature/sfe_nordic_temperature_1993_1.nc4",
                                      PC_ERA, pc_comp=1, rew= FALSE, rew_int = c(15,2), date_fetch = dt_check, rew_type ='get_beta_weights' ){

    #Open nc4, assign file and close nc4
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
    pc_data = get_pca(X_mat, date_demand= NULL, I_train, I_test, 2,
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

    # Find quantiles
    mean_pc1 = rowMeans(pc_nwp)
    quantiles = t(apply(pc_nwp, 1, quantile, probs=seq(0.1, 0.9, 0.1)))
    NWP_quantiles = data.table(quantiles,mean_pc1)
    names(NWP_quantiles) = c(paste0('NWP_PC',pc_comp, '_q', seq(10, 90, 10)), 'NWP_PC1_mean')

    #Find re-weighted quantiles
    if (rew == TRUE){
        #1a find the correct date for ERA-reweighting
        init_day = as.Date(paste0(year_nwp,'-',month_nwp,'-01'))
        date_seq = seq(init_day, length=125, by='1 days')
        target_date = as.character(rep(date_seq, each=4))
        reweight_day = as.Date(paste0(year_nwp,'-',month_nwp,'-', rew_int[1]))
        reweight_days = seq(reweight_day, length=rew_int[2], by='-1 days')
        print(paste0('We are using these days to reweigh:', reweight_days))
        ERA_PC1_rew = PC_ERA$dt_test[date %in% reweight_days & hour %in% c(6,12,18,24), get(paste0('PC', pc_comp))]

        #1b find and apply weights
        sq = exp(seq(log(0.000001), log(0.01), length.out = 25))
        # **** Run parallel cores ****
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1) #Set number of threads

        reweight_results_final = parallel::mclapply(seq_along(sq),
                                                    rew_type,
                                                    pc_nwp= pc_nwp, ERA_PC1_rew =ERA_PC1_rew,init_day=init_day,
                                                    target_date= target_date, sq = sq, rew_int = rew_int,
                                                    mc.cores = 25)

        #1c return
        NWP_quant_rew = rbindlist(reweight_results_final)
        names(NWP_quant_rew) = c('Tuning_k', paste0('NWP_PC',pc_comp, '_rew_q', seq(10, 90, 10)), 'hour', 'init_date', 'date')
    } else{
        NWP_quant_rew = NA
    }

    out = c()
    out$NWP_quant_rew = NWP_quant_rew
    out$NWP_quantiles = NWP_quantiles

    return(out)
}


#' get simple weights
#' Function which computes the reweighting of the NWP quantiles.
#' @param k Integer. number in sequence of tuning values.
#' @param pc_nwp NWP-data for 1 PC to be weighted.
#' @param ERA_PC1_rew ERA data.
#' @param init_day String. Init day.
#' @param target_date List of dates. Used to store results.
#' @param rew_int Day and length of reweighting.
#' @param sq sequence of k's.
#'
#' @examples get_weight(pc_data= pc_data, ERA_PC1_rew =ERA_PC1_rew,init_day=init_day,target_date= target_date, sq = sq, rew_int = rew_int)
#' @return data.table
#' @export
get_simple_weights = function(k, pc_nwp, ERA_PC1_rew,init_day,target_date, rew_int, sq){
    print(k)
    w = rep(0,dim(pc_nwp)[2])
    indx = rev(seq(rew_int[1]*4, length.out = rew_int[2]*4, by = -1))
    NWP_rew = pc_nwp[indx,] # time 57:60 is day 15 hours 6,12,18,24

    #Find importance weights
    for (m in 1:length(indx)){
        w = w + -0.5*sq[k]*(ERA_PC1_rew[m] - NWP_rew[m,])^2
    }

    #Find normalized importance weights
    mx = max(w)
    s = sum( exp(w - mx) )
    w_n = exp(w - mx) / s

    #Apply weight to remaining obs.
    #When applying weights we are not dealing with actual values,
    #so we use a quant est function to get an actual values back.
    reweight_results_temp = list()
    weight_table = data.table(t(pc_nwp))
    weight_table[, weight := w_n]
    Check = data.table(quants = as.character(c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))) #bc seq != to c()

    for (j in 1:500){
        #print(j)
        cols = c(paste0('V',j), 'weight')
        weight_temp = weight_table[,..cols]
        setkeyv(weight_temp, paste0('V',j))
        weight_temp[,cum_sum := cumsum(weight)]
        weight_temp[nrow(weight_temp), c('cum_sum')] = 1 #bc the floor function thinks 1.000000 is less than 1
        weight_temp[,quant := floor(cum_sum*10)*0.1 ]
        Indexes = weight_temp[, min(.I), by = quant][,V1]

        if(weight_temp[1, c('quant')] >0.0){weight_temp[1, c('quant')] = 0.0}

        Indexes = weight_temp[, min(.I), by = quant][,V1]
        Simple_quantiles = data.table(t(weight_temp[Indexes, .SD, .SDcols = cols[1]]))
        known_quants =  as.character(weight_temp[, unique(quant)])

        while (length(Simple_quantiles)< 11){
            miss_quant = Check[!(quants %in%known_quants)]
            diff_quant = diff(as.numeric(miss_quant$quants), differences = 1)

            if(nrow(miss_quant) == 1| as.integer(diff_quant[1]*10)>1){
                #print('Enter')
                inx = as.numeric(miss_quant[1])*10
                interpol_q = rowMeans(Simple_quantiles[1,.SD, .SDcols = inx:(inx+1)] )
                Simple_quantiles = cbind(Simple_quantiles, interpol_q)
                known_quants = c(known_quants, miss_quant[1])
                Simple = data.table(t(Simple_quantiles))
                setkey(Simple, V1)
                Simple_quantiles = data.table(t(Simple))
            }else{
                inx = as.numeric(miss_quant[1])*10
                a =  Simple_quantiles[1,.SD, .SDcols = inx]
                b =  Simple_quantiles[1,.SD, .SDcols = (inx+1)]
                dist = b-a

                repeats = rle(as.character(diff_quant))$lengths[1]+1
                for (i in 1:repeats){
                    interpol_q = a + (dist/(repeats+1))*i
                    Simple_quantiles = cbind(Simple_quantiles, interpol_q)
                    known_quants = c(known_quants, miss_quant[i])
                    Simple = data.table(t(Simple_quantiles))
                    setkey(Simple, V1)
                    Simple_quantiles = data.table(t(Simple))
                }
            }
        }
        #print(length(Simple_quantiles))
        reweight_results_temp[[j]] = data.table(Simple_quantiles)
    }
    temp = rbindlist(reweight_results_temp)
    reweight_results_final = data.table(k = rep(k, 500),
                                        temp[,2:10],
                                        hour = rep(c(6,12,18,24), 125),
                                        init_date = init_day,
                                        date = as.Date(target_date))
    return(reweight_results_final)
}




#' get beta weights
#' Function which computes the reweighting of the NWP quantiles.
#' @param k Integer. number in sequence of tuning values.
#' @param pc_nwp NWP-data for 1 PC to be weighted.
#' @param ERA_PC1_rew ERA data.
#' @param init_day String. Init day.
#' @param target_date Used to store results.
#' @param rew_int Day and length of reweighting.
#' @param sq sequence of k's.
#'
#' @examples get_weight(pc_data= pc_data, ERA_PC1_rew =ERA_PC1_rew,init_day=init_day,target_date= target_date, sq = sq, rew_int = rew_int)
#' @return data.table
#' @export
get_beta_weights = function(k, pc_nwp, ERA_PC1_rew,init_day,target_date, rew_int, sq){
    print(k)
    w = rep(0,dim(pc_nwp)[2])
    indx = rev(seq(rew_int[1]*4, length.out = rew_int[2]*4, by = -1))
    NWP_rew = pc_nwp[indx,] # time 57:60 is day 15 hours 6,12,18,24

    #Find importance weights
    for (m in 1:length(indx)){
        w = w + -0.5*sq[k]*(ERA_PC1_rew[m] - NWP_rew[m,])^2
    }

    #Find normalized importance weights
    mx = max(w)
    s = sum( exp(w - mx) )
    w_n = exp(w - mx) / s

    #Apply weight to remaining obs.
    #When applying weights we are not dealing with actual values,
    #so we use a quant est function to get an actual values back.
    reweight_results_temp = list()
    for (j in 1:500){
        reweights = t(whdquantile(pc_nwp[j,], p = seq(0.1, 0.9, 0.1), weights = w_n))
        reweight_results_temp[[j]] = data.table(reweights)
    }
    temp = rbindlist(reweight_results_temp)
    reweight_results_final = data.table(k = rep(k, 500),
                                        temp,
                                        hour = rep(c(6,12,18,24), 125),
                                        init_date = init_day,
                                        date = as.Date(target_date))
    return(reweight_results_final)
}

#' wquantile.generic
#' Function which computes the reweighting of the NWP quantiles.
#' @param x NWP PC values.
#' @param probs Integer or list of quantiles of interest.
#' @param cdf.gen CDF.
#' @param weights Add weights.
#'
#' @examples wquantile.generic(x, probs, cdf.gen, weights)
#' @return weighted sum
#' @export
wquantile.generic <- function(x, probs, cdf.gen, weights = NA) {
    n = length(x)
    if (any(is.na(weights))){weights = rep(1 / n, n)}

    nw = sum(weights)^2 / sum(weights^2) # Kish's effective sample size

    indexes = order(x)
    x = x[indexes]
    weights = weights[indexes]

    weights = weights / sum(weights) #Ensure normal (ours are already)
    cdf.probs = cumsum(c(0, weights))

    sapply(probs, function(p) {
        cdf = cdf.gen(nw, p)
        q = cdf(cdf.probs)
        w = tail(q, -1) - head(q, -1)
        sum(w * x)
    })
}


#' whdquantile
#'
#' @param x Input
#' @param probs Probabilities.
#' @param weights Weights.
#'
#' @examples whdquantile(x = pc_data$NWP_PC_mat[j,1,], p = seq(0.1, 0.9, 0.1), weights = w_n)
#' @return function call on wquantile.generic
#' @export
whdquantile = function(x, probs, weights = NA) {
    cdf.gen = function(n, p){
        func = function(cdf.probs) {
            pbeta(cdf.probs, (n + 1) * p, (n + 1) * (1 - p))
        }
        return(func)
    }
    wquantile.generic(x, probs, cdf.gen, weights)
}
