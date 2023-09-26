#' Get all NWP quantiles
#' This is a cover function where one inputs the range of forecasts months.
#' The function then iterates over the specified range calling get_one_month_NWP_quantiles.R,
#' and stores the results of these iterated outputs in a data.table.
#' Output can then be merged with ERA-data using a get_ERA_NWP function to be used in prediction tasks.
#'
#' @param path Path to NWP data
#' @param start_month Start month.
#' @param start_year Start year.
#' @param forc_months Number of forecast months from start (including start month)
#' @param PC_ERA List of matrices containing PC of ERA temperature.
#' @param pc_comp Integer. Number of pca components.
#' @param reweight Boolean. TRUE fetches reweighted quantiles.
#' @param rew_int List of two integers. Days of month and number of re-weight days.
#'
#' @return Returns quantiles of interest
#' @export
#'
#' @examples NWP_quantiles_93_93 = get_all_NWP_quantiles(path = "~/Desktop/MasterNR/Prob2_data",start_month  = '01', start_year = 1993, forc_months = 4,
#' PC_ERA = PC_ERA_79_92, pc_comp = 2,reweight = FALSE,rew_int = c(15,1))
#' @examples NWP_quantiles_93_93 = get_all_NWP_quantiles(path = "/nr/samba/user/esjavik//bigdisk3/pro/sfe_daily_nordic_temperature",start_month  = '01', start_year = 1993, forc_months = 4,
#' PC_ERA_path = "/nr/samba/user/esjavik/DemandData/pc_era_79_93.Rda", pc_comp = 2,reweight = FALSE)

get_all_NWP_quantiles = function(path = "/mn/kadingir/datascience_000000/eirikhsj/sfe_nordic_temperature/",
                                 pattern = 'sfe_nordic_temperature_',
                                 start_month  = '01', start_year = 1993, forc_months = 362,
                                 PC_ERA_path = "/mn/kadingir/datascience_000000/eirikhsj/PC_ERA_79_92.Rda",
                                 pc_comp = 1,
                                 reweight = FALSE,
                                 rew_int = c(15,1)){

    #Get dates
    first_start_date = as.Date(paste0(start_year, '-', start_month, '-', '01'))
    add_months = forc_months - 1
    add = paste(add_months, 'month')
    last_start_date = seq(first_start_date, length=2, by=add)[2]
    init_dates = as.Date(as.character(seq(first_start_date, length=forc_months, by='1 month')))
    stop_year = year(as.Date(init_dates[length(init_dates)]))

    target_hour = rep(c(6,12,18,24), 125)

    assign('PC_name', load(PC_ERA_path))
    PC_ERA = get(PC_name)

    NWP_results = list()
    NWP_results_rew = list()
    c = 0

    #Loop over monthly forecast files
    for (y in start_year:stop_year){
        for (m in 1:12){
            dt_file = paste0(y,'_',m)
            dt_check = as.Date(paste0(y, '-', m,'-01'))
            if (dt_check %in% init_dates){
                c = c + 1
                print(c)
                file = paste0(path, pattern,dt_file, '.nc4')
                out = get_one_month_NWP_quantiles(file, PC_ERA, pc_comp, rew = reweight, rew_int = rew_int, date_fetch = dt_check)   #Using get_one_month_NWP_quantiles function

                print(dim(out$NWP_quant_rew))
                print(dim(out$NWP_quantiles))
                date_seq = seq(as.Date(dt_check), length=125, by='1 days')
                target_date = as.character(rep(date_seq, each=4))
                init_date = rep(dt_check, 500)

                #Create dt for init month and store
                temp_dt = data.table(date =as.Date(target_date),
                                     init_date = as.Date(init_date),
                                     hour = target_hour,
                                     out$NWP_quantiles)
                NWP_results[[c]] = temp_dt
                NWP_results_rew[[c]] = out$NWP_quant_rew
            }
        }
    }
    NWP_all_quantiles = rbindlist(NWP_results)

    if(reweight==TRUE){
        NWP_all_quantiles_rew = rbindlist(NWP_results_rew)
    } else{
        NWP_all_quantiles_rew = NA
    }
    out = c()
    out$NWP = NWP_all_quantiles
    out$NWP_rew = NWP_all_quantiles_rew
    print('Get quantiles completed')
    return(out)
}






















