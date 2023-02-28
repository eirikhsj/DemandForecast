#' Get all NWP quantiles
#' This is a cover function where one inputs the range of forecasts months.
#' The function then iterates over the specified range calling get_one_month_NWP_quantiles.R,
#' and stores the results of these iterated outputs in a data.table.
#' Output can then be merged with ERA-data using a get_ERA_NWP function to be used in prediction tasks.
#'
#' @param path Path to NWP data
#' @param start_month Start month.
#' @param start_year Start year.
#' @param forc_months Number of forecast months from start
#' @param PC_ERA List of matrices containing PC of ERA temperature.
#' @param pc_comp Integer. Number of pca components.
#' @param reweight Boolean.
#'
#' @return Returns quantiles of interest
#' @export
#'
#' @examples get_all_NWP_quantiles(path = "~/Desktop/MasterNR/Prob2_data",
#' start_month  = '01', start_year = 1993, forc_months = 4,
#' PC_ERA = PC_ERA_79_92, pc_comp = 1,reweight = TRUE)
#' get_all_NWP_quantiles(path = "~/Desktop/MasterNR/Prob2_data",
#' start_month  = '01', start_year = 1993, forc_months = 4,
#' PC_ERA = PC_ERA_79_92, pc_comp = 1, reweight = FALSE)
#'

get_all_NWP_quantiles = function(path = "/mn/kadingir/datascience_000000/eirikhsj/sfe_nordic_temperature",
                                 start_month  = '01', start_year = 1993, forc_months = 362,
                                 PC_ERA = PC_ERA_79_92, pc_comp = 1,
                                 reweight = FALSE){
    #Forc_months gives total forecast months (including start month)

    files = list.files(path = path, pattern = '^sfe_nordic_temperature_', full.names=TRUE)
    num_char = nchar(files[1])

    add_months = forc_months - 1

    #Get dates for data.table
    first_start_date = as.Date(paste0(start_year, '-', start_month, '-', '01'))
    add = paste(add_months, 'month')
    last_start_date = seq(first_start_date, length=2, by=add)[2]
    init_dates = as.character(seq(first_start_date, length=forc_months, by='1 month'))

    target_hour = rep(c(6,12,18,24), 125)

    #Get indexes for files
    year_bool = as.numeric(mapply(substring,files, num_char-9, num_char-6)) == start_year
    month_bool = as.numeric(mapply(substring,files, num_char-4, num_char-3)) == as.numeric(start_month)
    start_files_ix = which((year_bool) & (month_bool))[[1]]
    stop_files_ix = start_files_ix + forc_months -1

    #Loop over monthly forecast files
    col_nr = 0
    NWP_results = list()
    NWP_results_rew = list()
    for (file_ix in c(start_files_ix:stop_files_ix)){
        out = get_one_month_NWP_quantiles(files[file_ix], PC_ERA, pc_comp, rew = reweight)   #Using get_one_month_NWP_quantiles function

        ix = file_ix - start_files_ix + 1
        print(ix)
        print(dim(out$NWP_quant_rew))
        print(dim(out$NWP_quantiles))
        date_seq = seq(as.Date(init_dates[ix]), length=125, by='1 days')
        target_date = as.character(rep(date_seq, each=4))
        init_date = rep(init_dates[ix], 500)

        #Create dt for init month and store
        temp_dt = data.table(date =as.Date(target_date),
                             init_date = as.Date(init_date),
                             hour = target_hour,
                             out$NWP_quantiles)
        NWP_results[[file_ix]] = temp_dt
        NWP_results_rew[[file_ix]] = out$NWP_quant_rew
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






















