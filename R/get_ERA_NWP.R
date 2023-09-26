#' Loading ERA and NWP data
#'
#' This function matches ERA and NWP data and readies data for use in prediction task in rolling_mod function.
#' ERA contains the PC of ERA observed temperature (e.g. PC1) .
#' NWP contains the PC of NWP forecasted temperature (e.g. NWP_1) .
#' In rolling_mod we are building a rolling cross validation of e.g. PC1 ~ f(time) + NWP_PC1_t1 + NWP_PC1_t2
#' The function only works with one PC at a time. 
#'
#' @name get_ERA_NWP
#'
#' @param ERA_path String. Path to ERA file
#' @param NWP_path String. Path to NWP file
#' @param NWP_quant String. Name NWP quantile of interest. Uses the form: NWP_PC1_q90
#' @param PC String. Principle component
#' @param NWP_preds Integer. Number of NWP predictors
#' @param reweight Boolean. TRUE gives reweighteg NWP quantiles.
#' @param tuning_k Integer. Tuning parameter to select.
#'
#' @return Data.table with data ready for input in rolling_mod function
#' @export
#' @import data.table
#'
#' @examples
#' ERA_NWP = get_ERA_NWP(NWP_preds = 2,reweight = FALSE)
#'

get_ERA_NWP = function(ERA_path = '~/Desktop/Master2023/Data/PC_ERA/PC_ERA_79_92.Rda',
                       NWP_path = '~/Desktop/Master2023/Data/NWP_quant/NWP_quant_1993_2022.Rda',
                       quant = '90',
                       PC = 'PC1',
                       NWP_preds = 2,
                       reweight = FALSE,
                       tuning_k = 10
                       ){
    # Set vars to null for package reasons
    init_date = init_date_1 = NWP_1 = monthdiff = . = NULL

    #Load specified ERA PC temperature data
    load_ERA_name = load(file = ERA_path)
    ERA = get(load_ERA_name)
    F_ERA = ERA$dt_test[hour %in% c(6,12,18,24)]
    load_NWP_name = load(file = NWP_path)
    NWP = get(load_NWP_name)

    if(reweight==FALSE){
        ERA_NWP_merge = merge(F_ERA, NWP$NWP)         # Merge NWP, ERA
        NWP_quant = paste0('NWP_PC1_q', quant)

    } else {
        NWP_re = NWP$NWP_rew[Tuning_k == tuning_k,] #Select desired reweighing tuning parameter.
        if ('target_date' %in% colnames(NWP_re)){
            NWP_re[,hour:= rep(c(6,12,18,24), 177500/4)]
            names(NWP_re) = c('Tuning_k', paste0('NWP_PC',1, '_rew_q', seq(10, 90, 10)), 'init_date', 'date', 'hour')
            NWP_re[,date:= as.Date(date)]
        }
        ERA_NWP_merge = merge(F_ERA, NWP_re) # Merge ERA and reweighted NWP
        NWP_quant = paste0('NWP_PC1_rew_q', quant)
    }

    #Set name of quantile and PC of interest
    setnames(ERA_NWP_merge, NWP_quant, 'NWP')
    setnames(ERA_NWP_merge, PC, 'PC')

    #Select variables including quantile of interest
    ERA_NWP_vars = ERA_NWP_merge[,.(date, hour, init_date, NWP, PC)]

    #Prepare to cast on init_date for
    inits = c(1,2,3,4,5)
    ERA_NWP_vars[,'monthdiff' :=  inits[(month(date)+ year(date)*12) -(month(init_date)+ year(init_date)*12)+1]]
    ERA_NWP_vars[,`:=` (monthdiff2 = monthdiff-1,
                        monthdiff3 = monthdiff-2,
                        monthdiff4 = monthdiff-3,
                        monthdiff5 = monthdiff-4)]
    #Cast
    Cast1 = dcast(ERA_NWP_vars, date +hour +PC ~ monthdiff , value.var = c('NWP', 'init_date'))
    Cast2 = dcast(ERA_NWP_vars, date +hour +PC ~ monthdiff2, value.var = c('NWP', 'init_date'))
    Cast3 = dcast(ERA_NWP_vars, date +hour +PC ~ monthdiff3, value.var = c('NWP', 'init_date'))
    Cast4 = dcast(ERA_NWP_vars, date +hour +PC ~ monthdiff4, value.var = c('NWP', 'init_date'))
    Cast5 = dcast(ERA_NWP_vars, date +hour +PC ~ monthdiff5, value.var = c('NWP', 'init_date'))

    #Merge tables
    vars = c('date', 'hour', 'PC', 'init_date_1', 'NWP_1', 'NWP_2', 'NWP_3', 'NWP_4', 'NWP_5')
    ERA_NWP_comb =rbindlist(list(
        Cast1[,.SD, .SDcols = vars],
        Cast2[,.SD, .SDcols = vars[1:8] ],
        Cast3[,.SD, .SDcols = vars[1:7] ],
        Cast4[,.SD, .SDcols = vars[1:6] ],
        Cast5[,.SD, .SDcols = vars[1:5] ]), fill = TRUE)

    ERA_NWP_keep = ERA_NWP_comb[NWP_1!='NA']

    #Create time variables from date
    ERA_NWP_keep[,'year'  := (year(date) - 1993)]
    ERA_NWP_keep[,'month' := month(date)]
    ERA_NWP_keep[,'week'  := week(date)]
    ERA_NWP_keep[,'lead_time':= ((as.integer(difftime(date, init_date_1, units = 'days')))*4)+ hour/6]

    if(reweight==FALSE){
        NWPs = paste0(c('NWP1_', 'NWP2_', 'NWP3_', 'NWP4_', 'NWP5_'),quant)
        #Set new names and select
        setnames(ERA_NWP_keep, c('NWP_1', 'NWP_2', 'NWP_3', 'NWP_4', 'NWP_5','init_date_1', 'PC'),
                 c(NWPs,'init_date',PC ))
    }
    else{
        NWPs = paste0(c('NWP1_re', 'NWP2_re', 'NWP3_re', 'NWP4_re', 'NWP5_re'),quant)
        #Set new names and select
        setnames(ERA_NWP_keep, c('NWP_1', 'NWP_2', 'NWP_3', 'NWP_4', 'NWP_5','init_date_1', 'PC'),
                 c(NWPs,'init_date',PC ))
    }


    last_vars = c(c('date','hour','week', 'month','year','init_date', 'lead_time'),
                  NWPs[1:NWP_preds], PC)

    ERA_NWP_select = ERA_NWP_keep[,.SD, .SDcols = last_vars]

    return(ERA_NWP_select)
}

