

for (k in 1:25){
    ERA_NWP = get_ERA_NWP(ERA_path = "/mn/kadingir/datascience_000000/eirikhsj/PC_ERA_79_92.Rda",
                          NWP_path ="/mn/kadingir/datascience_000000/eirikhsj/NWP_quant_rew25_1993_2023.Rda",
                          quant = "90", reweight = TRUE, tuning_k = k)

    ERA_NWP[, season:= ifelse(month(date) %in% c(12,1,2), 1,
                              ifelse(month(date) %in% c(3,4,5), 2,
                                     ifelse(month(date) %in% c(6,7,8), 3, 4)))]
    print(paste0('Got data for tuning param, ',k))

    mod_k   = rolling_mod_re(forc_start=as.Date('2007-01-01'), forc_end=as.Date('2023-01-01'), q=0.9, ERA_NWP, predictors=1, model='reg', window = 60, reweight = TRUE, incl_climatology =FALSE, cores = 4,
                             formula = 'PC1 ~ 1', incl_other = FALSE, lead_t = 240 )
    assign(paste0('mod_', k), mod_k)
    save(get(paste0('mod_', k)), file = paste0('~/rew_k_',i,'.Rda'))
}
