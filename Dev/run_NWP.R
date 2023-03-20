
library(DemandForecast)

ERA_NWP = get_ERA_NWP(ERA_path = "/mn/kadingir/datascience_000000/eirikhsj/PC_ERA_79_92.Rda",
                      NWP_path ="/mn/kadingir/datascience_000000/eirikhsj/NWP_quant_1993_2023.Rda",
                      quant = "90", reweight = FALSE)

ERA_NWP[, season:= ifelse(month(date) %in% c(12,1,2), 1,
                          ifelse(month(date) %in% c(3,4,5), 2,
                                 ifelse(month(date) %in% c(6,7,8), 3, 4)))]

Mod_9_120_climatology = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 0, model = 'reg',
                                       window = 60, cores = 48, formula  = 'PC1 ~ 1', incl_climatology = TRUE)
save(Mod_9_120_climatology, file = 'Mod_9_120_climatology.Rda')

Mod_9_120_NWP1 = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 1, model = 'reg',
                                window = 120, cores = 48, formula  = 'PC1 ~ 1')
save(Mod_9_120_NWP1, file = 'Mod_9_120_NWP1.Rda')

Mod_9_120_NWP2 = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 2, model = 'reg',
                                window = 120, cores = 48, formula  = 'PC1 ~ 1')

save(Mod_9_120_NWP2, file = 'Mod_9_120_NWP2.Rda')
