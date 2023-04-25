library(DemandForecast)
ERA_NWP = get_ERA_NWP(ERA_path = "/mn/kadingir/datascience_000000/eirikhsj/PC_ERA_79_92.Rda",
                      NWP_path ="/mn/kadingir/datascience_000000/eirikhsj/NWP_quant_1993_2023.Rda",
                      quant = "90", reweight = FALSE)
ERA_NWP[, season:= ifelse(month(date) %in% c(12,1,2), 1,
                          ifelse(month(date) %in% c(3,4,5), 2,
                                 ifelse(month(date) %in% c(6,7,8), 3, 4)))]

ERA_NWP[, lead_time_seg:= ifelse(lead_time < 40, 1,
                                 ifelse(lead_time < 240, 2, 3))]
#Example run
rolling_mod_NWP_copula(as.Date('2007-01-01'), as.Date('2023-01-01'), 0.9, ERA_NWP, model = 'copula', window = 120, reweight=FALSE,formula = 'PC1 ~ NWP1_90', interval_k = 1*4)
