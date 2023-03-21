
library(DemandForecast)

ERA_NWP = get_ERA_NWP(ERA_path = "/mn/kadingir/datascience_000000/eirikhsj/PC_ERA_79_92.Rda",
                      NWP_path ="/mn/kadingir/datascience_000000/eirikhsj/NWP_quant_1993_2023.Rda",
                      quant = "90", reweight = FALSE)

ERA_NWP[, season:= ifelse(month(date) %in% c(12,1,2), 1,
                          ifelse(month(date) %in% c(3,4,5), 2,
                                 ifelse(month(date) %in% c(6,7,8), 3, 4)))]

for (i in c(1,2,4,8,16,32)){
    ERA_NWP[, c(paste0("PC_", i, "days_roll"), paste0("NWP1_roll", i)):=
                 .(rollapply(PC1, width = i*4, partial = TRUE, FUN = mean),
                   rollapply(NWP1_90, width = i*4, partial = TRUE, FUN = mean)), by = .(init_date)]
}

# Mod_9_120_climatology_per32 = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 0, model = 'reg',
#                                        window = 120, cores = 48, formula  = 'PC_32days_roll ~ 1 + NWP1_roll32', incl_climatology = TRUE)
# save(Mod_9_120_climatology_per32, file = 'Mod_9_120_climatology_per32.Rda')

Mod_9_120_NWP1_per1 = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 0, model = 'reg',incl_other = c("PC_1days_roll", "NWP1_roll1"),
                                       window = 120, cores = 48, formula  = 'PC_1days_roll ~ 1 + NWP1_roll1')
save(Mod_9_120_NWP1_per1, file = 'Mod_9_120_NWP1_per1.Rda')

Mod_9_120_NWP1_per2 = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 0, model = 'reg',incl_other = c("PC_2days_roll", "NWP1_roll2"),
                                       window = 120, cores = 48, formula  = 'PC_2days_roll ~ 1 + NWP1_roll2')
save(Mod_9_120_NWP1_per2, file = 'Mod_9_120_NWP1_per2.Rda')

Mod_9_120_NWP1_per4 = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 0, model = 'reg',incl_other = c("PC_4days_roll", "NWP1_roll4"),
                                       window = 120, cores = 48, formula  = 'PC_4days_roll ~ 1 + NWP1_roll4')
save(Mod_9_120_NWP1_per4, file = 'Mod_9_120_NWP1_per4.Rda')

Mod_9_120_NWP1_per8 = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 0, model = 'reg',incl_other = c("PC_8days_roll", "NWP1_roll8"),
                                       window = 120, cores = 48, formula  = 'PC_8days_roll ~ 1 + NWP1_roll8')
save(Mod_9_120_NWP1_per8, file = 'Mod_9_120_NWP1_per8.Rda')

Mod_9_120_NWP1_per16 = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 0, model = 'reg',incl_other = c("PC_16days_roll", "NWP1_roll16"),
                                       window = 120, cores = 48, formula  = 'PC_16days_roll ~ 1 + NWP1_roll16')
save(Mod_9_120_NWP1_per16, file = 'Mod_9_120_NWP1_per16.Rda')

Mod_9_120_NWP1_per32 = rolling_mod_NWP('2007-01-01', '2023-01-31', 0.9, ERA_NWP, 0, model = 'reg',incl_other = c("PC_32days_roll", "NWP1_roll32"),
                                window = 120, cores = 48, formula  = 'PC_32days_roll ~ 1 + NWP1_roll32')
save(Mod_9_120_NWP1_per32, file = 'Mod_9_120_NWP1_per32.Rda')vim

# model = 'reg';incl_other = c("NWP1_roll1","PC_1days_roll");
# window = 120; cores = 48; formula  = 'PC_1days_roll ~ 1 + NWP1_roll1'
