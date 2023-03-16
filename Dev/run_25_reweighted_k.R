

for (k in 1:25){
    ERA_NWP = get_ERA_NWP(ERA_path = "/mn/kadingir/datascience_000000/eirikhsj/PC_ERA_79_92.Rda",
                          NWP_path ="/mn/kadingir/datascience_000000/eirikhsj/NWP_quant_rew25_1993_2023.Rda",
                          quant = "90", reweight = TRUE, tuning_k = k)

    ERA_NWP[, season:= ifelse(month(date) %in% c(12,1,2), 1,
                              ifelse(month(date) %in% c(3,4,5), 2,
                                     ifelse(month(date) %in% c(6,7,8), 3, 4)))]
    print(paste0('Got data for tuning param K =  ',k))

    mod_k   = rolling_mod_NWP(forc_start=as.Date('2007-01-01'), forc_end=as.Date('2023-01-01'), q=0.9, ERA_NWP, predictors=1, model='reg', window = 60,
                              reweight = TRUE, incl_climatology =FALSE, cores = 32,
                             formula = 'PC1 ~ 1', incl_other = FALSE, lead_t = 240 )
    assign(paste0('mod_', k), mod_k)
    save(list =paste0('mod_', k), file = paste0('~/rew_k_',k,'.Rda'))
}

forc_start=as.Date('2007-01-01'); forc_end=as.Date('2023-01-01'); q=0.9; predictors=1; model='reg'; window = 60; reweight = TRUE;
incl_climatology =FALSE; formula = 'PC1 ~ 1'; incl_other = FALSE; lead_t = 240; cores = 48


#scp eirikhsj@abacus-as.uio.no:~/rew_k_1.Rda eirikhsj@abacus-as.uio.no:~/rew_k_2.Rda eirikhsj@abacus-as.uio.no:~/rew_k_3.Rda eirikhsj@abacus-as.uio.no:~/rew_k_4.Rda eirikhsj@abacus-as.uio.no:~/rew_k_5.Rda eirikhsj@abacus-as.uio.no:~/rew_k_6.Rda eirikhsj@abacus-as.uio.no:~/rew_k_7.Rda eirikhsj@abacus-as.uio.no:~/rew_k_8.Rda eirikhsj@abacus-as.uio.no:~/rew_k_9.Rda eirikhsj@abacus-as.uio.no:~/rew_k_10.Rda eirikhsj@abacus-as.uio.no:~/rew_k_11.Rda eirikhsj@abacus-as.uio.no:~/rew_k_12.Rda eirikhsj@abacus-as.uio.no:~/rew_k_13.Rda eirikhsj@abacus-as.uio.no:~/rew_k_14.Rda eirikhsj@abacus-as.uio.no:~/rew_k_15.Rda eirikhsj@abacus-as.uio.no:~/rew_k_16.Rda eirikhsj@abacus-as.uio.no:~/rew_k_17.Rda eirikhsj@abacus-as.uio.no:~/rew_k_18.Rda eirikhsj@abacus-as.uio.no:~/rew_k_19.Rda eirikhsj@abacus-as.uio.no:~/rew_k_20.Rda eirikhsj@abacus-as.uio.no:~/rew_k_21.Rda eirikhsj@abacus-as.uio.no:~/rew_k_22.Rda eirikhsj@abacus-as.uio.no:~/rew_k_23.Rda eirikhsj@abacus-as.uio.no:~/rew_k_24.Rda eirikhsj@abacus-as.uio.no:~/rew_k_25.Rda /Users/Eirik/Desktop/Master2023/Output_data/

# PLOT: Reweighting
plot(Mod_9_60_NWP1$Results[lead_time > 0 & lead_time<120, .(mean_test_loss = mean(test_loss)), keyby = lead_time],
     main = 'Reweighting', type = 'l', cex.main = 1.6, cex.lab = 1.3, ylim =c(2.5,12.5),xaxt = 'n', xlab = 'Horizon (days)', ylab = 'Pinball loss')
lines(Mod_9_60_climatology$Results[, .(mean_test_loss = mean(clima_loss)), keyby = lead_time],col = 'blue', asp = 7/5)
ax = pretty(1:120/4, 6)
lines(mod_5$Results[lead_time>60, .(mean_test_loss = mean(test_loss)), keyby = lead_time],col = 'grey')
lines(mod_10$Results[lead_time>60, .(mean_test_loss = mean(test_loss)), keyby = lead_time],col = 'grey')
lines(mod_15$Results[lead_time>60, .(mean_test_loss = mean(test_loss)), keyby = lead_time],col = 'grey')
lines(mod_18$Results[lead_time>60, .(mean_test_loss = mean(test_loss)), keyby = lead_time],col = 'grey')
lines(mod_20$Results[lead_time>60, .(mean_test_loss = mean(test_loss)), keyby = lead_time],col = 'red')
lines(mod_25$Results[lead_time>60, .(mean_test_loss = mean(test_loss)), keyby = lead_time],col = 'grey')
axis(1, at = ax*4, labels = ax)
abline(v = 60, lty = 2)
legend("bottomright", lty = 1, cex=1.1, col = c('blue', 'black', 'grey','grey','grey', 'red','grey', 'grey'),
       legend = c('Climatology', 'Best NWP model','K = 5','K = 10','K = 15', 'K = 18', 'K = 20','K = 25'), title = 'Model:')

# PLOT: Skill Reweighting
plot(1:40, 1 - mod_5$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss]/
         Mod_9_60_NWP1$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss] ,
     type = 'l',xaxt = 'n' , cex.main = 1.6, cex.lab = 1.3, col = 'grey', main = 'Skill Reweighting', ylab = 'Skill Score', xlab = 'Horizon (days) from re-weighting', ylim = c(-0.2, 0.6))
lines(1:40, 1 - mod_10$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss]/
          Mod_9_60_NWP1$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss] ,
      col = 'grey')
lines(1:40, 1 - mod_15$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss]/
          Mod_9_60_NWP1$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss] ,
      col = 'grey')
lines(1:40, 1 - mod_18$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss]/
          Mod_9_60_NWP1$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss] ,
      col = 'red')
lines(1:40, 1 - mod_20$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss]/
          Mod_9_60_NWP1$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss] ,
      col = 'grey')
lines(1:40, 1 - mod_25$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss]/
          Mod_9_60_NWP1$Results[lead_time>60& lead_time<=100, .(mean_test_loss = mean(test_loss)), keyby = lead_time][,mean_test_loss] ,
      col = 'grey')
ax = pretty(1:40, 5)
axis(1, at = ax, labels = ax/4)
abline(h = 0)
legend("topright", lty = 1, cex=1.1, col = c('grey', 'grey', 'grey', 'red', 'grey', 'grey'),
       legend = c('K = 5', 'K = 10','K = 15', 'K = 18', 'K = 20', 'K = 25'), title = 'Model:')


# Pinball loss comparing K for different loss periods
l_5 = rep(NA,25)
l_7 = rep(NA,25)
l_10 = rep(NA,25)
l_15 = rep(NA,25)
for (k in 1:25){
    mod_name = paste0('mod_',k)
    mod = get(mod_name)
    l_5[k] = mod$Results[lead_time>60&lead_time<=80, mean(test_loss)]
    l_7[k] = mod$Results[lead_time>60&lead_time<=88, mean(test_loss)]
    l_10[k] = mod$Results[lead_time>60&lead_time<=100, mean(test_loss)]
    l_15[k] = mod$Results[lead_time>60&lead_time<=120, mean(test_loss)]
}


plot(l_5, ylim = c(6.4, 9.2), xlim = c(0, 27), main = 'Pinball loss comparing K for different loss periods',
     xlab = 'K - tuning parameter', cex.main = 1.6,  cex.lab = 1.3, ylab = 'Pinball Loss', xaxt = 'n')
points(l_7, col = 'purple')
points(l_10, col = 'red')
points(l_15, col = 'blue')
points(which.min(l_5), min(l_5), col = 'black', pch = 16)
points(which.min(l_7), min(l_7), col = 'purple', pch = 16)
points(which.min(l_10), min(l_10), col = 'red', pch = 16)
points(which.min(l_15), min(l_15), col = 'blue', pch = 16)

legend("bottomleft", lty = 1, cex=1.1, col = c('black', 'purple', 'red', 'blue'),
       legend = c('5 days', '7 days', '10 days', '15 days'), title = 'Loss period:')
ax = pretty(1:25, 5)
sq = exp(seq(log(0.000001), log(0.01), length.out = 25))
axis(1, at = ax, labels =     format(c(0, sq)[ax+1], scientific = TRUE, digits = 3))
