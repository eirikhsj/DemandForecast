

## COPULAS
#
# x1 = rnorm(100, 10, 1)
# x2 = rnorm(100, 12, 2)
# y1 = rnorm(100, x1 + 0.1*x2, 0.2)
# y2 = rnorm(100, x2, 2.5)
# new_data = data.table(x1 = c(10,20), x2 = c(13,10))
#
#
# tau_vals = 0.5
# qs = rq(y1 ~ x1 + x2, tau = tau_vals)
# predict(qs, new_data)
#
#
# sigma = cor(matrix(c(qnorm(c(0.1, 0.2, 0.6, 0.23)),qnorm(c(0.11, 0.12, 0.16, 0.93)), qnorm(c(0.19, 0.92, 0.96, 0.23))),ncol = 3, byrow = FALSE))
# Sim = data.table(mvrnorm(n = 1000, mu = c(0,0,0), Sigma = sigma))
# Sim[, lapply(.SD, pnorm)]
#
# predict(qs, new_data, interval = 'none')
#
# #Set quantile range
# m = 1000
# N = 1000
# tau_vals = seq(0, 1, 1/m)
#
# #Create predictions based on training interval
# mod_copula = rq(PC1 ~ NWP1_90, data = ERA_NWP[init_date> as.Date("2007-01-01)") & lead_time <=1*4,], tau = tau_vals)
# pred_mat_train = predict(mod_copula, interval = 'none')
# pred_mat_test = predict(mod_copula, ERA_NWP[init_date== as.Date("2007-01-01)") & lead_time <=1*4,], interval = "none")
# pred_mat_test = unname(pred_mat_test)
# matplot(t(pred_mat_test), type = 'l')
#
# #Check where the observed PC1 values land in the predictive cdf (i.e. find the quantile)
# check_test =  data.table(t(pred_mat_test > ERA_NWP[init_date== as.Date("2007-01-01)") & lead_time <=1*4,PC1]))
# check_train =  data.table(t(pred_mat_train > ERA_NWP[init_date< as.Date("2007-01-01)") & lead_time <=1*4,PC1]))
#
# q_train = t(check_train[, lapply(.SD, function(x) ifelse(is.na(match(TRUE, x)),(m-1),
#                                              ifelse(match(TRUE, x)>2,match(TRUE, x)-2, 1)) )]/m)
# q_test = t(check_test[, lapply(.SD, function(x) ifelse(is.na(match(TRUE, x)),(m-1),
#                                                   ifelse(match(TRUE, x)>2,match(TRUE, x)-2, 1)) )]/m)
# q_mat = data.table(matrix(q_train, ncol = 4, byrow = TRUE))
#
# #Normalize
# z_train = q_mat[, lapply(.SD, function(x)  qnorm(x))]
#
# #Find correlation between the time points e.g. 1-4
# sigma = cor(z_train)
#
# #Simulate from multivariate normal with mu = 0 and sigma = sigma
# z_sim = data.table(mvrnorm(n = N, mu = rep(0,dim(sigma)[1]), Sigma = sigma))
#
# #Find the quantile index
# q2 = z_sim[, lapply(.SD, function(x)  round(pnorm(x),log(m, 10))*m)]
#
# #q2 = round(pnorm(as.matrix((z_sim))),log(m, 10))*m
#
# q2 = as.matrix(q2)
#
# # Get back U values from q-index
# U = data.table()
# for (b in 1:(1*4)){
#     ind = q2[,b]
#     print(max(ind))
#     vec = pred_mat_test[b,ind+1]
#     U[, paste0("col", b)  := vec ]
# }
#
# #Find the simulated means
# W = rowMeans(U)
#
# #Utilize the actual predictions??
# pred_mat_09 = predict(mod_copula, ERA_NWP[init_date== as.Date("2007-01-01)") & lead_time <=1*4,], interval = "prediction", quantiles = 0.90)[,901]
# copula_quant =  (match(TRUE, quantile(W, probs = tau_vals) >mean(pred_mat_09))-2 )/m
#
# #Utilize the copula quantile
# copula_quant = quantile(W, probs = 0.9)
#
# test_loss = pinball_loss(0.9, copula_quant, ERA_NWP[init_date== as.Date("2007-01-01)") & lead_time <=1*4,mean(PC1)])
#
#


