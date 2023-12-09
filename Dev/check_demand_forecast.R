#X_mat
#date_demand

forc_start ='2014-01-01'
forc_end = '2023-01-01'
pred_win = 30
pred_lag= 15
train_y=5
reg_form = "volume ~  s(PC1)"
p_comps = 1
other_mods= NULL
comb = TRUE
custom = FALSE
incl_climatology = TRUE
no_pc = TRUE
gam_lasso = FALSE
cores = 29
Setup = "Rolling_test"


start = as.Date(forc_start)
end = as.Date(forc_end)
all_days = seq(start, end,  by = '1 days')
init_days_all = all_days[mday(all_days)==1]

# **** Run parallel cores ****
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1) #Set number of threads

init_days= init_days_all

## **** Step 1: Form the training and test datasets ****
init_day = init_days[i]
target_days = seq(init_day+ pred_lag, length.out = pred_win,  by = '1 days')
print(paste('Forecast made on:', init_day))

train_cutoff = seq(init_day,  length.out = 2, by = paste0('-',train_y, ' year'))[2]

I_train = date_demand[(date > train_cutoff) & (date < init_day), I]
I_test = date_demand[date %in% target_days, I]

if (p_comps > 0){
    pc_data  = get_pca(X_mat, I_train, I_test, p_comps) # **** SVD ****
    dt_train = pc_data$dt_train
    dt_test  = pc_data$dt_test
} else{
    dt_train = date_demand[I_train, ]
    dt_test  = date_demand[I_test, ]
}
dt_train= dt_train[!is.na(volume)]
print(paste('We are training on:', dim(dt_train)[1], 'observations'))
max_year_train = dt_train[, max(year)]
dt_test = dt_test[year > max_year_train, year := max_year_train] #pretend that new year is last year to avoid factor error
n = nrow(dt_test)


climatology = dt_train[, .(ave_vol = mean(volume)), by = .(paste0(yr_by_,i), 'hour')]



climatology[,.(paste0('yr_by_',i) =int, hour, ave_vol)]



climatology[int %in% dt_test[,get(paste0('yr_by_', i))],ave_vol]



climatology[,.(paste0('yr_by_',i)=int)]



for (j in 1:14){
    print(j)
    climatology = date_demand[, .(ave_vol = mean(volume, na.rm = TRUE)), by = c(paste0('yr_by_',j), "hour")]
    date_demand = merge(climatology, date_demand, by= c(paste0('yr_by_',j), 'hour') )
    setnames(date_demand, "ave_vol", paste0('clim_pred_',j))

    print(paste0("Climatology: ", date_demand[,sqrt(sum((volume- get(paste0('clim_pred_',j))) )^2)/n]))
}

a = Sys.time()
date_demand[,season := factor(season)]
print(Sys.time()-a)
