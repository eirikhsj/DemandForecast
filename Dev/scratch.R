arg1 = X_mat
arg2 = date_demand
arg3 = 30
arg4 = 15
arg5 = 5
arg6 = "volume ~ as.factor(hour)"
arg7 = 0
arg8 = FALSE
arg9 = FALSE
arg10 = FALSE
arg11 = FALSE
arg12 = 48
forc_start = '2016-01-01'
forc_end = '2023-01-01'

last_poss_pred = range(date_demand$date)[2] - 30 - 15

#1) Fix dates
start = as.Date(forc_start)
end = as.Date(forc_end)
all_days = seq(start, end,  by = '1 days')
init_days_all = all_days[mday(all_days)==1]

detailed_results = list()
for (i in seq_along(init_days_all)){
    print(i)
    detailed_results[[i]] = Rolling(i,X_mat = arg1, date_demand = arg2,init_days= init_days_all, pred_win = arg3, pred_lag= arg4, train_y=arg5,
            reg_form= arg6, p_comps= arg7, other_mods= arg8, comb = arg9, custom = arg10,
            incl_climatology = arg11)
    }
