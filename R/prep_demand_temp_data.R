#' Finds, checks and loads temperature and demand data.
#'
#' @param include_na_volume TRUE includes dates with na volume data. Use TRUE for creating historical PC data.
#' "Some" will not remove sequences of missing within the scope of the demand volume. Use some so that creation of PC will not be affected.
#' @param path String File path.
#'
#' @return Returns a matrix with grid temperature data and a data.table with demand data.
#' @import ncdf4
#' @import data.table

#' @examples prep = prep_demand_temp_data(include_na_volume = TRUE,path = "~/NR/ClimateFutures/RenewableEnergy/Data/DemandForecasting/")
#' @examples prep = prep_demand_temp_data(include_na_volume = TRUE,path = "/Users/Eirik/Desktop/Master2023/Data/")
#' @examples prep = prep_demand_temp_data(include_na_volume = TRUE,path = "/nr/samba/user/esjavik/EnergyDemandData/")
#'
#' @export
prep_demand_temp_data = function(include_na_volume = TRUE, path = "/mn/kadingir/datascience_000000/eirikhsj/"){

    ##----- 1) Load and extract temperature data -----
    #f_nc = paste0(path, "/temperature_data.nc4") #Old dataset 2013-2021 only
    #f_nc = paste0(path, "era_historical_data.nc4") #New dataset 1978-2023
    f_nc = paste0(path, "era_historical_data_1979_2024.nc") #New dataset 1978-2023
    nc = ncdf4::nc_open(f_nc)
    X = ncdf4::ncvar_get(nc, "2m_temperature")
    hours = 1:24
    days = as.Date(ncvar_get(nc, "day"), origin = "1970-01-01")
    nc_close(nc)

    X_temp = matrix(as.vector(X),
                    nrow = dim(X)[3] * dim(X)[4], #hours times days
                    ncol = dim(X)[1] * dim(X)[2], #lon time lat
                    byrow = TRUE)

    dt_date = data.table(expand.grid(hour = hours, date = days)) #track of dates for temp data

    ##----- 2) Load and extract demand volume data -----
    f_demand = paste0(path, "/nordpool_volume_2013_2023.csv") # 2013 - 2023
    dt_demand = fread(f_demand)
    dt_demand[, hour := hour + 1]


    ##----- 3) Debug Demand volume (Summer/Winter time problem)-----
    nas = which(is.na(dt_demand$volume))           #Find missing bc we skip 1 hour ahead -> summertime
    dt_demand[, rowCount:= .N, by = .(date, hour)] #Find doubles bc we jump back 1 h -> wintertime
    doubles= which(dt_demand$rowCount >1)[seq(from = 1, to = length(which(dt_demand$rowCount >1)), by = 2)]


    if(any(dt_demand[nas,hour]!=3)){
        print("We haves na's in the demand data not related to summer time/winter time.")

        if(any(diff(nas)==1)){
            print("The na's form a sequence.")
            sq_nas = which(diff(nas)==1)
            if(nas[length(nas)-1]+1==nas[length(nas)]){
                sq_nas = c(sq_nas, length(nas))
            }
            nas_sq = nas[sq_nas]
            na_dates = dt_demand[nas_sq,unique(date)]
            nas = nas[-sq_nas]
            missing_seq = TRUE
        } else{
            missing_seq = FALSE
        }
    }

    if ( length(nas) == length(doubles) ){ #if data-cutoff in winter
        volume = dt_demand$volume[-nas]
        date = dt_demand$date[-doubles]
        hour = dt_demand$hour[-doubles]
    } else if ( length(nas) == length(doubles) + 1 ){ #if data-cutoff in summer
        volume = dt_demand$volume[-nas]
        date = dt_demand$date[-c(doubles, length(dt_demand$date))]
        hour = dt_demand$hour[-c(doubles, length(dt_demand$date))]
        print("I guess we cut off demand observations sometime during summer time, and didn't go back to winter time.")
    } else{
        print("We have other sources of na's in the demand data than summer time/winter time.")
    }

    dt_demand = data.table(volume, date, hour)
    date_demand = dt_demand[dt_date, on=c("hour", "date")]  #Combine information

    ##----- 4) Debug Temperature Data (X) -----
    n_nas = rowSums(is.nan(X_temp))
    if(any(n_nas > 0)){
        w_missing = which(n_nas > 0)
        X_mat = X_temp[-w_missing, ]
        date_demand = date_demand[-w_missing, ]
    }else{
        X_mat = X_temp
    }

    if (nrow(date_demand) == nrow(X_mat)){
        print('We are good to go!')
    } else if(nrow(X_mat) > nrow(date_demand)){
        diff = nrow(X_mat)  - nrow(date_demand)
        diffseq = seq(nrow(X_mat), nrow(X_mat)-diff+1, -1)
        X_mat = X_mat[-(diffseq),]
    } else{
        print('We have insufficient temperature observations.')
    }

    date_demand[,I:= 1:.N] #set index

    ##----- 5) Specify time covariates-----

    date_demand[,week:= data.table::isoweek(date)] # jan 1st might get week nr 52 or 53
    date_demand[,week1 := cos(2*pi * week/52) ]
    date_demand[,week2 := sin(2*pi * week/52) ]

    date_demand[,month:= month(date)]
    date_demand[,month1 := cos(2*pi * month/12) ]
    date_demand[,month2 := sin(2*pi * month/12) ]

    date_demand[,hour1 := cos(2*pi * hour/24) ]
    date_demand[,hour2 := sin(2*pi * hour/24) ]

    date_demand[,year:= year(date)]

    date_demand[,w_day:= as.numeric(format(date, "%u"))] # ISO 8601
    date_demand[,w_day1 := cos(2*pi * w_day/7) ]
    date_demand[,w_day2 := sin(2*pi * w_day/7) ]

    date_demand[, season:= ifelse(month(date) %in% c(12,1,2), 1,
                                  ifelse(month(date) %in% c(3,4,5), 2,
                                         ifelse(month(date) %in% c(6,7,8), 3, 4)))]

    date_demand[,season1 := round(cos(2*pi * season/4), digits = 7)]
    date_demand[,season2 := round(sin(2*pi * season/4), digits = 7)]


    date_demand[, year_day := as.integer(format(date, format ="%j"))]
    date_demand[, paste0('yr_by_', seq(from = 1, to =14)) := lapply(1:14, function(c) floor(year_day/c))]


    date_demand[,hourly_mean_grid := rowMeans(X_mat)]

    ##----- 6) Do we want to return NA volume observations? -----
    out = list()
    if (include_na_volume== FALSE){
        d_d = date_demand[!is.na(volume)] #Select only rows with volume obs
        out$X_mat = X_mat[d_d$I,]
        out$date_demand = d_d[,I:= 1:.N] #reset index
    } else if (include_na_volume== 'Some'){
        d_d = date_demand[first(dt_demand[,date])<= date & last(dt_demand[,date])>= date] #Keep NAs if they are in time period with demand data.
        out$X_mat = X_mat[d_d$I,]
        out$date_demand = d_d[,I:= 1:.N]
    } else{
        out$X_mat = X_mat
        out$date_demand = date_demand
    }

    ##----- 7) Final Check and Print-----
    if (nrow(out$date_demand) == nrow(out$X_mat)){
        print(paste0('Demand volume data from ', dt_demand[1,date], ' to ', dt_demand[dim(dt_demand)[1],date],'.' ))
        print(paste0('Temperature data from ', out$date_demand[1,date], ' to ', out$date_demand[dim(out$date_demand)[1],date],'.'))
        print(paste0('We have ', dim(out$X_mat)[1], ' hourly observations at ', dim(out$X_mat)[2], ' grid points.'))
        if(missing_seq == TRUE){
            print(paste0('We have a sequence of ', length(nas_sq), ' missing energy volume values.'))
            if (include_na_volume ==FALSE){

                print(paste0('Since include_na_volume = FALSE, we have not included these or any other NAs in the volume data'))
            } else if( include_na_volume == 'Some'){
                print(paste0('Since include_na_volume = `Some`, we have included these, which are found on ',na_dates,' but not other NAs'))
            } else{
                print(paste0('Since include_na_volume = TRUE, we have included these missing data, which are found on ',na_dates))
            }
        }


    } else{
        print("Something went wrong")
    }
    return(out)
}
