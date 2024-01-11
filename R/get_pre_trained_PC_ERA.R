#' Gets PCA matrix of temperature data for a specified time period.
#'
#' @param file String. File to store and look for PCA matrix.
#' @param X_mat ERA temperature matrix.
#' @param date_demand Date and demand matrix. We only care about the dates in this application.
#' @param start_train String. Start date training.
#' @param stop_train String. End date training.
#' @param start_test String. Start date test.
#' @param stop_test String. End date test.
#' @param NWP NWP matrix.
#' @param run_again Boolean. If FALSE and file exist function loads pre-made pca matrix.
#'
#' @return A list of data.tables
#' @export
#'
#' @examples get_pre_trained_PC_ERA(file='~/Desktop/Master2023/Data/PC_ERA', X_mat,
#' start_train = "1979-01-01", stop_train="1992-12-31", start_test = "1993-01-01", stop_test="2023-01-31",
#' NWP=NA,run_again =TRUE)
#'
#'
get_pre_trained_PC_ERA = function(file='~/Desktop/Master2023/Data/PC_ERA', X_mat,date_demand = date_demand,
                                  start_train = '2013-01-01', stop_train='2013-12-31',
                                  start_test = '2014-01-01', stop_test='2014-12-31',
                                  NWP=NA,run_again =FALSE){

    if (start_train <date_demand[1,date]){
        print(paste0('Error: Start date provided is earlier than earliest temp observations: ',date_demand[1,date]))
        return('Either change start date or load appropriate X_mat temperature observations.')
        }

    train_start = date_demand[(date ==start_train & hour == 1),I]
    train_stop = date_demand[(date ==stop_train & hour == 24),I]
    I_train = seq(from = train_start, to = train_stop)

    test_start = date_demand[(date ==start_test & hour == 1),I]
    test_stop = date_demand[(date ==stop_test & hour == 24),I]
    I_test = seq(from = (test_start), to =(test_stop))

    name = paste0('PC_ERA_',substring(as.character(year(start_train)),3),'_',
                  substring(as.character(year(stop_train)),3))
    file_name = paste0(file,'/',name,'.Rda')
    print(file_name)
    if(file.exists(file_name) & run_again ==FALSE){
        print('Loading PC files')
        assign(name, get(load(file_name)))
        return(get(name))
    }else{
        print('Creating PC file')
        out = get_pca(X_mat, date_demand = date_demand, I_train, I_test, p_comps = 2, NWP = NWP)
        assign(name, out)
        save(list = name, file = file_name)
        return (get(name))
    }
}
