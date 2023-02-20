#' Multi-use function for getting pca-components.
#'
#' Fetches the desired number of pca-components and returns these within a data.table together with time covariates.
#' If given an NWP-matrix, the U and mu from the decomposition of ERA-temperature (X_mat) is used to compute and select factor loading of NWP forecasts.
#' If U and mu are pre-supplied the computation of NWP-forecasts is greatly sped up.
#' The function is used iteratively in demand_forecast.R
#' X_mat is supplied by prep_demand_temp_data.R, while I_train and I_test are normally chosen by the rolling cv function (e.g. demand_forecast.R).
#'
#' @param X_mat Temperature grid matrix.
#' @param I_train Indexes of training set.
#' @param I_test Indexes of test set.
#' @param p_comps Number of principle components to return.
#' @param NWP NWP matrix (dim  = 22 21 500 202)
#' @param U Use pre-trained U matrix to speed up time. If supplied we do not get dt_train or dt_test only the NWP PC matrix.
#' @param mu Pre-trained mean vector, corresponding to U.
#'
#' @return Training (dt_train) and test (dt_test) sets with pca components of ERA temperature field and matching time covariates.
#' U-matrix (462x462), and mu-vector (1xnrow(X_mat)). NWP-matrix with pca of NWP temperature forecasts.
#' @export
#'
#' @examples get_pca(X_mat, I_train, I_test, p_comps)
#'
#'

get_pca = function(X_mat, I_train, I_test, p_comps, NWP=NA,  U = NA, mu = NA){
    print('------- Running get_pca -------')
    if (length(U) == 1){

        X_train = X_mat[I_train,]
        ##--- Center the data ---
        mu = colMeans(X_train)
        X_center = matrix(NA, dim(X_train)[1], dim(X_train)[2])
        for(j in 1:dim(X_train)[2]){
            X_center[, j] = X_train[, j] - mu[j]
        }

        ##---- Perform the SVD ---
        Sigma = t(X_center) %*% X_center
        l_svd = svd(Sigma)
        U = l_svd$u[, 1:p_comps]

        ##--- Get the factor loadings for training dataset---
        Fact_train = X_center %*% U
        colnames(Fact_train) = paste0("PC", 1:dim(data.frame(U))[2])
        dt_train = data.table(date_demand[I_train, ], Fact_train)

        ##----- Form the test dataset ---
        X_test = X_mat[I_test,]
        X_test_center = X_test - mu
        Fact_test = X_test_center %*% U
        colnames(Fact_test) = paste0("PC", 1:dim(data.frame(U))[2])
        dt_test = data.table(date_demand[I_test, ], Fact_test)

    } else{  #We do not return dt_train and dt_test
        dt_train = NA
        dt_test = NA
    }

    ## ---- Form factor loading of NWP data ---
    if (is.null(dim(NWP)) == FALSE){
        members = dim(NWP)[4]
        NWP_PC_mat = array(NA, dim = c(500, p_comps, members))

        #Iterate through ensemble members, creating factor loadings
        for (ens_memb in 1:members){
            X_ensamble = NWP[,,,ens_memb]

            #Create 500x462 matrix:
            X_mat_500 = matrix(as.vector(X_ensamble),
                           nrow = dim(X_ensamble)[3], # nr of 6-hour spaced obs
                           ncol = dim(X_ensamble)[1] * dim(X_ensamble)[2], #lon time lat
                           byrow = TRUE)
            X_NWP_center = X_mat_500 - mu
            Fact_NWP = X_NWP_center %*% U        #Here we use U the information from the historical data
            NWP_PC_mat[,,ens_memb] = Fact_NWP
        }
    } else{
        NWP_PC_mat = NA
    }

    ##----- Return ---
    out = list()
    out$U = U
    out$mu = mu
    out$dt_train = dt_train
    out$dt_test = dt_test
    out$NWP_PC_mat = NWP_PC_mat
    return(out)
}
