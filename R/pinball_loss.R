#' @title Pinball loss
#' @name pinball_loss
#' @param quantile Quantile of interest
#' @param pred Predicted values
#' @param obs Observed values
#'
#' @return Pinball Loss
#' @export
#'
#' @examples
#' pinball_loss(quantile = 0.9, pred = 2, obs = 2)
#'
pinball_loss = function(quantile, pred, obs){
    #Input is quantile q and predicted values pred, and observed values obs
    l = rep(0, length(pred))

    for (i in 1:length(pred)){
        if (is.na(pred[i])| is.na(obs[i])){
            l[i]= NA
        } else if(pred[i]< obs[i] ){
            l[i] = (obs[i]-pred[i])*quantile
        } else{
            l[i] =  (pred[i]-obs[i])*(1-quantile)
        }
    }

    return (l)
}
