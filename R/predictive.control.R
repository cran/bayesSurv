predictive.control <- function(predict, store, quantile)
{
  if (is.null(predict$Et)) predict$Et <- FALSE 
  if (is.null(predict$t)) predict$t <- FALSE
  if (is.null(predict$Surv)) predict$Surv <- FALSE
  if (is.null(predict$hazard)) predict$hazard <- FALSE
  if (is.null(predict$cum.hazard)) predict$cum.hazard <- FALSE
  if (!(predict$Et || predict$t || predict$Surv || predict$hazard || predict$cum.hazard))
    stop("Nothing to be predicted.")


  if (is.null(store$Et)) store$Et <- FALSE
  if (is.null(store$t)) store$t <- FALSE
  if (is.null(store$Surv)) store$Surv <- FALSE
  if (is.null(store$hazard)) store$hazard <- FALSE
  if (is.null(store$cum.hazard)) store$cum.hazard <- FALSE

  if (!predict$Et) store$Et <- FALSE
  if (!predict$t) store$t <- FALSE
  if (!predict$Surv) store$Surv <- FALSE
  if (!predict$hazard) store$hazard <- FALSE
  if (!predict$cum.hazard) store$cum.hazard <- FALSE

  if (sum(quantile < 0 | quantile > 1)) stop("Quantiles must lie between 0 and 1.")

  back <- list(predict = predict, store = store)
  return(back)
}

