bayessurvreg1.checkStore <- function(store)
{
  if (is.null(store$y)) store$y <- FALSE
  if (is.null(store$r)) store$r <- FALSE
  if (is.null(store$b)) store$b <- FALSE
  if (is.null(store$u)) store$u <- FALSE
  if (is.null(store$MHb)) store$MHb <- FALSE
  if (is.null(store$regresres)) store$regresres <- FALSE

  return(store)  
}  
