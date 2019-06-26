

# penalty functions


#' expit function and integral of expit function
#' 
#' @param x argument to expit or intexpit function
expit <- function(x){ exp(x)/(1+exp(x)) }

#' @rdname expit
intexpit <- function(x){ log(1+exp(x)) }

#' penalize large values of parameter: penalty, 1st deriative, 2nd derivative
#'
#' @param x argument to penalty
#' @param tt scale parameter of penalty
#' @param aa location parameter of penalty
pen_hi <- function(x,tt,aa){ -tt*intexpit(x-aa) }

#' @rdname pen_hi
dpen_hi <- function(x,tt,aa){ -tt*expit(x-aa) }

#' @rdname pen_hi
ddpen_hi <- function(x,tt,aa){ -tt*expit(x-aa)/(1+exp(x-aa)) }
    
#' penalize small values of parameter: penalty, 1st deriative, 2nd derivative
#'
#' @param x argument to penalty
#' @param tt scale parameter of penalty
#' @param aa location parameter of penalty
pen_lo <- function(x,tt,aa){ -tt*intexpit(-x+aa) }

#' @rdname pen_lo
dpen_lo <- function(x,tt,aa){ +tt*expit(-x+aa) }

#' @rdname pen_lo
ddpen_lo <- function(x,tt,aa){ -tt*expit(-x+aa)/(1+exp(-x+aa)) }


#' penalize small values of log parameter: penalty, 1st deriative, 2nd derivative
#'
#' @param x argument to penalty
#' @param tt scale parameter of penalty
#' @param aa location parameter of penalty
pen_loglo <- function(x,tt,aa){ pen_lo(log(x),tt,aa) }

#' @rdname pen_loglo
dpen_loglo <- function(x,tt,aa){ dpen_lo(log(x),tt,aa)/x }

#' @rdname pen_loglo
ddpen_loglo <- function(x,tt,aa){ 
    ddpen_lo(log(x),tt,aa)/x^2 - dpen_lo(log(x),tt,aa)/x^2 
}
