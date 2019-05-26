

# penalty functions
intexpit <- function(x){ log(1+exp(x)) }
expit <- function(x){ exp(x)/(1+exp(x)) }

pen_hi <- function(x,tt,aa){ -tt*intexpit(x-aa) }
dpen_hi <- function(x,tt,aa){ -tt*expit(x-aa) }
ddpen_hi <- function(x,tt,aa){ -tt*expit(x-aa)/(1+exp(x-aa)) }
    
pen_lo <- function(x,tt,aa){ -tt*intexpit(-x+aa) }
dpen_lo <- function(x,tt,aa){ +tt*expit(-x+aa) }
ddpen_lo <- function(x,tt,aa){ -tt*expit(-x+aa)/(1+exp(-x+aa)) }

pen_loglo <- function(x,tt,aa){ pen_lo(log(x),tt,aa) }
dpen_loglo <- function(x,tt,aa){ dpen_lo(log(x),tt,aa)/x }
ddpen_loglo <- function(x,tt,aa){ 
    ddpen_lo(log(x),tt,aa)/x^2 - dpen_lo(log(x),tt,aa)/x^2 
}
