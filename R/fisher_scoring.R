
#' test likelihood object for NA or Inf values
#' 
#' @param likobj likelihood object
test_likelihood_object <- function(likobj){
    
    pass <- TRUE
    allvals <- c( likobj$loglik, likobj$grad, c(likobj$info) )
    if( sum(is.na(allvals)) > 0  ||  sum( abs(allvals) == Inf ) > 0 ){
        pass <- FALSE
    }
    return(pass)
}


#' compute condition number of matrix
#' 
#' @param info matrix
condition_number <- function(info){
    # assumes that information matrix has finite numbers in it
    if(max(diag(info))/min(diag(info)) > 1e6){
        return( max(diag(info))/min(diag(info)) )
    } else {
        ee <- eigen(info)
        return( max(ee$values)/min(ee$values) )
    }
}    


#' Fisher scoring algorithm
#' 
#' @param likfun likelihood function, returns likelihood, gradient, and hessian
#' @param start_parms starting values of parameters
#' @param link link function for parameters (used for printing)
#' @param silent TRUE/FALSE for suppressing output
#' @param convtol convergence tolerance on step dot grad
fisher_scoring <- function( likfun, start_parms, link, silent = FALSE, convtol = 1e-4 ){
    
    # link functions passed for printing purposes
    maxit <- 40
    
    # function for checking wolfe conditions
    wolfe_check <- function(likobj0,likobj1,logparms,step,both){
        c1 <- 1e-4
        c2 <- 0.9
        tol <- 1e-6
        ll0 <- likobj0$loglik
        gr0 <- likobj0$grad
        ll1 <- likobj1$loglik
        gr1 <- likobj1$grad
        if(!both){
            satfd <- ll1 <= ll0 + c1*crossprod(step,gr0) + tol
        } else {
            satfd <- ll1 <= ll0 + c1*crossprod(step,gr0) + tol &&
                 -crossprod(step,gr1) <= -c2*crossprod(step,gr0) + tol
        }
        return(satfd)
    }
    
    # evaluate function at initial values
    logparms <- start_parms
    likobj <- likfun(logparms)
    
    # test likelihood object    
    if( !test_likelihood_object(likobj) ){
        logparms <- 0.1*logparms
        likobj <- likfun(logparms)
    }
    
    # assign loglik, grad, and info
    loglik <- likobj$loglik        
    grad <- likobj$grad
    info <- likobj$info
    
    # add a small amount of regularization
    diag(info) <- diag(info) + 0.1*min(diag(info))

    # print some stuff out
    if(!silent){
        cat(paste0("Iter ",0,": \n"))
        cat("pars = ",  paste0(round(link(logparms),4)), "  \n" )
        cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
        cat("grad = ")
        cat(as.character(round(-grad,3)))
        cat("\n\n")
    }
    
    for(j in 1:maxit){
        
        likobj0 <- likobj
        
        # if condition number of info matrix large, then regularize
        tol <- 1e-8
        if (condition_number(info) > 1 / tol) {
            if (!silent) cat("Cond # of info matrix > 1/tol \n")
            eiginfo <- eigen(info)
            whichsmall <- which( eiginfo$values / max(eiginfo$values) < tol )
            eiginfo$values[whichsmall] <- tol*max(eiginfo$values)
            info <-
                eiginfo$vectors %*% diag(eiginfo$values) %*% t(eiginfo$vectors)
        }

        # calculate fisher step 
        step <- - solve(info, grad)
        
        # if step size large, then make it smaller
        if (mean(step^2) > 1) {
            if(!silent) cat("@@\n")
            step <- step/sqrt(mean(step^2))
        }
        
        # take step and calculate loglik, grad, and info
        newlogparms <- logparms + step
        likobj <- likfun(newlogparms)
        
        # check for Inf, NA, or NaN
        cnt <- 1
        while (!test_likelihood_object(likobj)) {
            if (!silent) cat("inf or na or nan in likobj\n")
            step <- 0.5 * step
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
            if (cnt == 10) { stop("could not avoid inf, na or nan\n") }
        }
        
        # Check the wolfe conditions
        # if not satisfied, shrink fisher step
        # after some iterations, switch to gradient
        cnt <- 1
        no_decrease <- FALSE
        both <- FALSE
        while (!wolfe_check(likobj0,likobj,logparms,newlogparms-logparms,both) &&
                !no_decrease ){
            mult <- 0.5 
            step <- mult * step
            if (cnt == 6) { # switch to gradient
                if(!silent) cat("**\n") 
                step <- - 0.1*grad/sqrt(sum(grad^2)) 
                both <- FALSE
            }
            if(!silent) cat("**\n") 
            if ( sqrt(sum(step^2)) < 1e-4 ){ no_decrease <- TRUE }  # maybe we should throw error here?
            newlogparms <- logparms + step
            likobj <- likfun(newlogparms)
            cnt <- cnt + 1
        }
        stepgrad <- c(crossprod(step,grad))
        
        # redefine logparms, loglik, grad, info
        logparms <- logparms + step
        loglik <- likobj$loglik        
        grad <- likobj$grad
        info <- likobj$info

        # print some stuff out
        if(!silent){
            cat(paste0("Iter ",j,": \n"))
            cat("pars = ",  paste0(round(link(logparms),4)), "  \n" )
            cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
            cat("grad = ")
            cat(as.character(round(grad,4)),"\n")
            cat("step dot grad = ",stepgrad,"\n")
            cat("\n")
        }
        
        # if gradient is small, then break and return the results        
        if( abs(stepgrad) < convtol || no_decrease ){
            break
        }
    }

    # collect information to return
    betahatinfo <- likobj        
    betahat <- as.vector(betahatinfo$betahat)
    betacov <- solve(betahatinfo$betainfo)
    sebeta <- sqrt(diag(betacov))
    tbeta <- betahat/sebeta

    ret <- list(
        covparms = link(logparms), 
        betahat = betahat, 
        sebeta = sebeta,
        betacov = betacov,
        tbeta = tbeta,
        loglik = loglik,
        no_decrease = no_decrease,
        grad = likobj$grad
    )
    return(ret)
}






