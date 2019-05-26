



fisher_scoring_old <- function( likfun, start_parms, link, dlink, invlink, silent = FALSE ){
    
    maxit <- 50
    logparms <- invlink(start_parms)
    if(!silent) print(logparms)
    likobj <- likfun(link(logparms))        
    
    for(j in 1:maxit){
        
        if(!silent) print(round(link(logparms),4))
        loglik <- likobj$loglik        
        grad <- -c(likobj$grad)*dlink(logparms)
        info <- likobj$info*outer(dlink(logparms),dlink(logparms))
        # testing a solution to singular information matrix
        #print(min(eigen(info)$values))

        if(!silent){
            cat(paste0("Iter ",j,": "))
            cat(paste0( "loglik = ", round(loglik,6), "  "))
            cat("grad = ")
            cat(as.character(round(grad,3)))
            cat("\n")
        }
        if( max(abs(info)) == Inf || sum(is.nan(info)) > 0 ) {
            fisherstep <- 0.01 * grad / sum(grad ^ 2)
        } else {
            ee <- eigen(info)
            # do a regular fisher step if info matrix is well behaved
            if (min(ee$values) / max(ee$values) > 1e-12) {
                fisherstep <- solve(info, grad)
                # scale it down if the step is too big
                if (sum(fisherstep ^ 2) > 1) {
                    fisherstep <- fisherstep / sqrt(sum(fisherstep ^ 2))
                }
                # otherwise move in the direction of the gradient, scaled by marg info
            } else {
                # try generalized inverse
                whichbig <- which(ee$values / max(ee$values) > 1e-12)
                cat("@@")
                print(whichbig)
                geninv <-
                    ee$vectors[, whichbig,drop=FALSE] %*% 
                    as.matrix(diag(1 / ee$values[whichbig])) %*%
                    t(ee$vectors[, whichbig,drop=FALSE])
                fisherstep <- c(geninv %*% grad)
            }
        }
        
        # compute the new likelihood, gradient, and fisher info
        likobj <- likfun(link(logparms - fisherstep))   
        # check for problems
        # check to see if gradient or loglik has bad values in it
        # maybe we should also check the information matrix
        # currently that's checked above.
        cnt <- 0
        while( sum(is.nan(likobj$grad)) > 0 || sum(is.na(likobj$grad)) > 0 
            || is.na(likobj$loglik) || is.nan(likobj$loglik) ){
            cnt <- cnt + 1
            if(cnt == 8){ stop("gradient or loglik has NAs or NaNs") }
            fisherstep <- 0.5*fisherstep
            likobj <- likfun(link(logparms - fisherstep))
            
        }
            
        # if loglikelihood is increased, go with the fisher step
        # if not, do a line search along the gradient
        
        # change this back!
        # maybe we don't need this at all
        if( likobj$loglik > loglik - 200 && !is.nan(likobj$loglik) ){
            if(likobj$loglik < loglik){ if(!silent){ cat("## went down, but accepted ## \n") }}
            logparms <- logparms - fisherstep
        } else {
            stepsize <- 1/sqrt(sum(grad^2))
            while( is.nan(likobj$loglik) || likobj$loglik <= loglik ){
                gradstep <- stepsize*grad
                likobj <- likfun(link(logparms - gradstep))        
                
                if(likobj$loglik > loglik && !is.nan(likobj$loglik)){
                    logparms <- logparms - gradstep
                }
                stepsize <- stepsize/2
                if(stepsize < 1e-9){ # if stepsize gets too small
                    cat("could not decrease along gradient \n")
                    break
                }
            }
        }
        # if gradient is small, then break and return the results        
        if( max(abs(grad)) < 5e-3 ){
            break
        }
    }

    # collect information to return
    betahatinfo <- likfun(link(logparms))        
    betahat <- as.vector(betahatinfo$betahat)
    names(betahat) <- colnames(X)
    betacov <- solve(betahatinfo$betainfo)
    sebeta <- sqrt(diag(betacov))
    tbeta <- betahat/sebeta
    #df <- length(y) - ncol(X) # assumes X is full rank!
    #pbeta <- 2*( 1 - pt( abs(tbeta), df ) )
    
    ret <- list(
        covparms = link(logparms), 
        betahat = betahat, 
        sebeta = sebeta,
        betacov = betacov,
        tbeta = tbeta,
        #pbeta = pbeta,
        loglik = loglik
    )
    return(ret)
    
}



test_likelihood_object <- function(likobj, logparms, dlink){
    
    pass <- TRUE
    allvals <- c( 
        likobj$loglik, 
        likobj$grad*dlink(logparms), 
        c( likobj$info*outer(dlink(logparms),dlink(logparms)) )
    )
    if( sum(is.na(allvals)) > 0  ||  sum( abs(allvals) == Inf ) > 0 ){
        pass <- FALSE
    }
    return(pass)
    
}


test_likelihood_object5 <- function(likobj){
    
    pass <- TRUE
    allvals <- c( likobj$loglik, likobj$grad, c(likobj$info) )
    if( sum(is.na(allvals)) > 0  ||  sum( abs(allvals) == Inf ) > 0 ){
        pass <- FALSE
    }
    return(pass)
}



condition_number <- function(info){
    # assumes that information matrix has finite numbers in it
    if(max(diag(info))/min(diag(info)) > 1e6){
        return( max(diag(info))/min(diag(info)) )
    } else {
        ee <- eigen(info)
        return( max(ee$values)/min(ee$values) )
    }
}    


fisher_scoring2 <- function( likfun, start_parms, link, dlink, invlink, silent = FALSE ){
    
    
    maxit <- 50
    
    # evaluate function at initial values
    logparms <- invlink(start_parms)
    likobj <- likfun(link(logparms))
    
    if( !test_likelihood_object(likobj,logparms,dlink) ){
        logparms <- 0.1*logparms
        likobj <- likfun(link(logparms))
    }
    
    loglik <- likobj$loglik        
    grad <- -c(likobj$grad)*dlink(logparms)
    info <- likobj$info*outer(dlink(logparms),dlink(logparms))
    
    # print some stuff out
    if(!silent){
        cat(paste0("Iter ",0,": \n"))
        cat("pars = ",  paste0(round(link(logparms),8)), "  \n" )
        cat(paste0("loglik = ", round(loglik,6),         "  \n"))
        cat("grad = ")
        cat(as.character(round(grad,6)))
        cat("\n\n")
    }


    for(j in 1:maxit){
        
        # transform grad and info
        loglik <- likobj$loglik        
        grad <- -c(likobj$grad)*dlink(logparms)
        info <- likobj$info*outer(dlink(logparms),dlink(logparms))

        # if condition number large, then regularize
        tol <- 1e-8
        if( condition_number(info) > 1/tol ){
            cat("&&\n")
            eiginfo <- eigen(info)
            whichsmall <- which( eiginfo$values/max(eiginfo$values) < tol )
            eiginfo$values[whichsmall] <- tol*mean(eiginfo$values)
            info <- eiginfo$vectors %*% diag(eiginfo$values) %*% t(eiginfo$vectors)
        }

        # calculate fisher step and criterion at new parameter        
        fisherstep <- solve(info, grad)
        if( sum(fisherstep^2) > 4 ){ 
            cat("**\n")
            fisherstep <- 2*fisherstep/sqrt(sum(fisherstep^2)) 
        }
        #cat("$$\n")
        #print(round(fisherstep,6))
        likobj <- likfun(link(logparms - fisherstep))   
        
        # check to make sure all values are reasonable
        # if not, regularize information further and compute new step
        cnt <- 1
        while( !test_likelihood_object(likobj, logparms-fisherstep, dlink) ){
            cat("inf or na or nan in likobj")
            info <- info + 2^cnt*0.001*max(info)*diag(nrow(info))
            cnt <- cnt+1         
            fisherstep <- solve(info, grad)
            likobj <- likfun(link(logparms - fisherstep))   
            if( cnt == 8 ){ stop("regularized fisher scoring failed") }
        }            
        
        logparms <- logparms - fisherstep
        
        # print some stuff out
        if(!silent){
            cat(paste0("Iter ",j,": \n"))
            cat("pars = ",  paste0(round(link(logparms),8)), "  \n" )
            cat(paste0("loglik = ", round(loglik,6),         "  \n"))
            cat("grad = ")
            cat(as.character(round(grad,6)))
            cat("\n\n")
        }

        
        # if gradient is small, then break and return the results        
        if( max(abs(grad)) < 5e-3 ){
            break
        }
    }

    # collect information to return
    betahatinfo <- likfun(link(logparms))        
    betahat <- as.vector(betahatinfo$betahat)
    betacov <- solve(betahatinfo$betainfo)
    sebeta <- sqrt(diag(betacov))
    tbeta <- betahat/sebeta
    #df <- length(y) - ncol(X) # assumes X is full rank!
    #pbeta <- 2*( 1 - pt( abs(tbeta), df ) )
    
    ret <- list(
        covparms = link(logparms), 
        betahat = betahat, 
        sebeta = sebeta,
        betacov = betacov,
        tbeta = tbeta,
        #pbeta = pbeta,
        loglik = loglik
    )
    return(ret)
}


fisher_scoring3 <- function( likfun, start_parms, link, dlink, invlink, silent = FALSE ){
    
    
    maxit <- 50
    
    wolfe_check <- function(likobj0,likobj1,logparms,step,both){
        c1 <- 1e-4
        c2 <- 0.5
        ll0 <- -likobj0$loglik
        gr0 <- -c(likobj0$grad)*dlink(logparms)
        ll1 <- -likobj1$loglik
        gr1 <- -c(likobj1$grad)*dlink(logparms + step)
        #print(ll1)
        #print(c(ll0+c1*crossprod(step,gr0)))
        #print(c(-crossprod(step,gr1)))
        #print(c(-c2*crossprod(step,gr0)))
        if(!both){
            satfd <- ll1 <= ll0 + c1*crossprod(step,gr0)
        } else {
            satfd <- ll1 <= ll0 + c1*crossprod(step,gr0) &&
                 -crossprod(step,gr1) <= -c2*crossprod(step,gr0)
        }
        return(satfd)
    }
    
    # evaluate function at initial values
    logparms <- invlink(start_parms)
    likobj <- likfun(link(logparms))
    
    if( !test_likelihood_object(likobj,logparms,dlink) ){
        logparms <- 0.1*logparms
        likobj <- likfun(link(logparms))
    }
    
    # transform grad and info
    loglik <- likobj$loglik        
    grad <- -c(likobj$grad)*dlink(logparms)
    info <- likobj$info*outer(dlink(logparms),dlink(logparms))

    # print some stuff out
    if(!silent){
        cat(paste0("Iter ",0,": \n"))
        cat("pars = ",  paste0(round(link(logparms),8)), "  \n" )
        cat(paste0("loglik = ", round(loglik,6),         "  \n"))
        cat("grad = ")
        cat(as.character(round(grad,3)))
        cat("\n\n")
    }
    
    for(j in 1:maxit){
        
        likobj0 <- likobj
        
        # if condition number large, then regularize
        tol <- 1e-8
        if( condition_number(info) > 1/tol ){
            if(!silent) cat("&&\n")
            eiginfo <- eigen(info)
            #whichsmall <- which( eiginfo$values/max(eiginfo$values) < tol )
            eiginfo$values <- eiginfo$values + tol*max(eiginfo$values)
            info <- eiginfo$vectors %*% diag(eiginfo$values) %*% t(eiginfo$vectors)
        }
        
        
            # calculate fisher step and criterion at new parameter        
            step <- solve(info, grad)
            if( sum(step^2) > 1 ){ 
                step <- step/sqrt(sum(step^2)) 
            }
            likobj <- likfun(link(logparms - step))   
        
            # check for Inf, NA, or NaN
            cnt <- 1
            while( !test_likelihood_object(likobj, logparms-step, dlink) ){
                if(!silent) cat("inf or na or nan in likobj")
                step <- 0.5*step;
                likobj <- likfun(link(logparms - step))   
                if( cnt == 8 ){ stop("could not avoid info, na or nan") }
            }
        
            # Check the wolfe conditions
            # if satisfied, move on
            # if not satisfied, try multiplying the fisher step by
            # 2^((-1)^(cnt+1)*cnt)
            cnt <- 1
            step0 <- step
            while( !wolfe_check(likobj0,likobj,logparms,-step,TRUE) ){
                mult <- 0.5^(cnt)
                step <- mult*step0
                likobj <- likfun(link(logparms - step))
                cnt <- cnt+1
                if(cnt==4){ step0 <- 0.5^(-cnt)*grad/max(info); }
                if(cnt==10){ break }
            } 
        #}
    
        
        # transform grad and info
        loglik <- likobj$loglik        
        grad <- -c(likobj$grad)*dlink(logparms)
        info <- likobj$info*outer(dlink(logparms),dlink(logparms))
        # update parameters
        logparms <- logparms - step

        # print some stuff out
        if(!silent){
            cat(paste0("Iter ",j,": \n"))
            cat("pars = ",  paste0(round(link(logparms),8)), "  \n" )
            cat(paste0("loglik = ", round(loglik,6),         "  \n"))
            cat("grad = ")
            cat(as.character(round(grad,6)))
            cat("\n\n")
        }

        # if gradient is small, then break and return the results        
        if( max(abs(grad)) < 5e-3 ){
            break
        }
    }

    # collect information to return
    betahatinfo <- likfun(link(logparms))        
    betahat <- as.vector(betahatinfo$betahat)
    betacov <- solve(betahatinfo$betainfo)
    sebeta <- sqrt(diag(betacov))
    tbeta <- betahat/sebeta
    #df <- length(y) - ncol(X) # assumes X is full rank!
    #pbeta <- 2*( 1 - pt( abs(tbeta), df ) )
    
    ret <- list(
        covparms = link(logparms), 
        betahat = betahat, 
        sebeta = sebeta,
        betacov = betacov,
        tbeta = tbeta,
        #pbeta = pbeta,
        loglik = loglik
    )
    return(ret)
}





fisher_scoring4 <- function( likfun, start_parms, link, dlink, invlink, silent = FALSE ){
    
    
    maxit <- 50
    
    wolfe_check <- function(likobj0,likobj1,logparms,step,both){
        c1 <- 1e-4
        c2 <- 0.9
        ll0 <- -likobj0$loglik
        gr0 <- -c(likobj0$grad)*dlink(logparms)
        ll1 <- -likobj1$loglik
        gr1 <- -c(likobj1$grad)*dlink(invlink(link(logparms + step)))
        if(!both){
            satfd <- ll1 <= ll0 + c1*crossprod(step,gr0)
        } else {
            satfd <- ll1 <= ll0 + c1*crossprod(step,gr0) &&
                 -crossprod(step,gr1) <= -c2*crossprod(step,gr0)
        }
        return(satfd)
    }
    
    # evaluate function at initial values
    logparms <- invlink(start_parms)
    likobj <- likfun(link(logparms))
    
    if( !test_likelihood_object(likobj,logparms,dlink) ){
        logparms <- invlink( link( 0.1*logparms ) )
        likobj <- likfun(link(logparms))
    }
    
    # transform grad and info
    loglik <- likobj$loglik        
    grad <- -c(likobj$grad)*dlink(logparms)
    info <- likobj$info*outer(dlink(logparms),dlink(logparms))

    # print some stuff out
    if(!silent){
        cat(paste0("Iter ",0,": \n"))
        cat("pars = ",  paste0(round(link(logparms),8)), "  \n" )
        cat(paste0("loglik = ", round(loglik,6),         "  \n"))
        cat("grad = ")
        cat(as.character(round(grad,3)))
        cat("\n\n")
    }
    
    for(j in 1:maxit){
        
        likobj0 <- likobj
        
        # if condition number large, then regularize
        tol <- 1e-8
        if (condition_number(info) > 1 / tol) {
            if (!silent) cat("&&\n")
            eiginfo <- eigen(info)
            whichsmall <- which( eiginfo$values / max(eiginfo$values) < tol )
            eiginfo$values[whichsmall] <- tol*max(eiginfo$values)
            info <-
                eiginfo$vectors %*% diag(eiginfo$values) %*% t(eiginfo$vectors)
        }

        # calculate fisher step and criterion at new parameter
        step <- -solve(info, grad)
        if (mean(step ^ 2) > 1) {
            if(!silent) cat("@@\n")
            step <- step / sqrt(mean(step ^ 2))
        }
        # make sure the step keeps you inside the box
        while(sum(abs( invlink(link(logparms+step)) - logparms-step )) > 1e-9){
            if(!silent) cat("## outside box \n")
            if(!silent) cat(round(invlink(link(logparms+step)),6),"\n")
            if(!silent) cat(round(logparms+step,6),"\n")
            step <- 0.8*step
        }
        newlogparms <- logparms + step
        #newlogparms <- invlink(link(logparms+step))
        if(!silent) cat("^^", round(step,6), "   ", round(newlogparms,6),"\n")
        likobj <- likfun(link(newlogparms))
        
        # check for Inf, NA, or NaN
        cnt <- 1
        while (!test_likelihood_object(likobj, newlogparms, dlink)) {
            if (!silent) cat("inf or na or nan in likobj")
            step <- 0.5 * step
            #newlogparms <- invlink(link(logparms+step))
            newlogparms <- logparms + step
            if(!silent) cat("&&", round(step,6), "   ", round(newlogparms,6),"\n")
            likobj <- likfun(link(newlogparms))
            if (cnt == 8) { stop("could not avoid inf, na or nan") }
        }
        
        # Check the wolfe conditions
        # if not satisfied, shrink fisher step
        # after some interations, switch to gradient and shrink that
        cnt <- 1
        #step0 <- newlogparms - logparms
        # try simply making max cnt smaller,
        # so that the parameters get kicked to a new spot
        # regardless of whether its a better spot
        no_decrease <- FALSE
        both <- TRUE
        while (!wolfe_check(likobj0,likobj,logparms,newlogparms-logparms,both)){
            mult <- 0.5 
            step <- mult * step
            if (cnt == 6) { # switch to gradient
                if(!silent) cat("\n") 
                step <- - grad / sqrt(mean(grad^2)) 
                both <- FALSE
                # maybe we should make sure that norm of step of smaller than 1
                while(sum(abs( invlink(link(logparms+step)) - logparms-step )) > 1e-9){
                    if(!silent){
                    cat("## outside box \n")
                    cat(round(invlink(link(logparms+step)),6),"\n")
                    cat(round(logparms+step,6),"\n")
                    }
                    step <- 0.5*step
                }
            }
            if (cnt == 20) { no_decrease <- TRUE; break }  # maybe we should throw error here?
            newlogparms <- logparms + step
            #newlogparms <- invlink(link(logparms+step))
            if(!silent) cat("**", round(step,6), "   ", round(link(newlogparms),6),"\n")
            likobj <- likfun(link(newlogparms))
            cnt <- cnt + 1
        }

        # transform grad and info
        #logparms <- invlink(link(logparms + step))
        logparms <- logparms + step
        loglik <- likobj$loglik        
        grad <- -c(likobj$grad)*dlink(logparms)
        info <- likobj$info*outer(dlink(logparms),dlink(logparms))
        # update parameters

        # print some stuff out
        if(!silent){
            cat(paste0("Iter ",j,": \n"))
            cat("pars = ",  paste0(round(link(logparms),6)), "  \n" )
            cat(paste0("loglik = ", round(loglik,6),         "  \n"))
            cat("grad = ")
            cat(as.character(round(grad,6)))
            cat("\n\n")
        }

        # if gradient is small, then break and return the results        
        if( max(abs(grad)) < 5e-3 || no_decrease ){
            break
        }
    }

    # collect information to return
    betahatinfo <- likfun(link(logparms))        
    betahat <- as.vector(betahatinfo$betahat)
    betacov <- solve(betahatinfo$betainfo)
    sebeta <- sqrt(diag(betacov))
    tbeta <- betahat/sebeta
    #df <- length(y) - ncol(X) # assumes X is full rank!
    #pbeta <- 2*( 1 - pt( abs(tbeta), df ) )
    
    ret <- list(
        covparms = link(logparms), 
        betahat = betahat, 
        sebeta = sebeta,
        betacov = betacov,
        tbeta = tbeta,
        #pbeta = pbeta,
        loglik = loglik,
        no_decrease = no_decrease
    )
    return(ret)
}






# put link functions in fit_model, rather than here
fisher_scoring <- function( likfun, start_parms, link, silent = FALSE ){
    
    # link functions passed for printing purposes
    maxit <- 40
    
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
    
    if( !test_likelihood_object5(likobj) ){
        logparms <- 0.1*logparms
        likobj <- likfun(logparms)
    }
    
    # transform grad and info
    loglik <- likobj$loglik        
    grad <- likobj$grad
    info <- likobj$info
    # test this for regularity?
    diag(info) <- diag(info) + 0.5*min(diag(info))

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
        
        # if condition number large, then regularize
        tol <- 1e-8
        if (condition_number(info) > 1 / tol) {
            if (!silent) cat("Cond # of info matrix > 1/tol \n")
            eiginfo <- eigen(info)
            whichsmall <- which( eiginfo$values / max(eiginfo$values) < tol )
            eiginfo$values[whichsmall] <- tol*max(eiginfo$values)
            info <-
                eiginfo$vectors %*% diag(eiginfo$values) %*% t(eiginfo$vectors)
        }

        # calculate fisher step and criterion at new parameter
        step <- - solve(info, grad)
        if (mean(step^2) > 1) {
            if(!silent) cat("@@\n")
            step <- step/sqrt(mean(step^2))
        }
        newlogparms <- logparms + step
        #if(!silent) cat("^^", round(step,4), "   ", round(link(newlogparms),4),"\n")
        likobj <- likfun(newlogparms)
        
        # check for Inf, NA, or NaN
        cnt <- 1
        while (!test_likelihood_object5(likobj)) {
            if (!silent) cat("inf or na or nan in likobj")
            step <- 0.5 * step
            newlogparms <- logparms + step
            #if(!silent) cat("!!", round(step,4), "   ", round(link(newlogparms),4),"\n")
            likobj <- likfun(newlogparms)
            if (cnt == 10) { stop("could not avoid inf, na or nan") }
        }
        
        # Check the wolfe conditions
        # if not satisfied, shrink fisher step
        # after some interations, switch to gradient and shrink that
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
            #if(!silent) cat("**", round(step,4), "   ", round(link(newlogparms),4),"\n")
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
        if( abs(stepgrad) < 1e-4 || no_decrease ){
            break
        }
    }

    # collect information to return
    betahatinfo <- likobj        
    betahat <- as.vector(betahatinfo$betahat)
    betacov <- solve(betahatinfo$betainfo)
    sebeta <- sqrt(diag(betacov))
    tbeta <- betahat/sebeta
    #df <- length(y) - ncol(X) # assumes X is full rank!
    #pbeta <- 2*( 1 - pt( abs(tbeta), df ) )
    
    ret <- list(
        covparms = link(logparms), 
        betahat = betahat, 
        sebeta = sebeta,
        betacov = betacov,
        tbeta = tbeta,
        #pbeta = pbeta,
        loglik = loglik,
        no_decrease = no_decrease
    )
    return(ret)
}






