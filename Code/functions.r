## function to get fst given host Ne * m and transmission params
#  based on approximation
getFst <- function(tau , rho , M) {
    (tau^2 * (3 - 2 * tau) * (2 - rho) + rho) /
        (tau^2 * (3 - 2 * tau) * (2 - rho) + rho * (M + 1))
}

## function to get symbiont theta from host theta estimator
getThetaS <- function (tau , rho , thetaH) {
    (thetaH * rho) / (tau^2 * (3 - 2 * tau) * (2 - rho) + rho)
}

## function to get population from sample number in ms data
getPop <- function (x , ranges , pops) {
    for (i in 1:nrow(ranges)) {
        if (x %in% ranges[i, ]) {
            return(pops[i])
        } else {
            NA
        }
    } 
}

## returns harmonic mean
harMean <- function (x) {
    1 / mean(1 / x)
}

msCall <- function (tau , rho , thetaH , M , npop , nchrom , samp) {
    const <- rho / (tau ^ 2 * (3 - 2 * tau) * (2 - rho) + rho)
    paste(
        c(
            "ms" , 
            as.character(nchrom * npop) ,
            as.character(samp) , '-t' ,
            as.character(thetaH * const) ,
            '-I' ,
            as.character(npop) , 
            as.character(rep(nchrom , times = npop)) ,
            as.character(M * const) 
            ) ,
        collapse = ' '
        )
}

