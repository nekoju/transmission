library(reshape2)
library(dplyr)
library(abc)
library(msr)
source("functions.r")


## parameters
npop <- 30
npopSamp <- 10
nchrom <- 10
samp <- 100
tau <- 1
rho <- 1
nh <- 1000
u <- 6 * 10 ^ (-8) * 1500
m <- 0.005
thetaH <- 2 * nh * u
M <- 2 * nh * m
const <- rho / (tau ^ 2 * (3 - 2 * tau) * (2 - rho) + rho)
tauPriors <- c(1 , 1)
rhoPriors <- c(10 , 10)

## ms calls
msTest <- call_ms(
    npop * nchrom , samp , 
    t = thetaH * const ,
    I = c(npop , rep(nchrom , times = npop) , M * const) 
    ) %>%
parse_ms()

## data tidying
snps <- as.data.frame(
    melt(
        msTest$gametes[!sapply(msTest$gametes , is.null)]
        , varnames = c("chrom" , "pos") , value.name = "gt"
        )
    )
names(snps)[4] <- "samp"
snps$pop <- vapply(
    snps$chrom , 
    function (x) getPop(x , matrix(1:(nchrom * npop) , nrow = npop , byrow = TRUE) , 1:npop) , factor(1)
    )
snps <- snps[c(1 , 2 , 4 , 5 , 3)]
snpsSubset <- filter(snps , pop %in% sample(1:npop , npopSamp))

ht <- group_by(snpsSubset , pos , samp)  %>%
    summarize(
        pBar = sum(gt) / (nchrom * npopSamp) ,
        qBar = 1 - pBar ,
        ht = ((nchrom * npopSamp / (nchrom * npopSamp - 1)) * 2 * pBar * qBar) ,
    ) %>%
    group_by(samp) %>%
    summarize(ht = mean(ht) , nseg = length(unique(pos)))
hs <- group_by(snpsSubset , pos , samp , pop) %>%
    summarize(p = sum(gt) / nchrom , q = 1 - p) %>%
    group_by(pos , samp) %>%
    summarize(hs = (nchrom / (nchrom - 1)) * sum(2 * p * q) / npopSamp) %>%
    group_by(samp) %>%
    summarize(hs = mean(hs) , nseg = max(pos))
fst <- merge(ht , hs)[ht$ht != 0, ]
fst <- mutate(fst , fst = 1 - hs / ht , fstStd = fst / (1 - hs))
fstWt <- weighted.mean(fst$fst , fst$nseg)
fstWtVar <- sum(fst$nseg * (fst$fst - fstWt) ^ 2) / (sum(fst$nseg) - 1)
htBar <- weighted.mean(fst$ht , fst$nseg) 
hsBar <- weighted.mean(fst$hs , fst$nseg)
fstBar <- 1 - hsBar / htBar
fstBarStd <- fstBar / (1 - hsBar)

# nTilde <- nchrom
# p <- group_by(snpsSubset , pos , samp , pop) %>%
#     summarize(p = sum(gt) / nchrom , q = 1 - p)
# p <- melt(p , measure.vars = c("p" , "q") , variable.name = "allele" , value.name = "freq")
# h0 <- group_by(p , pos , samp) %>%
#     summarize(h0 = 1 - sum(freq ^ 2 / npopSamp))
# hs <- group_by(p , pos , samp , allele) %>%
#     summarize(xSquaredBar = sum(freq ^ 2) / npopSamp) %>%
#     group_by(pos , samp) %>%
#     summarize(xSquaredBarSum = sum(xSquaredBar))
# ht <- group_by(p , pos , samp , allele) %>%
#     summarize(xBarSquared = (sum(freq) / npopSamp) ^ 2) %>%
#     group_by(pos , samp) %>%
#     summarize(xBarSquaredSum = sum(xBarSquared))

# fstClassic <- filter(p , allele == "p" & !(freq %in% c(0 , 1))) %>%
#     group_by(samp , pos) %>%
#     summarize(
#         pSquaredBar = mean(freq ^2) ,
#         pBar = mean(freq) ,
#         pBarSquared = pBar ^ 2 ,
#         fst = (pSquaredBar - pBarSquared) / (pBar * (1 - pBar))
#         )

    # fstUnbiased <- merge(merge(h0 , hs) , ht)
    # fstUnbiased <- mutate(
    #     fstUnbiased ,
    #     hs = nTilde / (nTilde - 1) * (1 - xSquaredBarSum - h0 / (2 * nTilde)) ,
    #     ht = 1 - xBarSquaredSum + hs / (nTilde * npopSamp) -h0 / (2 * nTilde * npopSamp) , 
    #     fst = 1 - hs / ht
    #     )
    # fstUnbiased <- fstUnbiased[!is.na(fstUnbiased$fst), ]


    # thetaS must be divided by number of populations because theta hat is for the
    # metapopulation while thetaH is for the subpopulations
    thetaS <- group_by(snps , samp) %>%
        summarize(K = length(unique(pos)) , a = sum(1 / 1 : (nchrom * npop - 1)) , thetaS = K / (npop * a)) %>%
        group_by(samp) %>%
        summarize(K = mean(K) , a = mean(a) , thetaS = mean(thetaS))

    dist <- group_by(snps , pop , pos) %>%
        summarize(p = sum(gt) / length(gt) , q = 1 - p)

    nsim <- 1000000
    sims <- data.frame(
        tau = rbeta(nsim , tauPriors[1] , tauPriors[2]) , 
        rho = rbeta(nsim , rhoPriors[1] , rhoPriors[2]) * (2 - 0) + 0
        )
    sims <- mutate(sims , fst = getFst(tau , rho , M) , thetaS = getThetaS(tau , rho , thetaH))

    abcTest <- abc(
        target = fstBar,
        param = sims[ ,c(1 , 2)] , 
        sumstat = sims[,c(3)] , 
        tol = 0.010 , method = "rejection"
        )


    par(mfrow = c(2 , 2))
    contPointsTau <- seq(0 , 1 , by = 0.01)
    contPointsRho <- seq(0 , 2 , by = 0.02)
    contPointsF <- outer(contPointsTau , contPointsRho , getFst , M = M)
    contPointsF[is.na(contPointsF)] <- 0
    plot(rho ~ tau , sims[sample(1:nrow(sims) , 10000), ] , xlim = c(0 , 1) , ylim = c(0 , 2))
    with(abcTest , points(unadj.values[,1] , unadj.values[,2] , col = "red"))
    contour(
        contPointsTau , contPointsRho , contPointsF ,
        levels = c(fstBar , getFst(tau , rho , M)) ,
        xlab = "tau" , ylab = "rho" ,
        add = TRUE , labels = c("fst" , "fst-pred") , col = "green" , lty = c(1 , 2)
        )
    # contour(
    #     contPointsTau , contPointsRho , contPointsF ,
    #     levels = fst$fst[sample(1:length(fst$fst) , 10)] , 
    #     xlab = "tau" , ylab = "rho"
    #     )

    # contour(
    #     contPointsTau , contPointsRho , contPointsT ,
    #     levels = c(harMean(thetaS$thetaS) , getThetaS(tau , rho , thetaH))  ,
    #     xlab = "" , ylab = "" , 
    #     add = TRUE , labels = "theta" , col = "blue" , lty = c(1 , 2)
    # )
    # abline(v = c(qbeta(0.50 , tauPriors[1] , tauPriors[2]) , median(abcTest$unadj.values[,1])) , lty = c(1 , 2) , col = "purple")
    # times 2 for nonstandard beta distribution
    # abline(h = c(qbeta(0.50 , rhoPriors[1] , rhoPriors[2]) * 2 , median(abcTest$unadj.values[,2])) , lty = c(1 , 2) , col = "purple")
    # abline(v = tau , h = rho , col = "blue")
    # plot(thetaS ~ fst , sims)

    plot(
        seq(0 , 1 , by = 0.05) , dbeta(seq(0 , 1 , by = 0.05) , tauPriors[1] , tauPriors[2]) ,
        type = 'l' , col = "red" , 
        xlab = "tau" , ylab = "density" , 
        ylim = c(
            0 ,
            max(
                c(
                    max(density(abcTest$unadj.values[,1])$y) ,
                    max(dbeta(seq(0 , 1 , by = 0.05) , tauPriors[1] , tauPriors[2]))
                    )
                )
            )
        )
    lines(density(abcTest$unadj.values[,1]))
    plot(
        seq(0 , 2 , by = 0.05) , dbeta(seq(0 , 1 , by = 0.025) ,rhoPriors[1] , rhoPriors[2]) , 
        type = 'l' , col = "red" , 
        xlab = "rho" , ylab = "density" , 
        ylim = c(
            0 ,
            max(
                c(
                    max(density(abcTest$unadj.values[,2])$y) ,
                    max(dbeta(seq(0 , 1 , by = 0.05) , rhoPriors[1] , rhoPriors[2]))
                    )
                )
            )
        )
    lines(density(abcTest$unadj.values[,2]))
    par(mfrow = c(1 , 1))

    priorSims <- read.csv("~/Data/MixedTransmission/msOut.fst" , header = TRUE) 
    priorParams <- read.csv("~/Data/MixedTransmission/msHyperParams.sh" , header = TRUE)
    priors <- data.frame(cbind(priorSims , priorParams))
    simData <- cbind(priorParams , priorSims)
    simData <- mutate(
        simData ,
        fstScaled = scale(mean.fst , center = fstWt) ,
        sdScaled = scale(sqrt(var.fst) , center = sqrt(fstWtVar)) 
        )

    abcMs <- abc(
        target = c(0 , 0) ,
        param = simData[ , c("tau" , "rho")] ,
        sumstat = simData[ , c("fstScaled" , "sdScaled")] ,
        method = "rejection" , 
        tol = 0.1
        )
