library(spatstat.geom)
library(pracma)
library(nloptr)


library(TMB)
## Compling TMB implementation.
compile("tmb/saddlepoint.cpp")
dyn.load(dynlib("tmb/saddlepoint"))

## KGF and MGF functions.
kgf.x.given.s <- function(t, s, traps, lambda0, sigma, nonzero){
    d <- sqrt(rowSums(cbind(s[1] - traps[, 1], s[2] - traps[, 2])^2))
    lambdas <- lambda0*exp(-(d^2)/(2*sigma^2))
    kgf.nonzero <- sum(kgf.pois(t, lambdas[nonzero == 1]))
    kgf.zero <- sum(dpois(0, lambdas[nonzero == 0], log = TRUE))
    kgf.nonzero + kgf.zero
}

mgf.x.given.s <- function(t, s, traps, lambda0, sigma, nonzero){
    exp(kgf.x.given.s(t, s, traps, lambda0, sigma, nonzero))
}

mgf.x <- function(t, m, traps, lambda0, sigma, nonzero){
    m.area <- (m[1, 1] - m[2, 1])^2
    n.mask <- nrow(m)
    mgf.s <- numeric(n.mask)
    for (i in 1:nrow(m)){
        mgf.s[i] <- mgf.x.given.s(t, m[i, ], traps, lambda0, sigma, nonzero)
    }
    ##m.area*sum(mgf.s)/(m.area*nrow(m))
    sum(mgf.s)/nrow(m)
}
kgf.x <- function(t, m, traps, lambda0, sigma, nonzero){
    log(mgf.x(t, m, traps, lambda0, sigma, nonzero))
}

kgf.c <- function(t, m, traps, lambda0, sigma, EN, nonzero){
    kgf.pois(kgf.x(t, m, traps, lambda0, sigma, nonzero), EN)
}
mgf.c <- function(t, m, traps, lambda0, sigma, EN, nonzero){
    exp(kgf.c(t, m, traps, lambda0, sigma, EN, nonzero))
}

kgf.pois <- function(t, lambda){
    lambda*(exp(t) - 1)
}
mgf.pois <- function(t, lambda){
    exp(kgf.pois(t, lambda))
}

nlog.saddlepoint.density <- function(par, C, m, traps, trace = TRUE){
    EN <- par[1]
    lambda0 <- par[2]
    sigma <- par[3]
    nonzero <- ifelse(C != 0, 1, 0)
    ## Finding optimal t.
    f <- function(t, m, traps, lambda0, sigma, EN, C, nonzero){
        kgf.c(t, m, traps, lambda0, sigma, EN, nonzero) - (t %*% C[nonzero == 1])
    }
    opt <- nlminb(rep(0.1, n.traps)[nonzero == 1], f, m = m, traps = traps,
                  lambda0 = lambda0, sigma = sigma, EN = EN, C = C, nonzero = nonzero)
    t.hat <- opt$par
    ## KGF first derivatives.
    k.prime <- jacobian(kgf.c, t.hat, m = m, traps = traps, lambda0 = lambda0, sigma = sigma, EN = EN, nonzero = nonzero)
    ## KGF Hessian.
    kgf.hess <- hessian(kgf.c, t.hat, m = m, traps = traps, lambda0 = lambda0, sigma = sigma, EN = EN, nonzero = nonzero)
    ## Negative log-likelihood.
    out <- -(opt$objective - log(sqrt(det(2*pi*kgf.hess))))
    if (trace){
        cat("par:", par, "LL:", -out, "\n")
    }
    out
}

nlog.saddlepoint.density.tmb <- function(par, C, m, traps, trace = TRUE){
    EN <- par[1]
    lambda0 <- par[2]
    sigma <- par[3]
    nonzero <- ifelse(C != 0, 1, 0)
    ## Data for TMB.
    data <- list(c = C, traps = traps, m = m, EN = EN, lambda0 = lambda0,
                 sigma = sigma, nonzero = nonzero, is_kgf = 0)
    pars <- list(t = rep(0.1, n.traps)[nonzero == 1])
    ## Making function to obtain t-hat.
    obj.solve <- MakeADFun(data, pars, DLL = "saddlepoint")
    data$is_kgf <- 1
    ## Making function for the KGF.
    obj.kgf <- MakeADFun(data, pars, DLL = "saddlepoint")
    ## Finding optimal t.
    opt <- nlminb(rep(0.1, n.traps)[nonzero == 1], obj.solve$fn, obj.solve$gr)
    t.hat <- opt$par
    ## KGF first derivatives.
    k.prime <- obj.kgf$gr(t.hat)
    ## KGF Hessian.
    kgf.hess <- optimHess(t.hat, obj.kgf$fn, obj.kgf$gr)
    ## ## Negative log-likelihood.
    out <- -(opt$objective - log(sqrt(det(2*pi*kgf.hess))))
    if (trace){
        cat("par:", par, "LL:", -out, "\n")
    }
    out
}
             

## Simulating a single data set.
set.seed(1111)
traps <- as.matrix(expand.grid(1:5, 1:5))
n.traps <- nrow(traps)
EN <- 4
sigma <- 0.5
lambda0 <- 10
xlim <- range(traps[,1]) + c(-0.5, 0.5)
ylim <- range(traps[,2]) + c(-0.5, 0.5)
m <- as.matrix(expand.grid(seq(xlim[1], xlim[2], 0.25), seq(ylim[1], ylim[2], 0.25)))
N <- EN
S <- cbind(runif(N, xlim[1], xlim[2]), runif(N, xlim[1], xlim[2]))
dists <- crossdist(S[, 1], S[, 2], traps[, 1], traps[, 2])
lambdas <- lambda0*exp(-dists^2/(2*sigma^2))
capt <- matrix(rpois(N*n.traps, lambdas), nrow = N, ncol = n.traps)
C <- colSums(capt)

## Plotting data.
#plot(traps, type = "n", xlim = xlim, ylim = ylim)
#text(traps, labels = C)
#points(S)

t <- runif(sum(C > 0), 0, 0.1)
s <- c(3, 3)
nonzero <- ifelse(C != 0, 1, 0)
## Using TMB.
data <- list(c = C, traps = traps, m = m, EN = EN, lambda0 = lambda0, sigma = sigma, nonzero = nonzero, is_kgf = 0)
pars <- list(t = t)
obj.solve <- MakeADFun(data, pars, DLL = "saddlepoint")
data$is_kgf <- 1
obj.kgf <- MakeADFun(data, pars, DLL = "saddlepoint")

f <- function(t, m, traps, lambda0, sigma, EN, C, nonzero){
    kgf.c(t, m, traps, lambda0, sigma, EN, nonzero) - (t %*% C[nonzero == 1])
}
opt <- nlminb(rep(0.1, n.traps)[nonzero == 1], f, m = m, traps = traps,
              lambda0 = lambda0, sigma = sigma, EN = EN, C = C, nonzero = nonzero)
opt
t.hat <- opt$par
## KGF first derivatives.
k.prime <- jacobian(kgf.c, t.hat, m = m, traps = traps, lambda0 = lambda0, sigma = sigma, EN = EN, nonzero = nonzero)
## KGF Hessian.
kgf.hess <- hessian(kgf.c, t.hat, m = m, traps = traps, lambda0 = lambda0, sigma = sigma, EN = EN, nonzero = nonzero)

opt.tmb <- nlminb(rep(0.1, n.traps)[nonzero == 1], obj.solve$fn, obj.solve$gr)
t.hat.tmb <- opt.tmb$par
k.prime.tmb <- obj.kgf$gr(t.hat.tmb)
kgf.hess.tmb <- optimHess(t.hat.tmb, obj.kgf$fn, obj.kgf$gr)

f(t, m, traps, lambda0, sigma, EN, C, nonzero)
obj.solve$fn(t)

debug(nlog.saddlepoint.density)
debug(nlog.saddlepoint.density.tmb)

## Single log-likelihood evaluation.
nlog.saddlepoint.density(c(20, 15, 0.25), C = C, m = m, traps = traps) # 40.43862
nlog.saddlepoint.density.tmb(c(20, 15, 0.25), C = C, m = m, traps = traps)

## Fitting a model.
fit <- bobyqa(c(4, 10, 0.50), nlog.saddlepoint.density, lower = c(0, 0, 0) + 1e-5, upper = c(Inf, Inf, Inf), C = C, m = m, traps = traps) # par: 1.589707 30.61105 0.4734628 LL: -34.2657 

## Simulation.
n.sims <- 100
traps <- as.matrix(expand.grid(1:5, 1:5))
n.traps <- nrow(traps)
N <- EN <- 4
sigma <- 0.5
lambda0 <- 10
xlim <- range(traps[,1]) + c(-0.5, 0.5)
ylim <- range(traps[,2]) + c(-0.5, 0.5)
m <- as.matrix(expand.grid(seq(xlim[1], xlim[2], 0.25), seq(ylim[1], ylim[2], 0.25)))
pars <- matrix(0, nrow = n.sims, ncol = 3)
for (i in 1:n.sims){
    cat(i, "\n")
    S <- cbind(runif(N, xlim[1], xlim[2]), runif(N, xlim[1], xlim[2]))
    dists <- crossdist(S[, 1], S[, 2], traps[, 1], traps[, 2])
    lambdas <- lambda0*exp(-dists^2/(2*sigma^2))
    capt <- matrix(rpois(N*n.traps, lambdas), nrow = N, ncol = n.traps)
    C <- colSums(capt)
    fit <- bobyqa(c(4, 10, 0.50), nlog.saddlepoint.density, lower = c(0, 0, 0) + 1e-5,
                  upper = c(Inf, Inf, Inf), C = C, m = m, traps = traps, trace = FALSE)
    pars[i, ] <- fit$par
    par(mfrow = c(3, 1))
    boxplot(pars[1:i, 1]); abline(h = EN, lty = "dotted", main = "EN")
    boxplot(pars[1:i, 2]); abline(h = lambda0, lty = "dotted", main = "EN")
    boxplot(pars[1:i, 3]); abline(h = sigma, lty = "dotted", main = "EN")
}
