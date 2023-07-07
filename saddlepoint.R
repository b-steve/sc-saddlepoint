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
    n.traps <- nrow(traps)
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

nlog.saddlepoint.density.tmb <- function(par, C, m, traps, version = "original", trace = TRUE){
    EN <- par[1]
    lambda0 <- par[2]
    sigma <- par[3]
    nonzero <- ifelse(C != 0, 1, 0)
    n.traps <- nrow(traps)
    ## Data for TMB.
    data <- list(c = C, traps = traps, m = m, EN = EN, lambda0 = lambda0,
                 sigma = sigma, nonzero = nonzero, is_kgf = 0)
    pars <- list(t = rep(0.1, n.traps)[nonzero == 1])
    ## Making function to obtain t-hat.
    obj.solve <- MakeADFun(data, pars, DLL = "saddlepoint", silent = TRUE)
    data$is_kgf <- 1
    ## Making function for the KGF.
    obj.kgf <- MakeADFun(data, pars, DLL = "saddlepoint", silent = TRUE)
    ## Finding optimal t.
    opt <- nlminb(rep(0.1, n.traps)[nonzero == 1], obj.solve$fn, obj.solve$gr)
    t.hat <- opt$par
    ## KGF first derivatives.
    k.prime <- obj.kgf$gr(t.hat)
    ## KGF Hessian.
    kgf.hess <- optimHess(t.hat, obj.kgf$fn, obj.kgf$gr)
    dim.hess <- dim(kgf.hess)[1]
    ## Negative log-likelihood.
    if (version == "original"){
        out <- -(opt$objective - 0.5*determinant(2*pi*kgf.hess, logarithm = TRUE)$modulus)
    } else if (version == "gaussian"){
        out <- -(opt$objective - 0.5*determinant(diag(dim.hess) + 2*pi*kgf.hess, logarithm = TRUE)$modulus)
    } else {
        if (version == "original2"){
            g <- function(lambda) -0.5*log(2*pi*lambda)
        } else if (version == "gaussian2"){
            g <- function(lambda) -0.5*log(1 + 2*pi*lambda)
        } else if (version == "bessel"){
            g <- function(lambda) log(besselI(x=lambda, nu=0, expon.scaled = TRUE))
        } else if (version == "poisson"){
            g <- function(lambda) -lambda + lambda*log(lambda) - lgamma(1+lambda)
        }
        out <- -(opt$objective + sum(g(eigen(kgf.hess)$values)))
    }
    if (trace){
        cat("par:", par, "LL:", -out, "\n")
    }
    out
}

