library(spatstat.geom)
library(pracma)
library(nloptr)

## KGF and MGF functions.
kgf.x.given.s <- function(t, s, traps, lambda0, sigma, nz){
    d <- sqrt(rowSums(cbind(s[1] - traps[, 1], s[2] - traps[, 2])^2))
    lambdas <- lambda0*exp(-(d^2)/(2*sigma^2))
    kgf.nz <- sum(kgf.pois(t, lambdas[nz]))
    kgf.z <- sum(dpois(0, lambdas[-nz], log = TRUE))
    kgf.nz + kgf.z
}

mgf.x.given.s <- function(t, s, traps, lambda0, sigma, nz){
    exp(kgf.x.given.s(t, s, traps, lambda0, sigma, nz))
}

mgf.x <- function(t, m, traps, lambda0, sigma, nz){
    m.area <- (m[1, 1] - m[2, 1])^2
    n.mask <- nrow(m)
    mgf.s <- numeric(n.mask)
    for (i in 1:nrow(m)){
        mgf.s[i] <- mgf.x.given.s(t, m[i, ], traps, lambda0, sigma, nz)
    }
     m.area <- (m[1, 1] - m[2, 1])^2
    m.area*sum(mgf.s)/(m.area*nrow(m))
}
kgf.x <- function(t, m, traps, lambda0, sigma, nz){
    log(mgf.x(t, m, traps, lambda0, sigma, nz))
}

kgf.c <- function(t, m, traps, lambda0, sigma, EN, nz){
    kgf.pois(kgf.x(t, m, traps, lambda0, sigma, nz), EN)
}
mgf.c <- function(t, m, traps, lambda0, sigma, EN, nz){
    exp(kgf.c(t, m, traps, lambda0, sigma, EN, nz))
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
    nz <- which(C > 0)
    ## Finding optimal t.
    f <- function(t, m, traps, lambda0, sigma, EN, C, nz){
        kgf.c(t, m, traps, lambda0, sigma, EN, nz) - (t %*% C[nz])
    }
    opt <- nlminb(rep(0.1, n.traps)[nz], f, m = m, traps = traps,
                  lambda0 = lambda0, sigma = sigma, EN = EN, C = C, nz = nz)
    t.hat <- opt$par
    ## KGF first derivatives.
    k.prime <- jacobian(kgf.c, t.hat, m = m, traps = traps, lambda0 = lambda0, sigma = sigma, EN = EN, nz = nz)
    ## KGF Hessian.
    kgf.hess <- hessian(kgf.c, t.hat, m = m, traps = traps, lambda0 = lambda0, sigma = sigma, EN = EN, nz = nz)
    ## Negative log-likelihood.
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
plot(traps, type = "n", xlim = xlim, ylim = ylim)
text(traps, labels = C)
points(S)

## Single log-likelihood evaluation.
nlog.saddlepoint.density(c(20, 15, 0.25), C = C, m = m, traps = traps)

## Fitting a model.
fit <- bobyqa(c(4, 10, 0.50), nlog.saddlepoint.density, lower = c(0, 0, 0) + 1e-5, upper = c(Inf, Inf, Inf), C = C, m = m, traps = traps)

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
