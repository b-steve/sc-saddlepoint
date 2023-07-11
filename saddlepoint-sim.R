library(secr)
source("saddlepoint.R")
source("quasifunctions.R")

## Simulation.
set.seed(1234)
n.sims <- 100
traps <- as.matrix(expand.grid(1:10, 1:10))
n.traps <- nrow(traps)
N <- EN <- 50
lambda0 <- 3
sigma <- 1
xlim <- range(traps[,1]) + c(-4, 4)
ylim <- range(traps[,2]) + c(-4, 4)
m <- as.matrix(expand.grid(x = seq(xlim[1], xlim[2], 0.25), y = seq(ylim[1], ylim[2], 0.25)  ))

## For David's function.
simtraps <- read.traps(data = data.frame(x = traps[, 1], y = traps[, 2]), type = "proximity")
simmask <- read.mask(data = m)
dists <- edist(simtraps,simmask)
dists.sq <- dists^2
a <- summary(simmask)$cellarea*10^4
A <- length(simmask$x)*a
# Compile Nimble code
calcEns <- nimEns(simtraps, simmask)
calcEns.c <- compileNimble(calcEns)

## Matrices for estimates.
pars1 <- pars2 <- pars3 <- pars4 <- pars5 <- pars6 <- pars7 <- matrix(0, nrow = n.sims, ncol = 3)
for (i in 1:n.sims){
    cat(i, "\n")
    S <- cbind(runif(N, xlim[1], xlim[2]), runif(N, xlim[1], xlim[2]))
    dists <- crossdist(S[, 1], S[, 2], traps[, 1], traps[, 2])
    lambdas <- lambda0*exp(-dists^2/(2*sigma^2))
    capt <- matrix(rpois(N*n.traps, lambdas), nrow = N, ncol = n.traps)
    C <- colSums(capt)
    fit1 <- bobyqa(c(30, 10, 0.50), nlog.saddlepoint.density.tmb, lower = c(0, 0, 0) + 1e-5,
                   upper = c(Inf, Inf, Inf), C = C, m = m, traps = traps, version = "original", trace = FALSE)
    fit2 <- bobyqa(c(30, 10, 0.50), nlog.saddlepoint.density.tmb, lower = c(0, 0, 0) + 1e-5,
                   upper = c(Inf, Inf, Inf), C = C, m = m, traps = traps, version = "original2", trace = FALSE)
    fit3 <- bobyqa(c(30, 10, 0.50), nlog.saddlepoint.density.tmb, lower = c(0, 0, 0) + 1e-5,
                   upper = c(Inf, Inf, Inf), C = C, m = m, traps = traps, version = "gaussian", trace = FALSE)
    fit4 <- bobyqa(c(30, 10, 0.50), nlog.saddlepoint.density.tmb, lower = c(0, 0, 0) + 1e-5,
                   upper = c(Inf, Inf, Inf), C = C, m = m, traps = traps, version = "gaussian2", trace = FALSE)
    fit5 <- bobyqa(c(30, 10, 0.50), nlog.saddlepoint.density.tmb, lower = c(0, 0, 0) + 1e-5,
                   upper = c(Inf, Inf, Inf), C = C, m = m, traps = traps, version = "bessel", trace = FALSE)
    fit6 <- bobyqa(c(30, 10, 0.50), nlog.saddlepoint.density.tmb, lower = c(0, 0, 0) + 1e-5,
                   upper = c(Inf, Inf, Inf), C = C, m = m, traps = traps, version = "poisson", trace = FALSE)
    fit7 <- optim(log(c(D = 30/A, lambda0 = 10, sigmasq = 0.5^2)),
                        quasifn, y = C, dists2 = dists.sq, A = A, a = a)
    pars1[i, ] <- fit1$par
    pars2[i, ] <- fit2$par
    pars3[i, ] <- fit3$par
    pars4[i, ] <- fit4$par
    pars5[i, ] <- fit5$par
    pars6[i, ] <- fit6$par
    pars7[i, ] <- c(exp(fit7$par[1])*A, exp(fit7$par[2]), sqrt(exp(fit7$par[3])))
    ## Live boxplot.
    par(mfrow = c(1, 3))
    true.pars <- c(EN, lambda0, sigma)
    par.names <- c("N", "lambda0", "sigma")
    for (j in 1:3){
        boxplot(pars2[1:i, j], pars4[1:i, j], pars5[1:i, j], pars6[1:i, j],
                names = c("Original", "Gaussian", "Bessel", "Poisson"),
                main = par.names[j])
    }
}
