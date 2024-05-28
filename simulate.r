library(secr)
source("quasifunctions.R")
source("saddlepoint.R")

## Paul's unknown ID SCR data simulation function
simSCR <- function(N = 50, sigma = 0.5, lambda = 0.5, StudyPeriod = 25, traps, xlim, ylim)
{
    locs <- cbind(x = runif(N, xlim[1], xlim[2]), 
                  y = runif(N, ylim[1], ylim[2]))
    J <- nrow(traps)
    capthist <- NULL
	ids <- NULL
    for(i in 1:N)
    {
        d2 <- (locs[i,1] - traps[,1])^2 + (locs[i,2] - traps[,2])^2
        Hkj <- lambda*exp(-d2/(2*sigma^2))*StudyPeriod
		nkj <- rpois(J, Hkj)
		if(sum(nkj) > 0) capthist <- rbind(capthist, nkj)
    }
	return(colSums(capthist))
}


## Paul's simulation setup code.
traps <- expand.grid(x = 1:10,y = 1:10)
N <- 50
sigma <- 1
lambda <- 3
xlim <- range(traps[,1]) + c(-4, 4)
ylim <- range(traps[,2]) + c(-4, 4)
mask <- expand.grid(x = seq(xlim[1], xlim[2], 0.25), y = seq(ylim[1], ylim[2], 0.25)  )
nmask <- nrow(mask)
J <- nrow(traps)
A <- 0.25^2
area <- nrow(mask) * A
D <- N/area
d2mask <- as.matrix(t(apply(mask, 1, FUN = function(x){(x[1]-traps[,1])^2 + (x[2]-traps[,2])^2})))


# Turn traps and mask into secr objects because the function I have for covariance 
# assumes that they are. Should no doubt change this in due course ...
simtraps = read.traps(data=traps, type="proximity")
simmask = read.mask(data=mask)
plot(simmask, border=0)
plot(simtraps,add=TRUE)

# calculate the things that are  constant over simulations
dists = edist(simtraps,simmask)
dists2 = dists^2
ntraps = nrow(simtraps)
nmask = nrow(simmask)
a = summary(simmask)$cellarea*10^4
A = length(simmask$x)*a
D=N/A

# set up starting parameter values
# (uses sigma^2, not sigma; might be sensible to change that ...)
pars = log(c(D=D, lambda0=lambda, sigmasq=sigma^2))

# Compile Nimble code
calcEns <- nimEns(simtraps, simmask)
calcEns.c <- compileNimble(calcEns)

# Generate some data and check that the objective function evaluates and fit seems to work
# ------------------
set.seed(9)
count <- simSCR(N, sigma, lambda, 1, traps, xlim, ylim)
plot.nu(count,mask=simmask,traps=simtraps,what="both")
# evaluate quasi-likelihood function
quasifn(pars, y=count, dists2, A, a)
# ... and that fit seems to do the right thing
qtime <- system.time({fit=optim(pars, quasifn, y=count, dists2=dists2, A=A, a=a, control=list(trace=5))})

sp.time <- system.time({fit.sp <- bobyqa(c(50, 3, 1), nlog.saddlepoint.density.tmb, lower = c(0, 0, 0) + 1e-5,
                          upper = c(Inf, Inf, Inf), C = count, m = as.matrix(mask), traps = as.matrix(traps),
                          trace = FALSE)})
sp.gauss.time <- system.time({fit.sp.gauss <- bobyqa(c(50, 3, 1), nlog.saddlepoint.density.tmb, lower = c(0, 0, 0) + 1e-5,
                                                     upper = c(Inf, Inf, Inf), C = count, m = as.matrix(mask), traps = as.matrix(traps),
                                                     version = "gaussian", trace = FALSE)})
## Comparison:
exp(fit$par)
fit.sp$par

exp(fit$par)[1]
fit.sp$par[1]/A

# Try simulating:
# --------------
nsim=100
ests = matrix(rep(0,3*nsim),nrow=nsim) # matrix to hold simulation results
ests.sp <- matrix(rep(0,3*nsim),nrow=nsim)
ests.sp.gauss <- matrix(rep(0,3*nsim),nrow=nsim)
set.seed(12321)
stime = date()
time.quasi <- numeric(nsim)
time.sp <- numeric(nsim)
time.sp.gauss <- numeric(nsim)
par(mfrow = c(1, 3))
for(i in 1:nsim) {
    cat(i, "\n")
    count <- simSCR(N, sigma, lambda, 1, traps, xlim, ylim)
    time.quasi[i] <- system.time({fit=optim(pars, quasifn, y=count, dists2=dists2, A=A, a=a, control=list(trace=5))})[3]
    time.sp[i] <- system.time({fit.sp <- bobyqa(c(50, 3, 1), nlog.saddlepoint.density.tmb,
                                                lower = c(0, 0, 0) + 1e-5,
                                                upper = c(Inf, Inf, Inf),
                                                C = count, m = as.matrix(mask),
                                                traps = as.matrix(traps),
                                                trace = FALSE)})[3]
    time.sp.gauss[i] <- system.time({fit.sp.gauss <- bobyqa(c(50, 3, 1), nlog.saddlepoint.density.tmb,
                                                            lower = c(0, 0, 0) + 1e-5,
                                                            upper = c(Inf, Inf, Inf),
                                                            C = count, m = as.matrix(mask),
                                                            traps = as.matrix(traps),
                                                            version = "gaussian",
                                                            trace = FALSE)})[3]
    ests[i,] = exp(fit$par) # save estimates on natural scale
    ests.sp[i, ] <- fit.sp$par
    ests.sp.gauss[i, ] <- fit.sp.gauss$par
    boxplot(ests[1:i, 1], ests.sp[1:i, 1]/A, ests.sp.gauss[1:i, 1]/A); abline(h = D, lty = "dotted")
    boxplot(ests[1:i, 2], ests.sp[1:i, 2], ests.sp.gauss[1:i, 2]); abline(h = lambda, lty = "dotted")
    boxplot(ests[1:i, 3], ests.sp[1:i, 3], ests.sp.gauss[1:i, 3]); abline(h = sigma, lty = "dotted")
}
etime = date()
simest = data.frame(D=ests[,1], lambda0=ests[,2], sigma=sqrt(ests[,3])) # convert sigma^2 to sigma
simestmean = apply(simest,2,mean) # mean estimates for each parameter
truth = c(D=D, lambda0=lambda, sigma=sigma) # true parameters
pcbias = 100*(simestmean-truth)/truth # % bias
pdf(file="QuasiSim.pdf",h=4,w=12)
par(mfrow=c(1,3))
for(i in 1:3) {
  title = paste(names(simestmean)[i],": Bias = ",signif(pcbias[i],1),"%",sep="")
  boxplot(ests[,i],main=title)
  abline(h=truth[i],col="red")
}
dev.off()

quasisim1 = list(truth=truth, simest=simest, seed=12321)
saveRDS(quasisim1,file="QuasiSim1.Rds")


pdf("sc-comparison.pdf")
boxplot(ests[, 1], ests.sp[, 1]/A, main = "Density", names = c("Quasi", "Saddlepoint"))
abline(h = D, lty = "dotted")
points(c(1, 2), c(mean(ests[, 1]), mean(ests.sp[, 1])/A), col = "red", pch = 16)

boxplot(ests[, 2], ests.sp[, 2], main = "Lambda0", names = c("Quasi", "Saddlepoint"))
abline(h = lambda, lty = "dotted")
points(c(1, 2), c(mean(ests[, 2]), mean(ests.sp[, 2])), col = "red", pch = 16)

boxplot(ests[, 3], ests.sp[, 3], main = "Lambda0", names = c("Quasi", "Saddlepoint"))
abline(h = sigma, lty = "dotted")
points(c(1, 2), c(mean(ests[, 3]), mean(ests.sp[, 3])), col = "red", pch = 16)
dev.off()

plot(ests, ests.sp)
