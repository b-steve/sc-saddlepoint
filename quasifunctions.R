library(nimble)
                                        # function to calcuate E(n) using numerical approximation with mask
# If by=="trap", calculates E(n) per trap, else per trap and mask cell
mu.En.hn = function(pars,dists2,a,by="trap") {
  # half-normal encounter function
  hn = function(d2,l0,sig2) return(l0*exp(-0.5*d2/sig2))
  # unpack parameters
  D = exp(pars[1])
  l0 = exp(pars[2])
  sig2 = exp(pars[3])
  # calculate the encounter function at each mask point (column) for each trap (row)
  Ens = dists2*0
  for(i in 1:dim(Ens)[1]) Ens[i,] = hn(dists2[i,],l0,sig2)*a
  if(by=="trap") return(D*apply(Ens,1,sum))
  else return(Ens)
}

# Function to calculate variance-covariance matrix
# Input is 
# Ens: Expected number of counts at each trap (row) for AC at each mask point (column),
# mu_N: Expected number of animals within mask
make.cov = function(Ens, mu_N) {
  ndet = dim(Ens)[1]
  E.lam = apply(Ens,1,mean)
  E.lamii = apply(Ens^2,1,mean)
  cov = matrix(rep(0,ndet^2),nrow=ndet)
  diag(cov) = mu_N*(E.lamii + E.lam) # variances
  for(i in 1:ndet) 
    for(j in 1:ndet) 
      cov[j,i] = mu_N*mean(Ens[i,]*Ens[j,]) # covariances
  return(cov)
}

# Nimble code to calculate derivatives and Ens
nimEns <- nimbleFunction(
  setup = function(traps = double(2), mask = double(2)){
    dists = edist(traps,mask)
    dists2 = dists^2
    ntraps = nrow(traps)
    nmask = nrow(mask)
    a = summary(mask)$cellarea*10^4
    A = length(mask$x)*a
  },
  run = function(pars = double(1)) {
    # unpack parameters
    D = exp(pars[1])
    l0 = exp(pars[2])
    sig2 = exp(pars[3])
    #D = pars[1]
    #l0 = pars[2]
    #sig2 = pars[3]
    Ens <- numeric(value = 0, length = ntraps)
    # calcuate the encounter function at each mask point (column) for each trap (row)
    for(i in 1:ntraps) {
      #      haz <- 0.0
      for(j in 1:nmask) {
        Ens[i] <- Ens[i] +  l0*exp(-0.5*dists2[i,j]/sig2)
      }
      #  Ens[i]<- D*Ens[i]/nmask*a
      Ens[i]<- D*Ens[i]*a
    }
    return(Ens)
    returnType(double(1))
  },
  methods = list(
    grEns = function(pars=double(1)) {
      #      ans = nimDeriv(run(pars), wrt=1:length(pars), order=0:2)
      #      retrunType(ADNimbleList)
      ans = derivs(run(pars), wrt=1:length(pars), order=0:2)
      return(ans)
      returnType(ADNimbleList())
    }),
  buildDerivs = list(run=list(ignore=c('i', 'j')))
)

# Quasi-likelihood function (used in GEEs) written as an objective function.
# Parameter estimates are parameter values when the (vector) quasi-likelihood 
# function evaluates to (vector) zero, so I make the objective function
# the sum of the squared LH sides of the equation, as this achieves its 
# minimum at zero.
quasifn = function(pars,y,dists2,A,a) {
  D = exp(pars[1])
  Ens = mu.En.hn(pars,dists2,a,by="mask")
  vcv = make.cov(Ens, mu_N=D*A)
  Nimobj = calcEns.c$grEns(pars)
  Nimobj$value
  Nimobj$jacobian
  fn = sum((t(Nimobj$jacobian)%*%solve(vcv)%*%(y-Nimobj$value))^2)
  return(fn)
}



#' @title Plots numbers of detections at each trap.
#'
#' @description
#' Plots a histogram of detection counts by trap or count-sized circles on a plot of traps on 
#' the mask, or both.
#'
#' @param nu Number of detections at each trap. A vector of same length as 
#' object \code{traps}, and in the same order as the traps in this object.
#' @param mask A \code{mask} object.
#' @param traps A \code{traps} object.
#' @param what Type of plot to draw. Must be "both","hist", or "traps".
#' 
#' @export
plot.nu = function(nu,mask=NULL,traps=NULL,what=c("both","hist","traps"),border=0,...) {
  what = match.arg(what)
  if(what=="both") {
    if(is.null(traps)) stop("You have to pass traps.")
    par(mfrow=c(1,2))
    hist(nu,main="",xlab="Number detections",ylab="Trap frequency")
    if(!is.null(mask)) {
      plot(mask,border=border)
      plot(traps,add=TRUE)
    } else plot(traps,border=border)
    points(traps$x, traps$y, cex=nu/max(nu)*3)
    title(paste0("Counts from ",min(nu)," to ",max(nu)),...)
  } else if(what=="hist") {
    hist(nu,main="",xlab="Number detections",ylab="Trap frequency")
  } else {
    if(is.null(traps)) stop("You have to pass traps.")
    if(!is.null(mask)) {
      plot(mask,border=border)
      plot(traps,add=TRUE)
    } else plot(traps,border=border)
    points(traps$x, traps$y, cex=nu/max(nu)*3)
    title(paste0("Counts from ",min(nu)," to ",max(nu)),...)
  }
}
