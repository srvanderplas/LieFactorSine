library(pscl)
library(LearnBayes)
library(multicore)

### Read in Data
turkdata.full <- read.csv("./data/turkdataclean.csv", stringsAsFactors=FALSE)
turkdata <- ddply(subset(turkdata, len>10), .(ip.id, test_param), transform, n=length(test_param))
turkdata <- subset(turkdata, n>4)
turkdata$ip.id.old <- turkdata$ip.id
turkdata$ip.id <- as.numeric(factor(turkdata$ip.id))


### Setup - fixing hyperparameters
# solves for invgamma parameters for priors
invgammaparametersolver <- function(qu,p,init=c(1, 1)) {
  # I think that the p-th quantile occurs at qu
  qu <- qu
  p <- p
  igammaoptim = function(param) { 
    if(param[1]<0) return(Inf)
    if(param[2]<0) return(Inf)
    # print(param)
    q1 <- qu[1]
    q2 <- qu[2]
    p1 <- p[1]
    p2 <- p[2]
    
    (pigamma(q1,param[1],param[2])-p1)^2 + (pigamma(q2,param[1],param[2])-p2)^2
  }
  
  r = optim(init,igammaoptim, control=list(reltol=1e-10, abstol=1e-10))
  v = unlist(r)
  t = c(v[1],v[2])
  print(t)
}

### Prior parameters

gammapars <- invgammaparametersolver(c(.001, .25), c(.1, .8), c(1, 1))
alpha.g <- gammapars[1] # 1.776
beta.g <- gammapars[2] # 0.169
tau2 <- 2
mu0 <- 1.5
eta2 <- 1
sigma.range <- seq(.01, 2, .001)
theta.range <- seq(1, 3, .001)

### Log conditional posteriors
logf.theta <- function(theta, mu, sigma2, x, j){
  sapply(theta, function(k){
    1/tau2*sum((k-mu)^2)-
      sum(1/(2*sigma2[j])*(x-k)^2)-
      sum(log(1-pnorm((1-k)/sigma[j])))
  })
}

logf.mu <- function(mu, theta, sigma2){
  dnorm(mu, mean=((mu0/eta2+sum(theta)/tau2)/(1/eta2+length(theta)/tau2)), sd=1/sqrt(1/eta2+length(theta)/tau2), log=T)
}

logf.sigma2 <- function(sigma2, theta, x, j){
  # thetaj is a list of all thetas, x is a list of all x's, and j is a mapping of which x goes to which thetaj that is the same length as x.
  sapply(sigma2, function(k){
    (-alpha.g-1)*(sum(log(k)))-
      sum(.5*log(2*pi*k)+log(1-pnorm((1-theta[j])/sqrt(k))))-
      beta.g*sum(1/k)-
      sum(1/(2*k^2)*(x-theta[j])^2)
  })
}

### Functions to draw from conditional posteriors
r.mu <- function(theta, sigma2){
  dnorm(1, mean=((mu0/eta2+sum(theta)/tau2)/(1/eta2+length(theta)/tau2)), sd=1/sqrt(1/eta2+length(theta)/tau2))
}

r.theta <- function(mu, sigma2, x, j){
  samples <- unlist(mclapply(unique(j), function(k){
    logs <- logf.theta(theta.range, mu, sigma2[k], x[which(j==k)], rep(1, sum(j==k)))
    logs <- logs-max(logs, na.rm=TRUE)
    logs <- exp(logs)
    logs <- logs/sum(logs)
    
    ints <- which(!is.na(logs))
    theta.range[ints[sample.int(length(ints), size=1, prob=logs[!is.na(logs)])]]
  }, mc.cores=8))
  samples
}

r.sigma2 <- function(theta, x, j){
  # thetaj is a list of all thetas, x is a list of all x's, and j is a mapping of which x goes to which thetaj that is the same length as x.
  samples <- unlist(mclapply(unique(j), function(k){
    logs <- logf.sigma2(sigma.range, theta[k], x[which(j==k)], rep(1, sum(j==k)))
    logs <- logs-max(logs, na.rm=TRUE)
    logs <- exp(logs)
    logs <- logs/sum(logs)
    
    ints <- which(!is.na(logs))
    sigma.range[ints[sample.int(length(ints), size=1, prob=logs[!is.na(logs)])]]
  }, mc.cores=8))
  samples
}

### Gibbs Sampling
N <- 5000
mu.vec <- rep(0, N)
theta.vec <- matrix(0, ncol=length(unique(turkdata$ip.id)), nrow=N)
sigma2.vec <- matrix(0, ncol=length(unique(turkdata$ip.id)), nrow=N)
mu.vec[1] <- 1
theta.vec[1,] <- 1.3
sigma2.vec[1,] <- .5
for(i in 1:(N-1)){
  if(i%%10==0) print(i)
  mu.vec[i+1] <- r.mu(theta.vec[i], sigma2.vec[i,])
  theta.vec[i+1,] <- r.theta(mu.vec[i+1], sigma2.vec[i,], turkdata$ans.liefactor, turkdata$ip.id)
  sigma2.vec[i+1,] <- r.sigma2(theta.vec[i+1,], turkdata$ans.liefactor, turkdata$ip.id)
}