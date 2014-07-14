library(pscl)
library(LearnBayes)
library(multicore)
library(compiler)
library(plyr)
library(reshape2)
library(ggplot2)
library(msm)

### Read in Data
turkdata.full <- read.csv("./data/turkdataclean.csv", stringsAsFactors=FALSE)
turkdata <- ddply(subset(turkdata.full, len>10), .(ip.id, test_param), transform, n=length(test_param))
turkdata <- subset(turkdata, n>4)
turkdata <- subset(turkdata, test_param=="Sin")
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

gammapars <- invgammaparametersolver(c(.001, 3), c(.1, .8), c(2, 2))
alpha.g <- gammapars[1] # 3.438
beta.g <- gammapars[2] # 0.466
tau2 <- 2
mu0 <- 1.5
eta2 <- 1

### Log conditional posteriors
logf.thetaj <- cmpfun(function(theta, mu, sigma2, x){
  # for a single theta, sigma2, and corresponding subset of x values
  nj <- length(x)
  -1/(2*tau2)*(theta-mu)^2-
    sum(1/(2*sigma2)*(x-theta)^2)-
    nj*pnorm((1-theta)/sigma2, lower.tail=FALSE, log=TRUE)-
    c(0,Inf)[(sigma2<0)+1] - 2*pi*tau2 -c(0,Inf)[(theta<1)+1]
})

logf.theta.vec <- cmpfun(function(theta, mu, sigma2, x, j){
  ddply(data.frame(theta=theta[j], mu=mu, sigma2 = sigma2[j], x=x, j=j), .(j), summarise, f=logf.thetaj(theta[1], mu[1], sigma2[1], x))$f
})

logf.mu <- cmpfun(function(mu, theta, sigma2){
  dnorm(mu, mean=((mu0/eta2+sum(theta)/tau2)/(1/eta2+length(theta)/tau2)), sd=1/sqrt(1/eta2+length(theta)/tau2), log=T)
})

logf.sigma2j <- cmpfun(function(sigma2, theta, x, j){
  # for a single theta, sigma2, and corresponding subset of x values
  nj <- length(x)
  (-alpha.g-1-nj/2)*log(sigma2)-
    nj*pnorm((1-theta)/sigma2, lower.tail=FALSE, log=TRUE)-
    1/sigma2*(beta.g+.5*sum((x-theta)^2))-c(0,Inf)[(sigma2<0)+1] + 
    alpha.g*log(beta.g)-log(gamma(alpha.g)) - nj/2*log(2*pi)
})

logf.sigma.vec <- cmpfun(function(theta, sigma2, x, j){
  ddply(data.frame(theta=theta[j], sigma2 = sigma2[j], x=x, j=j), .(j), summarise, f=logf.sigma2j(theta[1], sigma2[1], x))$f
})

logf <- cmpfun(function(theta, mu, sigma2, x, j){
  nj <- tapply(j, j, length)
  xj <- tapply(x, j, c)
  sum(log(densigamma(sigma2, alpha.g, beta.g)))+sum(dnorm(theta, mu, sqrt(tau2), log=TRUE)) + dnorm(mu, mu0, sqrt(eta2), log=TRUE) + sum(sapply(1:length(x), function(k) dtnorm(x, theta[j[k]], sqrt(sigma2[j[k]]), lower=1, log=TRUE)))
})

### Functions to draw from conditional posteriors

var.par.2 <- .1
r.mu <- cmpfun(function(theta, mu.old, sigma2, x, j){
  mu.old.mean <- as.numeric((mu.old/eta2 + sum(theta)/sigma2)/(1/eta2+sum(theta)/tau2))
  mu.sd <- as.numeric(1/(1/eta2+sum(theta)/tau2))
  if(sum(mu.sd<=0)>0){
    mu.sd <- pmax(0.0001, th.sd)
    message("sd values negative")
  } 
  mu.sd <- sqrt(mu.sd)/var.par.2
  mu.new <- rnorm(1, mean=mu.old.mean, sd=mu.sd)
  
#   mu.new.mean <- as.numeric((theta.new/tau2 + tapply(x, j, sum)/sigma2)/(1/tau2+tapply(x, j, length)/sigma2))
#   
  logf.mu.old <- logf(theta, mu.old, sigma2, x, j)
  logf.mu.new <- logf(theta, mu.new, sigma2, x, j)
  
  
  r <- (logf.mu.new)/(logf.mu.old)
  u <- runif(length(r))
  return(data.frame(accept=(r>=u), var=as.numeric(mu.new*(r>=u)+mu.old+(r<u))))
})

var.par <- .25

r.theta <- cmpfun(function(theta.old, mu, sigma2, x, j){
  th.old.mean <- as.numeric((theta.old/tau2 + tapply(x, j, sum)/sigma2)/(1/tau2+tapply(x, j, length)/sigma2))
  th.sd <- as.numeric(1/(1/tau2+tapply(x, j, length)/sigma2))
  if(sum(th.sd<=0)>0){
    th.sd <- pmax(0.0001, th.sd)
    message("sd values negative")
  } 
  th.sd <- sqrt(th.sd)/var.par
  theta.new <- rnorm(length(theta.old), mean=th.old.mean, sd=th.sd)
  
  th.new.mean <- as.numeric((theta.new/tau2 + tapply(x, j, sum)/sigma2)/(1/tau2+tapply(x, j, length)/sigma2))
  
  logf.theta.old <- logf(theta.old, mu, sigma2, x, j)
  logf.theta.new <- logf(theta.new, mu, sigma2, x, j)


  r <- (logf.theta.new)/(logf.theta.old)
  u <- runif(length(r))
  return(data.frame(accept=(r>=u), var=as.numeric(theta.new*(r>=u)+theta.old+(r<u))))
})

r.sigma2 <- cmpfun(function(sigma.old, mu, theta, x, j){
  nj <- as.numeric(tapply(x, j, length))
  old.sum <- ddply(data.frame(x=x, theta=theta[j], j=j), .(j), summarise, z=sum((x-theta)^2))$z
  sig.alpha <- alpha.g+nj/2
  if(sum(sig.alpha<=2)>0) message("alpha values < 2")
  sig.alpha <- pmax(sig.alpha, 2)
  sig.beta <- beta.g+.5*old.sum
  sigma.new <- rigamma(length(sigma.old), sig.alpha, sig.beta)
  
  logf.sigma.old <- logf(theta, mu, sigma.old, x, j)
  logf.sigma.new <- logf(theta, mu, sigma.new, x, j)
  jnew.old <- sapply(1:length(sigma.new), function(i) densigamma(sigma.new[i], sig.alpha[i], sig.beta[i]))
  jold.new <- sapply(1:length(sigma.old), function(i) densigamma(sigma.old[i], sig.alpha[i], sig.beta[i]))
  
  r <- (logf.sigma.new/jnew.old)/(logf.sigma.old/jold.new)
#   r <- logf.sigma.new/logf.sigma.old
  u <- runif(length(r))
  return(data.frame(accept=(r>=u), var=as.numeric(sigma.new*(r>=u)+sigma.old+(r<u))))
})
  

### Check posterior conditionals
# vars <- expand.grid(theta=seq(.5, 3, .05), mu=seq(0, 2, .05), sigma=seq(.01, 1, .01))
# # vars$f.theta <- sapply(1:nrow(vars), function(k) logf.thetaj(vars$theta[k], vars$mu[k], vars$sigma[k], turkdata$ans.liefactor[which(turkdata$ip.id==1)]))
# # vars$f.mu <- sapply(1:nrow(vars), function(k) logf.mu(vars$mu[k], vars$theta[k], vars$sigma[k]))
# # vars$f.sigma <- sapply(1:nrow(vars), function(k) logf.sigma2j(vars$sigma[k], vars$theta[k], turkdata$ans.liefactor[which(turkdata$ip.id==1)]))
# 
# # qplot(data=vars, x=theta, y=f.theta, group=interaction(vars$mu, vars$sigma), geom="line", alpha=I(.2), colour=sigma)
# 
# vars$f <- unlist(mclapply(1:nrow(vars), function(k) logf(vars$theta[k], vars$mu[k], vars$sigma[k], turkdata$ans.liefactor[which(turkdata$ip.id==1)], rep(1, sum(turkdata$ip.id==1)))))
# vars$f[!is.finite(vars$f)] <- NA
# vars$f <- vars$f-max(vars$f[is.finite(vars$f)])
# vars$f <- exp(vars$f)
# vars$f <- vars$f/sum(vars$f)
# 
# 
# 
# qplot(data=subset(vars, mu==1.5), x=sigma, y=f, group=theta, geom="line", alpha=I(.2), colour=theta)
# qplot(data=vars, x=theta, y=sigma, z=f.theta, group=mu, geom="contour")+xlim(c(1, 10))
# qplot(data=subset(vars, mu==1.5), x=theta, y=sigma, z=f, geom="contour")
# qplot(data=subset(vars, round(sigma,2)==.05), x=theta, y=mu, fill=f, geom="tile")
# 
# qplot(data=vars[which(is.finite(vars$f.sigma)),], x=theta, y=sigma, geom="jitter", alpha=I(.01))
# 
# qplot(data=vars[which(is.finite(vars$f.sigma)),], x=theta, y=mu, geom="jitter", alpha=I(.01))
# qplot(data=vars[which(is.finite(vars$f.sigma)),], x=mu, y=sigma, geom="jitter", alpha=I(.01))

### Gibbs Sampling
N <- 100
mu.vec <- rep(0, N)
mu.acc <- 0
theta.vec <- matrix(0, ncol=length(unique(turkdata$ip.id)), nrow=N)
theta.acc <- rep(0, length(unique(turkdata$ip.id)))
sigma2.vec <- matrix(0, ncol=length(unique(turkdata$ip.id)), nrow=N)
sigma.acc <- rep(0, length(unique(turkdata$ip.id)))
mu.vec[1] <- 1.388 # initialize to avg turk lie factor
theta.vec[1,] <- unlist(dlply(turkdata, .(ip.id), function(i) mean(i$ans.liefactor)))
sigma2.vec[1,] <- jitter(as.numeric(tapply(turkdata$ans.liefactor, turkdata$ip.id, sd)))
for(i in 1:(N-1)){
  if(i%%10==0) print(i)
  temp <- r.mu(theta = theta.vec[i,],
                      mu = mu.vec[i], 
                      sigma2 = sigma2.vec[i,], 
                      x = turkdata$ans.liefactor, 
                      j = turkdata$ip.id)
  mu.vec[i+1] <- temp$var
  mu.acc <- mu.acc+temp$accept
  temp <- r.theta(theta.old = theta.vec[i,],
                   mu = mu.vec[i+1], 
                   sigma2 = sigma2.vec[i,], 
                   x = turkdata$ans.liefactor, 
                   j = turkdata$ip.id)
  theta.vec[i+1,] <- temp$var
  theta.acc <- theta.acc+temp$accept
  temp <- r.sigma2(sigma.old = sigma2.vec[i,], 
                   mu = mu.vec[i+1],
                   theta = theta.vec[i+1,], 
                   x = turkdata$ans.liefactor, 
                   j = turkdata$ip.id)
  sigma2.vec[i+1,] <- temp$var
  sigma.acc <- sigma.acc + temp$accept
}

theta.df <- melt(theta.vec, varnames=c("i", "ip.id"), value.name="theta.i")
sigma.df <- melt(sigma2.vec, varnames=c("i", "ip.id"), value.name="sigma.i")
pars <- merge(theta.df, sigma.df)
pars$mu <- mu.vec[pars$i]
pars <- pars[order(pars$i, pars$ip.id),]

qplot(data=pars, x=i, y=theta.i, geom="line", group=ip.id, alpha=I(.25))+ylim(c(0,2.5))
qplot(data=pars, x=i, y=sigma.i, geom="line", group=ip.id, alpha=I(.25))
qplot(data=subset(pars, ip.id==1), x=i, y=mu, geom="line")
