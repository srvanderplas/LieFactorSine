library(reshape2)
library(plyr)
library(xtable)

createSine <- function(n=200, len=1, f=f, fprime=fprime, f2prime=f2prime, a=0, b=2*pi) {
  x <- seq(a, b, length=n+2)[(2:(n+1))]
  ell <- rep(len, length=length(x))
  fx <- f(x)
  ystart <- fx - .5*ell
  yend <- fx + .5*ell
  
  # now correct for line illusion in vertical direction
  dy <- diff(range(fx))
  dx <- diff(range(x))
  # fprime works in framework of dx and dy, but we represent it in framework of dx and dy+len
  # needs to be fixed by factor a:  
  a <- dy/(dy + len) 
  # ellx is based on the "trig" correction
  ellx <- ell / cos(atan(abs(a*fprime(x))))
  # ellx2 is based on linear approximation of f  
  ellx2 <- ell * sqrt(1 + a^2*fprime(x)^2)
  
  # make this a data frame - ggplot2 doesn't do well with floating vectors
  dframe <- data.frame(x=x, xstart=x, xend=x, y=fx, ystart=ystart, yend=yend, ell=ell, ellx = ellx, ellx2=ellx2)
  
  # third adjustment is based on quadratic approximation of f.
  # this needs two parts: correction above and below f(x)  
  
  fp <- a*fprime(x)
  f2p <- a*f2prime(x)
  lambdap <- (sqrt((fp^2+1)^2-f2p*fp^2*ell) + fp^2 + 1)^-1    
  lambdam <- -(sqrt((fp^2+1)^2+f2p*fp^2*ell) + fp^2 + 1)^-1    
  
  
  dframe$ellx4.l <- (4*abs(lambdap)*sqrt(1+fp^2))^-1
  dframe$ellx4.u <- (4*abs(lambdam)*sqrt(1+fp^2))^-1
  
  dframe
}

weightYTrans <- function(df, w){
  df$elltrans <- w*df$ellx2/2 + (1-w)*df$ell/2
  df$ystart <- df$y - df$elltrans
  df$yend <- df$y + df$elltrans
  df$w <- w
  df
}

f <- function(x) 2*sin(x)
fprime <- function(x) 2*cos(x)
f2prime <- function(x) -2*sin(x)
orig <- createSine(40, 1, f, fprime, f2prime, 0, 2*pi)

w.all <- matrix(c(  0, .2,  .4,  .8,  1.25, 1.4,       0,  .15, .35, .8, 1.2,  1.4, 
                    0, .2,  .4,  .6,   .8,  1.0,      .1,  .3,  .5,  .7,  .9,  1.1,
                  .05, .3,  .5,  .65,  .8,  1.0,      .1,  .3,  .55, .7,  .85, 1.0,
                   .4, .6,  .7,  .8,   .9,  1.05,     .35, .65, .75, .85, .95, 1.05,
                  .35, .5,  .6,  .7,   .8,   .95,     .4,  .55, .65, .75, .85, 1.0,
                   .5, .65, .75, .8,   .9,  1.0,      .5,  .6,  .7,  .75, .85, 1.0), nrow=12, ncol=6, byrow=TRUE)


w.long <- melt(w.all)
names(w.long) <- c("Stimuli", "Plot", "w")


w.long$d <- sapply(w.long$w, function(i) {
                                temp <- weightYTrans(df=createSine(40, 1, f, fprime, f2prime, 0, 2*pi), w=i)$elltrans
                                abs(diff(range(temp)))/min(temp)
                              }
                   )

w.long$d2 <- sapply(w.long$w, function(i) {
                                  temp <- weightYTrans(df=createSine(40, 1, f, fprime, f2prime, 0, 2*pi), w=i)$elltrans
                                  max(temp)/min(temp)
                                }
)
sine.long <- w.long
sine.long$type="Sine"

w.wide.w <- dcast(w.long, Stimuli ~ Plot, value.var="w")
w.wide.w$difficulty <- c("test", "test", "1", "1", "2", "2", "3", "3", "4", "4", "5", "5")
w.wide.w$blank <- ""
w.wide.w <- w.wide.w[,c(1, 8, 2:7, 9)]
w.wide.d <- dcast(w.long, Stimuli ~ Plot, value.var="d2")


sine.table <- cbind(w.wide.w, round(w.wide.d[,-1], 2))
print(xtable(sine.table), include.rownames=FALSE)

f <- function(x) exp(x/2)
fprime <- function(x) 1/2*exp(x/2)
f2prime <- function(x) 1/4*exp(x/2)
orig <- createSine(40, 1, f, fprime, f2prime, -pi, pi)

w.all <- matrix(c(  0, .2,  .4,  .8,  1.25, 1.4,       0,  .15, .35, .8, 1.2,  1.4, 
                    0, .2,  .4,  .6,   .8,  1.0,      .1,  .3,  .5,  .7,  .9,  1.1,
                    .05, .3,  .5,  .65,  .8,  1.0,      .1,  .3,  .55, .7,  .85, 1.0,
                    .4, .6,  .7,  .8,   .9,  1.05,     .35, .65, .75, .85, .95, 1.05,
                    .35, .5,  .6,  .7,   .8,   .95,     .4,  .55, .65, .75, .85, 1.0,
                    .5, .65, .75, .8,   .9,  1.0,      .5,  .6,  .7,  .75, .85, 1.0), nrow=12, ncol=6, byrow=TRUE)


w.long <- melt(w.all)
names(w.long) <- c("Stimuli", "Plot", "w")


w.long$d <- sapply(w.long$w, function(i) {
  temp <- weightYTrans(df=createSine(40, 1, f, fprime, f2prime, -pi, pi), w=i)$elltrans
  abs(diff(range(temp)))/min(temp)
}
)

w.long$d2 <- sapply(w.long$w, function(i) {
  temp <- weightYTrans(df=createSine(40, 1, f, fprime, f2prime, -pi, pi), w=i)$elltrans
  max(temp)/min(temp)
}
)

exp.long <- w.long
exp.long$type="Exponential"
w.long <- ddply(w.long, .(Stimuli), transform, minlie=min(d2))
w.long$liefactor <- w.long$d2/w.long$minlie

w.wide.w <- dcast(w.long, Stimuli ~ Plot, value.var="w")
w.wide.w$difficulty <- c("test", "test", "1", "1", "2", "2", "3", "3", "4", "4", "5", "5")
w.wide.w$blank <- ""
w.wide.w <- w.wide.w[,c(1, 8, 2:7, 9)]
w.wide.d <- dcast(w.long, Stimuli ~ Plot, value.var="d2")


exp.table <- cbind(w.wide.w, round(w.wide.d[,-1], 2))
print(xtable(exp.table), include.rownames=FALSE)


f <- function(x) 5/6*1/x
fprime <- function(x) -5/6*x^(-2)
f2prime <- function(x) 2*5/6*x^(-3)
orig <- createSine(40, 1, f, fprime, f2prime, 1/2, 3.5)

w.all <- matrix(c(  0, .2,  .4,  .8,  1.25, 1.4,       0,  .15, .35, .8, 1.2,  1.4, 
                    0, .2,  .4,  .6,   .8,  1.0,      .1,  .3,  .5,  .7,  .9,  1.1,
                    .05, .3,  .5,  .65,  .8,  1.0,      .1,  .3,  .55, .7,  .85, 1.0,
                    .4, .6,  .7,  .8,   .9,  1.05,     .35, .65, .75, .85, .95, 1.05,
                    .35, .5,  .6,  .7,   .8,   .95,     .4,  .55, .65, .75, .85, 1.0,
                    .5, .65, .75, .8,   .9,  1.0,      .5,  .6,  .7,  .75, .85, 1.0), nrow=12, ncol=6, byrow=TRUE)


w.long <- melt(w.all)
names(w.long) <- c("Stimuli", "Plot", "w")


w.long$d <- sapply(w.long$w, function(i) {
  temp <- weightYTrans(df=createSine(40, 1, f, fprime, f2prime, 1/2, 3.5), w=i)$elltrans
  abs(diff(range(temp)))/min(temp)
}
)

w.long$d2 <- sapply(w.long$w, function(i) {
  temp <- weightYTrans(df=createSine(40, 1, f, fprime, f2prime, 1/2, 3.5), w=i)$elltrans
  max(temp)/min(temp)
}
)
inv.long <- w.long
inv.long$type <- "Inverse"

w.long <- ddply(w.long, .(Stimuli), transform, minlie=min(d2))
w.long$liefactor <- w.long$d2/w.long$minlie

w.wide.w <- dcast(w.long, Stimuli ~ Plot, value.var="w")
w.wide.w$difficulty <- c("test", "test", "1", "1", "2", "2", "3", "3", "4", "4", "5", "5")
w.wide.w$blank <- ""
w.wide.w <- w.wide.w[,c(1, 8, 2:7, 9)]
w.wide.d <- dcast(w.long, Stimuli ~ Plot, value.var="d2")


inv.table <- cbind(w.wide.w, round(w.wide.d[,-1], 2))
print(xtable(inv.table), include.rownames=FALSE)

library(ggplot2)
w.long <- rbind(sine.long, inv.long, exp.long)

qplot(data=w.long, x=w, y=d, color=type) + facet_wrap(~Stimuli)

qplot(data=w.long, y=Stimuli, x=d) + facet_wrap(~type)
