# Goal is to get an idea of the size of the correction most people prefer. Uses quadratic correction in Y.

library(ggplot2)
library(reshape2)
library(plyr)
source("./code/themeStimuli.R")
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
orig <- createSine(50, 1, f, fprime, f2prime, 0, 2*pi)

w.all <- rbind(seq(0, 1.3, len=6), 
           c(0, seq(.2, .6, len=4), 1.3), 
           c(0, seq(.5, .9, len=4), 1.3), 
           c(0, seq(.8, 1.2, len=4), 1.3),
           c(0, .4, .5, .6, .7, 1.3), c(0, .5, .6, .7, .8, 1.3), c(0, .6, .7, .8, .9, 1.3), 
           c(0, .7, .8, .9, 1, 1.3), c(0, .8, .9, 1, 1.1, 1.3), c(.4, .7, .75, .8, .85, 1))

data.file <- NULL
for(j in 1:10){
  w <- w.all[j,]
  diff.w <- min(ceiling((diff(w)/.05)))
  frameorder <- sample(w, 6)
  
  data <- rbind.fill(ldply(w, function(i) weightYTrans(orig, i)))
  data$set <- sapply(data$w, function(i) which(w %in% i))
  data$display <- sapply(data$w, function(i) which(frameorder %in% i))
  ans <- which.min(abs(frameorder-.99))
  qplot(data=data, x=xstart, y=ystart, xend=xend, yend=yend, geom="segment") + facet_wrap(~display) + coord_equal(ratio=1) + theme_stimuli()
  ggsave(filename=paste("YTransPic", j, ".png", sep=""), path="./stimuli", width=8, height=6, units="in")
  data.file <- rbind(data.file, data.frame(
    pic_id=j, 
    sample_size=50, 
    test_param="weight of y-transformation", 
    param_value=paste("w", median(w), diff.w, sep="-"), 
    p_value=0,
    obv_plot_location=ans,
    pic_name=paste("YTransPic", j, ".png", sep=""),
    experiment="turk8",
    difficulty=diff.w))
}
write.csv(data.file, "./stimuli/datafile.csv", row.names=FALSE)

# examples
w <- c(0, 1, 1.2)
frameorder <- c(1, 0, 1.2)
data <- rbind.fill(ldply(w, function(i) weightYTrans(orig, i)))
data$set <- sapply(data$w, function(i) which(w %in% i))
data$display <- sapply(data$w, function(i) which(frameorder %in% i))

qplot(data=data, x=xstart, y=ystart, xend=xend, yend=yend, geom="segment") + facet_wrap(~display) + theme_stimuli()
ggsave("./stimuli/Examples/Example2.png",width=6, height=3)

x <- seq(0, 2*pi, .1)
data <- rbind(data.frame(x=x, y=x, ystart=x-.5, yend=x+.5, display=1),
              data.frame(x=x, y=x, ystart=.5*abs(x-pi/2)+.2, yend=1.5*abs(x-pi)-.2, display=2), 
              data.frame(x=x, y=x, ystart=x-sin(x), yend=x+sin(x), display=3))
qplot(data=data, x=x, y=ystart, xend=x, yend=yend, geom="segment") + facet_wrap(~display) + theme_stimuli()
ggsave("./stimuli/Examples/Example1.png",width=6, height=3)

x <- seq(0, 2*pi, .1)
data <- rbind(data.frame(x=x, xend=x, xstart=x, y=x, ystart=2*abs(x-pi)-.25, yend=2*abs(x-pi)+.25, display=3),
              data.frame(x=x, xend=x, xstart=x, y=x, ystart=x-.25+cos(x)*.25, yend=x+.25+cos(x)*.25 + x*.075, display=2),
              data.frame(x=x, xend=x, xstart=x, y=x, ystart=x-.25+.5*log(x+1), yend=x+.25+.5*log(2*x+1), display=1))
qplot(data=data, x=x, y=ystart, xend=x, yend=yend, geom="segment") + facet_wrap(~display) + theme_stimuli()
ggsave("./stimuli/Examples/Example3.png",width=6, height=3)

