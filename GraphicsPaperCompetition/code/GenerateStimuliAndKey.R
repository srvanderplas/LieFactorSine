# Goal is to get an idea of the size of the correction most people prefer. 
# Uses quadratic correction in Y.

library(ggplot2)
library(reshape2)
library(plyr)
# source("./code/themeStimuli.R")

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

getYlim <- function(w, orig, f, fprime, f2prime){
  temp <- melt(ldply(c(0, 1.4), function(i) weightYTrans(orig, i)[,c(1,5, 6)]), id.vars="x", value.name="y", variable.name="var")
  dy <- diff(apply(temp[,c(1,3)], 2, function(k) diff(range(k))))
  dx <- 0
  if(dy>0) {
    dx <- dy
    dy <- 0
  }
  return(list(dx=range(temp$x)+c(-1, 1)*dx/2, dy=range(temp$y)+c(-1, 1)*-dy/2))
}

#---------------------------------------------------------------------------------------------------------------
w.all <- matrix(c( 0,  .2,  .4,  .8, 1.25, 1.4,    0,  .15, .35, .8,  1.2,  1.4, 
                   0,  .2,  .4,  .6,  .8,  1.0,   .1,  .3,  .5,  .7,   .9,  1.1,
                   0,  .5,  .7,  .8,  .9,  1.0,    0,  .45, .65, .75,  .85, 1.0,
                   .4, .7,  .8,  .9, 1.0,  1.3,   .3,  .65, .75, .85,  .95, 1.2,
                   .5, .6,  .7,  .8,  .9,  1.0,   .55, .65, .75, .85,  .95, 1.0,
                   .6, .7,  .75, .85, .9,  1.0,   .6,  .7,  .75, .8,   .9,  1.0), nrow=12, ncol=6, byrow=TRUE)

w.all <- matrix(c(  0, .2,  .4,  .8,  1.25, 1.4,       0,  .15, .35, .8, 1.2,  1.4, 
                    0, .2,  .4,  .6,   .8,  1.0,      .1,  .3,  .5,  .7,  .9,  1.1,
                  .05, .3,  .5,  .65,  .8,  1.0,      .1,  .3,  .55, .7,  .85, 1.0,
                   .4, .6,  .7,  .8,   .9,  1.05,     .35, .65, .75, .85, .95, 1.05,
                  .35, .5,  .6,  .7,   .8,   .95,     .4,  .55, .65, .75, .85, 1.0,
                   .5, .65, .75, .8,   .9,  1.0,      .5,  .6,  .7,  .75, .85, 1.0), nrow=12, ncol=6, byrow=TRUE)
data.file <- NULL
#---------------------------------------------------------------------------------------------------------------

set.seed(314159)

f <- function(x) 2*sin(x)
fprime <- function(x) 2*cos(x)
f2prime <- function(x) -2*sin(x)
orig <- createSine(50, 1, f, fprime, f2prime, 0, 2*pi)
lims <- getYlim(c(0, 1.4), orig, f, fprime, f2prime)

weight.key <- NULL
for(j in 1:12){
  w <- w.all[j,]
  diff.w <- ceiling(j/2)-1
  frameorder <- sample(w, 6)
  
  data <- rbind.fill(ldply(w, function(i) weightYTrans(orig, i)))
  data$set <- sapply(data$w, function(i) which(w %in% i))
  data$display <- sapply(data$w, function(i) which(frameorder %in% i))
  ans <- which.min(abs(frameorder-.99))
  if(diff.w==0) ans <- which.min(abs(frameorder-.8))
  qplot(data=data, x=xstart, y=ystart, xend=xend, yend=yend, geom="segment") + facet_wrap(~display) + coord_equal(ratio=1) + theme_stimuli() + xlim(lims$dx) + ylim(lims$dy)
#   ggsave(filename=paste("YtransSin", j, ".png", sep=""), path="./stimuli", width=5.37, height=3.88, units="in", dpi=100)
  data.file <- rbind(data.file, data.frame(
    pic_id=j, 
    sample_size=25, 
    test_param=0, 
    param_value=paste("w", median(w), diff.w, sep="-"), 
    p_value=0,
    obs_plot_location=ans,
    pic_name=paste("YtransSin", j, ".png", sep=""),
    experiment="turk8",
    difficulty=diff.w))
  temp <- unique(data[,14:15] )
  weight.key <- rbind(weight.key, 
                      data.frame(pic_name=paste("YtransSin", j, ".png", sep=""), 
                                 w=w[temp$set], set=temp$set, response_no=temp$display,
                                 lmax=sapply(temp$display, function(i) max(subset(data ,display==i)$elltrans))))
}

w.all <- w.all[3:12,]
f <- function(x) exp(x/2)
fprime <- function(x) 1/2*exp(x/2)
f2prime <- function(x) 1/4*exp(x/2)
orig <- createSine(50, 1, f, fprime, f2prime, -pi, pi)
lims <- getYlim(c(0, 1.2), orig, f, fprime, f2prime)
  
for(j in 1:10){
  w <- w.all[j,]
  diff.w <- ceiling(j/2)
  frameorder <- sample(w, 6)
  
  data <- rbind.fill(ldply(w, function(i) weightYTrans(orig, i)))
  data$set <- sapply(data$w, function(i) which(w %in% i))
  data$display <- sapply(data$w, function(i) which(frameorder %in% i))
  ans <- which.min(abs(frameorder-.99))
  if(diff.w==0) ans <- which.min(abs(frameorder-.5))
  qplot(data=data, x=xstart, y=ystart, xend=xend, yend=yend, geom="segment") + facet_wrap(~display) + coord_equal(ratio=1) + theme_stimuli() + xlim(lims$dx) + ylim(lims$dy)
#   ggsave(filename=paste("YTransExp", j, ".png", sep=""), path="./stimuli", width=5.37, height=3.88, units="in", dpi=100)
  data.file <- rbind(data.file, data.frame(
    pic_id=j, 
    sample_size=13, 
    test_param=1, 
    param_value=paste("w", median(w), diff.w, sep="-"), 
    p_value=0,
    obs_plot_location=ans,
    pic_name=paste("YTransExp", j, ".png", sep=""),
    experiment="turk8",
    difficulty=diff.w))
  temp <- unique(data[,14:15] )
  weight.key <- rbind(weight.key, 
                      data.frame(pic_name=paste("YTransExp", j, ".png", sep=""), 
                                 w=w[temp$set], set=temp$set, response_no=temp$display,
                                 lmax=sapply(temp$display, function(i) max(subset(data ,display==i)$elltrans))))
}


f <- function(x) 5/6*1/x
fprime <- function(x) -5/6*x^(-2)
f2prime <- function(x) 2*5/6*x^(-3)
orig <- createSine(50, 1, f, fprime, f2prime, 1/2, 3.5)
lims <- getYlim(c(0, 1.4), orig, f, fprime, f2prime)

for(j in 1:10){
  w <- w.all[j,]
  diff.w <- ceiling(j/2)
  frameorder <- sample(w, 6)
  data <- rbind.fill(ldply(w, function(i) weightYTrans(orig, i)))
  data$set <- sapply(data$w, function(i) which(w %in% i))
  data$display <- sapply(data$w, function(i) which(frameorder %in% i))
  ans <- which.min(abs(frameorder-.99))
  if(diff.w==0) ans <- which.min(abs(frameorder-.5))
  qplot(data=data, x=xstart, y=ystart, xend=xend, yend=yend, geom="segment") + facet_wrap(~display) + coord_equal(ratio=1) + theme_stimuli()  + xlim(lims$dx) + ylim(lims$dy)
#   ggsave(filename=paste("YTransInv", j, ".png", sep=""), path="./stimuli", width=5.37, height=3.88, units="in", dpi=100)
  data.file <- rbind(data.file, data.frame(
    pic_id=j, 
    sample_size=12, 
    test_param=2, 
    param_value=paste("w", median(w), diff.w, sep="-"), 
    p_value=0,
    obs_plot_location=ans,
    pic_name=paste("YTransInv", j, ".png", sep=""),
    experiment="turk8",
    difficulty=diff.w))
  temp <- unique(data[,14:15] )
  weight.key <- rbind(weight.key, 
                      data.frame(pic_name=paste("YTransInv", j, ".png", sep=""), 
                                 w=w[temp$set], set=temp$set, response_no=temp$display,
                                 lmax=sapply(temp$display, function(i) max(subset(data ,display==i)$elltrans))))
}

write.csv(data.file, "./stimuli/datafile.csv", row.names=FALSE)
write.csv(weight.key, "./data/pictureKey.csv", row.names=FALSE)
write.csv(weight.key, "/home/susan/Dropbox/GraphicsGroup/TurkLieFactorSine/pictureKey.csv", row.names=FALSE)

# examples
w <- c(0, .9, 1.4)
frameorder <- c(.9, 0, 1.4)
f <- function(x) 2*sin(x)
fprime <- function(x) 2*cos(x)
f2prime <- function(x) -2*sin(x)
orig <- createSine(50, 1, f, fprime, f2prime, 0, 2*pi)
lims <- getYlim(c(0, 1.4), orig, f, fprime, f2prime)

data <- rbind.fill(ldply(w, function(i) weightYTrans(orig, i)))
data$set <- sapply(data$w, function(i) which(w %in% i))
data$display <- sapply(data$w, function(i) which(frameorder %in% i))

qplot(data=data, x=xstart, y=ystart, xend=xend, yend=yend, geom="segment") + facet_wrap(~display) + theme_stimuli() + xlim(lims$dx) + ylim(lims$dy)
# ggsave("./stimuli/Example2.png",width=5.37, height=2.10, units="in", dpi=100)

x <- orig$x
data <- rbind(data.frame(x=x, y=x, ystart=x-.5, yend=x+.5, display=1),
              data.frame(x=x, y=x, ystart=.5*abs(x-pi/2)+.2, yend=1.5*abs(x-pi)-.2, display=2), 
              data.frame(x=x, y=x, ystart=x-sin(x), yend=x+sin(x), display=3))
qplot(data=data, x=x, y=ystart, xend=x, yend=yend, geom="segment") + facet_wrap(~display) + theme_stimuli() + xlim(c(-3/8+min(orig$x), 3/8+max(orig$x)))
# ggsave("./stimuli/Example1.png",width=5.37, height=2.10, units="in", dpi=100)

x <- orig$x
data <- rbind(data.frame(x=x, xend=x, xstart=x, y=x, ystart=1.5*abs(x-pi)-.25, yend=1.5*abs(x-pi)+.75, display=3),
              data.frame(x=x, xend=x, xstart=x, y=x, ystart=.75*x-.25+cos(x)*.5, yend=.75*x+.25+cos(x)*.5 + x*.1, display=2),
              data.frame(x=x, xend=x, xstart=x, y=x, ystart=.75*x-.25+.5*log(.25*x+1), yend=.75*x+.25+.5*log(4*x+1), display=1))
qplot(data=data, x=x, y=ystart, xend=x, yend=yend, geom="segment") + facet_wrap(~display) + theme_stimuli() + ylim(range(c(data$ystart, data$yend)) + c(-.183, .183))
# ggsave("./stimuli/Example3.png",width=5.37, height=2.10, units="in", dpi=100)
