# Code to generate 2x2 grid of pictures, 3 or 4 of which have similar or identical slopes.
library(ggplot2)
library(plyr)
library(reshape2)
source("./code/themeStimuli.R")

# Input number of line segments and a vector of nonnegative slope magnitudes. 
# Output a dataset of x, y, slope, and segment number
getDataSet <- function(nsegment=4, slope=c(1))
{
  segmentLength=2:(nsegment+1)
  initialSlope <- rbind(rep(c(-1, 1), times=ceiling(nsegment/2))[1:nsegment], 
                        rep(c(1, -1), times=ceiling(nsegment/2))[1:nsegment])
  x <- seq(1, sum(segmentLength), length=10*sum(segmentLength))
  slope <- rep(slope, times=10*ceiling(nsegment/length(slope)))
  genY <- data.frame(slope = sample(slope, nsegment, replace=FALSE)*initialSlope[sample(c(1, 2), 1), ], 
                     times = sample(segmentLength, nsegment, replace=FALSE))
  slope <- unlist(sapply(1:nsegment, function(i) rep(genY$slope[i], times=10*genY$times[i])))
  segment <- unlist(sapply(1:nsegment, function(i) rep(i, times=10*genY$times[i])))
  data <- data.frame(x=x, y=cumsum(slope)/10, slope=slope, segmentLength=segment)
  data$ynorm <- (data$y-mean(data$y))/sd(data$y)
  data
}

get4 <- function(order){
  ans <- which(order==1)
#   a <- .75# *(order==1) + .95*(order!=1)
#   b <- 2# *(order==1) + 1.05*(order!=1)
  nslopes <- runif(6, min=1, max=5)
  nseg <- 8
  data <- do.call("rbind", lapply(1:length(order), function(i) data.frame(
    getDataSet(nseg, slope=nslopes),
    set=i, correct=ans)))
  data$sd <- .125*with(data, (set==ans) + (set!=ans)*sqrt(1+slope^2))
#   data <- ddply(data, .(set), transform, sd = sd-mean(sd)+.25)
  data$sd[which(data$set==ans)] <- mean(data$sd[which(data$set!=ans)])
  data <- ddply(data, .(set), transform, sd2 = loess(sd~x, span=.2, degree=1, surface="interpolate")$fitted)
#   data <- ddply(data, .(set), transform, sd2 = sd2-mean(sd2)+.125)
}

data <- get4(sample(1:4, 4))
data.points <- ddply(data, .(set, x, segmentLength, correct, slope), summarise, ynorm=rnorm(10, ynorm, sd2))

# Ribbon plot
qplot(data=data, x=x, ymin=ynorm-sd2, ymax=ynorm+sd2, geom="ribbon") + 
  facet_wrap(~set) + theme_stimuli() 

qplot(data=data, x=x, y=sd2, geom="line") + facet_wrap(~set) #+ theme_stimuli() 
# Points
qplot(data=data.points, x=x, y=ynorm, geom="point", colour=I("grey10")) + 
  geom_line(data=data, aes(x=x, y=ynorm), color="black") +
  facet_wrap(~set) + theme_stimuli()

data$correct[1]