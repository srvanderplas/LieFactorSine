# Code to generate 2x2 grid of pictures, 3 or 4 of which have similar or identical slopes.
library(ggplot2)
library(plyr)
library(reshape2)
source("./code/themeStimuli.R")

# Input number of line segments and a vector of nonnegative slope magnitudes. 
# Output a dataset of x, y, slope, and segment number
getDataSet <- function(nsegment=4, slope=c(1)){
  segmentLength=2:(nsegment+1)
  initialSlope <- rbind(rep(c(-1, 1), times=ceiling(nsegment/2))[1:nsegment], 
                        rep(c(1, -1), times=ceiling(nsegment/2))[1:nsegment])
  x <- seq(1, sum(segmentLength), length=2*sum(segmentLength))
  slope <- rep(slope, times=2*ceiling(nsegment/length(slope)))
  genY <- data.frame(slope = sample(slope, nsegment, replace=FALSE)*initialSlope[sample(c(1, 2), 1), ], 
                     times = sample(segmentLength, nsegment, replace=FALSE))
  slope <- unlist(sapply(1:nsegment, function(i) rep(genY$slope[i], times=2*genY$times[i])))
  segment <- unlist(sapply(1:nsegment, function(i) rep(i, times=2*genY$times[i])))
  data <- data.frame(x=x, y=cumsum(slope)/2, slope=slope, segmentLength=segment)
  data$ynorm <- (data$y-mean(data$y))/sd(data$y)
  data
}

get4 <- function(order){
  ans <- which(order==1)
#   a <- .75# *(order==1) + .95*(order!=1)
#   b <- 2# *(order==1) + 1.05*(order!=1)
  nslopes <- sample(c(runif(2*length(order), .5, 1), runif(length(order), 1, 1.5)), 1.5*length(order))
  nseg <- 5
  data <- do.call("rbind", lapply(1:length(order), function(i) data.frame(
    getDataSet(nseg, slope=nslopes),
    set=i, correct=ans)))
  data$sd <- with(data, (set==ans) + (set!=ans)*sqrt(1+slope^2))
#   data <- ddply(data, .(set), transform, sd = sd-mean(sd)+.25)
  data$sd[which(data$set==ans)] <- mean(ddply(data[which(data$set!=ans),], .(set), summarise, meansd=mean(sd))$meansd)
  data <- ddply(data, .(set), transform, sd2 = loess(sd~x, span=.1, degree=1, surface="interpolate")$fitted)
#   data <- ddply(data, .(set), transform, sd2 = sd2-mean(sd2)+.125)
}

for(j in 1:10){
  data <- get4(sample(1:6, 6))

  # Line Segment plot
  qplot(data=data, x=x, xend=x, y=ynorm-.5*sd2, yend=ynorm+.5*sd2, geom="segment") + 
    facet_wrap(~set) + theme_stimuli() 
  ggsave(filename=paste("ZigZagSegment_Pic_", j, ".png", sep=""), path="./stimuli/ZigZags/", width=8, height=6, units="in")
  
  key <- ddply(data, .(set, correct), summarise, max.ratio.slope=max(abs(slope))/min(abs(slope)), max.ratio.sd=max(sd2)/min(sd2))
  write.table(key, paste("./stimuli/ZigZags/ZigZags_Table_", j, ".txt", sep=""), sep="\t", row.names=FALSE)
  write.csv(data, paste("./stimuli/ZigZags/ZigZags_Data_", j, ".csv", sep=""), row.names=FALSE)
}