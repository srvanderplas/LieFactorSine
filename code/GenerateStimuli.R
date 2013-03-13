# Code to generate 2x2 grid of pictures, 3 or 4 of which have similar or identical slopes.
require(ggplot2)
require(plyr)
require(reshape2)

# Input number of line segments and a vector of nonnegative slope magnitudes. 
# Output a dataset of x, y, slope, and segment number
getDataSet <- function(nsegment=4, slope=c(1))
{
  segmentLength=2:(nsegment+1)
  initialSlope <- rbind(rep(c(-1, 1), times=ceiling(nsegment/2))[1:nsegment], 
                        rep(c(1, -1), times=ceiling(nsegment/2))[1:nsegment])
  x <- seq(1, sum(segmentLength), 1)
  slope <- rep(slope, times=ceiling(nsegment/length(slope)))
  genY <- data.frame(slope = sample(slope, nsegment, replace=FALSE)*initialSlope[sample(c(1, 2), 1), ], 
                     times = sample(segmentLength, nsegment, replace=FALSE))
  slope <- unlist(sapply(1:nsegment, function(i) rep(genY$slope[i], times=genY$times[i])))
  segment <- unlist(sapply(1:nsegment, function(i) rep(i, times=genY$times[i])))
  data <- data.frame(x=x, y=cumsum(slope), slope=slope, segmentLength=segment)
  data$ynorm <- (data$y-mean(data$y))/sd(data$y)
  data
}


theme_stimuli <- function(base_size = 12, base_family = ""){
                   theme_grey(base_size = base_size, base_family = base_family) +
                   theme(legend.background = element_blank(), 
                         legend.key = element_blank(), 
                         panel.background = element_rect(fill="white",color="black"), 
                         panel.border = element_blank(), 
                         strip.background = element_rect(fill="grey",color="black"), 
                         plot.background = element_blank(),
                         panel.grid.major = element_line(color="grey75"),
                         axis.text=element_blank())
                  }

slopes <- list(c(1), c(1), c(1), c(1.5, 1))
data <- do.call("rbind", lapply(sample(1:4, 4), function(i) cbind(getDataSet(8, slope=unlist(slopes[i])), set=i)))


qplot(data=data, x=x, ymin=ynorm-.5, ymax=ynorm+.5, geom="ribbon") + facet_wrap(~set) + theme_stimuli()