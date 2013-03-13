# Code to generate 2x2 grid of pictures, 3 or 4 of which have similar or identical slopes.

# Input number of line segments and a vector of nonnegative slope magnitudes. 
# Output a dataset of x, y, slope, and segment number
getDataSet <- function(nsegment=4, slopes=1)
{
  segmentLength=2:(nsegment+1)
  initialSlope <- rbind(rep(c(-1, 1), times=ceiling(nsegment/2))[1:nsegment], 
                        rep(c(1, -1), times=ceiling(nsegment/2))[1:nsegment])
  x <- seq(1, sum(segmentLength), 1)
  genY <- data.frame(slope = sample(slopes, nsegment, replace=TRUE)*initialSlope[sample(c(1, 2), 1), ], 
                     times = sample(segmentLength, 4, replace=FALSE))
  slope <- unlist(sapply(1:4, function(i) rep(genY$slope[i], times=genY$times[i])))
  segment <- unlist(sapply(1:4, function(i) rep(i, times=genY$times[i])))
  data <- data.frame(x=x, y=cumsum(slope), slope=slope, segmentLength=segment)
  data
}



