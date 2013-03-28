if(sum(grepl("multicore", installed.packages()[,1]))>0) library(multicore) else mclapply <- lapply
library(ggplot2)
library(reshape2)
library(plyr)

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

trans3d <- function(x,y,z, plot) {
  tmat <- t(cbind(x,y,z,1) %*% plot)
  list(x=tmat[1,]/tmat[4,], y=tmat[2,]/tmat[4,])
} 

x <- seq(0, 2*pi, length=42)[2:41]
data <- do.call("rbind", lapply(seq(-.5, .5, 1), function(i) data.frame(x=x, y=2*sin(x), z=i)))

data.persp <- acast(data, x~z, value.var="y")
x <- sort(unique(data$x))
y <- sort(unique(data$y))
z <- sort(unique(data$z))

#' Using the perspective plot, we can get the length of the lines in perspective space and 
#' look at the variation in line length for various phi and d values. 
plot.opts <- expand.grid(phi=seq(-90, 90, 15), d=seq(1, 50, length=10))
data2 <- rbind.fill(mclapply(1:nrow(plot.opts), function(k){
  p <- persp(x, z, data.persp, xlab="", ylab="", zlab="", theta=0, phi=plot.opts$phi[k], border="black", 
             shade=.35, col="white", xlim=c(-pi/12, 2*pi+pi/12), ylim=c(-1.75, 1.75), scale=FALSE, box=FALSE, d=plot.opts$d[k]) # , ltheta=0, lphi=-15
  
  data$transx <- trans3d(data$x, data$z, data$y, p)$x
  data$transy <- trans3d(data$x, data$z, data$y, p)$y
  data.frame(phi=plot.opts$phi[k], d=plot.opts$d[k], ddply(data, .(x), summarize, transdist=sqrt(diff(transx)^2 + diff(transy)^2)))
}))

data3 <- ddply(data2, .(phi, d),transform, transdist=transdist, transdist.center=transdist-mean(transdist, na.rm=TRUE), transratio=max(transdist)/min(transdist))

#' Actual line distance in 3d projection
qplot(data=data3, x=x, y=transdist, geom="line", colour=d, group=d) + facet_wrap(~phi, scales="free_y")
#' Centered line distance (i.e. effect of  x|d, phi)
qplot(data=data3, x=x, y=transdist.center, geom="line", colour=d, group=d) + facet_wrap(~phi)


#' Ratio of max/min line distance for 3d projection, over multiple perspective phi angles and d values
qplot(data=subset(data3, phi!=0), x=abs(phi), y=transratio, geom="line", colour=d, group=d)



#' Examining changes in phi on the length and length ratio
plot.opts <- expand.grid(phi=seq(5, 90, 5), d=c(.5, 1, 2, 5, 10, 15, 25, 50, 100))
data2 <- rbind.fill(mclapply(1:nrow(plot.opts), function(k){
  p <- persp(x, z, data.persp, xlab="", ylab="", zlab="", theta=0, phi=plot.opts$phi[k], border="black", 
             shade=.35, col="white", xlim=c(-pi/12, 2*pi+pi/12), ylim=c(-1.75, 1.75), scale=FALSE, box=FALSE, d=plot.opts$d[k]) # , ltheta=0, lphi=-15
  
  data$transx <- trans3d(data$x, data$z, data$y, p)$x
  data$transy <- trans3d(data$x, data$z, data$y, p)$y
  data.frame(phi=plot.opts$phi[k], d=plot.opts$d[k], ddply(data, .(x), summarize, transdist=sqrt(diff(transx)^2 + diff(transy)^2)))
}))
    
data3 <- ddply(data2, .(phi, d),transform, transdist=transdist, transdist.center=transdist-mean(transdist, na.rm=TRUE), transratio=max(transdist)/min(transdist))

#' Actual line distance in 3d projection
qplot(data=data3, x=x, y=transdist, geom="line", colour=d, group=d) + facet_wrap(~phi, scales="free_y")
#' Centered line distance (i.e. effect of  x|d, phi)
qplot(data=data3, x=x, y=transdist.center, geom="line", colour=d, group=d) + facet_wrap(~phi)


#'-----------------------------------------------------------------------------------------
#' Conclusions 
#'   - phi=60 seems to be the best rotation angle for showing the line difference effect
#'   - phi=45 shows even line distances, but also connects back to the illusion 
#'     (i.e. it's visible if shading is off)
#'     
#'-----------------------------------------------------------------------------------------

p <- persp(x, z, data.persp, xlab="", ylab="", zlab="", theta=0, phi=45, border="black", 
           shade=.35, col="white", xlim=c(-pi/12, 2*pi+pi/12), ylim=c(-1.75, 1.75), 
           scale=FALSE, box=FALSE, d=500, expand=3/(pi)) # , ltheta=0, lphi=-15

#' If we back-transform the 3d transformation (as it's a linear transformation, we can do 
#' this with lm()) we can transform the scaled persp() values back into our sine-transformation
#' values, which will allow us to compare 3d to our actual transformation. 

data$transx <- trans3d(data$x, data$z, data$y, p)$x
data$transy <- trans3d(data$x, data$z, data$y, p)$y
data.2d <- ddply(data, .(x,y), summarise, xstart=min(transx), xend=max(transx), ystart=min(transy), yend=max(transy), len=sqrt(diff(transx)^2+diff(transy)^2))

#'-----------------------------------------------------------------------------------------
#' Linear relationship between x and xstart, xend
# qplot(data=data.2d, x=x, y=xstart, geom="line") + geom_line(aes(y=xend))
#' Linear relationship between y and ystart, yend
# qplot(data=data.2d, x=y, y=ystart, geom="line") + geom_line(aes(y=yend))
#'-----------------------------------------------------------------------------------------

data.2d$xstart <- predict(lm(data=data.2d, x~xstart))
data.2d$xend <- predict(lm(data=data.2d, x~xend))
data.2d$ystart <- lm(data=data.2d, y~ystart)$coefficients[2]*data.2d$ystart
data.2d$yend <- lm(data=data.2d, y~yend)$coefficients[2]*data.2d$yend

trans.resid <- data.frame(x = data.2d$x,
  xresid1 = data.2d$xstart-createSine(40, 1, f=f, fprime=fprime, f2prime)$xstart,
  xresid2 = data.2d$xend-createSine(40, 1, f=f, fprime=fprime, f2prime)$xend,
  yresid1 = data.2d$ystart-createSine(40, 1, f=f, fprime=fprime, f2prime)$ystart,
  yresid2 = data.2d$yend-createSine(40, 1, f=f, fprime=fprime, f2prime)$yend
  )

trans.resid.long <- melt(trans.resid, id.vars=1, value.name="resid", variable.name="type")
trans.resid.long$axis <- substr(trans.resid.long$type, 1, 1)
trans.resid.long$end <- factor(as.numeric(substr(trans.resid.long$type, 7, 7)), 
                               labels=c("start", "end"), ordered=TRUE)
qplot(data=trans.resid.long, x=x, y=resid, geom="line", linetype=end) + 
  facet_wrap(~axis) + ylab("Difference in endpoints, Persp-createSine()") + 
  scale_linetype_discrete("Endpoint") + 
  scale_x_continuous(breaks=seq(0, 2*pi, by=pi/2), 
                     labels=c("0", expression(paste(pi,"/2")), expression(pi), 
                              expression(paste("3",pi, "/2")), expression(paste("2",pi))))

#' Very little difference with d=500 and original createSine() data. 
#' Slight periodicity may be useful as a correction factor later?