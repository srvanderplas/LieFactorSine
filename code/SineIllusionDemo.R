library(ggplot2)
library(reshape2)
trans3d <- function(x,y,z, plot) {
  tmat <- t(cbind(x,y,z,1) %*% plot)
  list(x=tmat[1,]/tmat[4,], y=tmat[2,]/tmat[4,])
} 

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


getSecantSegment <- function(x0, df, f, fprime, f2prime){
  ell     <- sapply(x0, function(i) df$ell[which.min(abs(i-df$x))]/2)
  
  dy <- diff(range(df$y))
  dx <- diff(range(df$x))
  a <- dx/(dy + 2*ell) 
  
  fp <- a*fprime(x0)
  f2p <- a*f2prime(x0)
  lambdap <- (sqrt((fp^2+1)^2-f2p*fp^2*ell) + fp^2 + 1)^-1    
  lambdam <- -(sqrt((fp^2+1)^2+f2p*fp^2*ell) + fp^2 + 1)^-1    
  
  #---- Approximation
  
  x2 <- lambdap*fprime(x0)+x0
  x1 <- lambdam*fprime(x0)+x0
  y2 <- f(x0)-lambdap
  y1 <- f(x0)-lambdam
  #----
  
  df2 <- data.frame(x=x0, y=f(x0), deriv=fprime(x0),
                    sec.xstart=x1, sec.xend = x2, 
                    sec.ystart=y1, sec.yend = y2,
                    ell = 2*ell)
  
  df2$sec.ellp <- (4*abs(lambdap)*sqrt(1+fp^2))^-1
  df2$sec.ellm <- (4*abs(lambdam)*sqrt(1+fp^2))^-1
  #   df2$sec.ellp <- with(df2, sqrt((sec.yend-y)^2+(sec.xend-x)^2))
  #   df2$sec.ellm <- with(df2, sqrt((y-sec.ystart)^2+(x-sec.xstart)^2))
  df2$type <- "Perceived Width"
  df2$a <- a
  return(df2)
}

f <- function(x) 2*sin(x)
fprime <- function(x) 2*cos(x)
f2prime <- function(x) -2*sin(x)

# Attempt at making "self-luminous" images as in Gregory (1968) "Perceptual Illusions and Brain Models." Proc Roy. Soc. B 171 (279-296).
# Not sure what's required for "self luminous" images, but this doesnt really fix it for me... 
qplot(x=x, xend=xend, y = ystart, yend=yend, geom="segment", data=createSine(40, 1, f=f, fprime=fprime, f2prime), colour=I("white"), size=I(2)) +
  geom_line(aes(x=x, y=ystart), colour="white") + 
  geom_line(aes(x=x, y=yend), colour="white") + 
  theme(panel.grid.major=element_blank(), panel.background = element_rect(fill = "black", colour = "black"),
        panel.grid.minor=element_blank(), panel.background=element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), 
        axis.text = element_blank()) + coord_equal(ratio=1)


# 3-dimensionality seems to be a function of the curve mismatch or the gradient (or both) - it's not evident here.
x <- seq(-6, 6, .1)
y <- cumsum(sign(cos(x))*rep(.1, length(x)))
qplot(x,y, geom="line") + geom_line(aes(x=x, y=y+1))

x <- seq(0, 2*pi, length=40)
data <- do.call("rbind", lapply(seq(-.5, .5, .5), function(i) data.frame(x=x, y=2*sin(x), z=i)))

# Two-dimensional plot, where y_plot = y+z
qplot(data=data, x, y=y+z, geom="line", colour=z, group=z) + scale_colour_continuous(guide="none", high="grey30", low="grey60") + 
  geom_segment(data=createSine(40, 1, f=f, fprime=fprime, f2prime, a=0, b=2*pi), aes(x=x, y=ystart, xend=xend, yend=yend), inherit.aes=FALSE) + 
  theme(panel.grid.major=element_blank(), panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor=element_blank(), panel.background=element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), 
        axis.text = element_blank()) + coord_equal(ratio=1)

data.persp <- acast(data, x~z, value.var="y")
x <- sort(unique(data$x))
y <- sort(unique(data$y))
z <- sort(unique(data$z))


persp(x, z, data.persp, ylab="z", zlab="y", theta=-10, phi=35, border=NA, shade=.75) 
# Extant width makes perfect sense here - it's exactly what you'd be judging if you were judging width for a perspective view, 
# except that you would have more information and/or be able to assume some sort of width constancy behavior. 
# It's just when we look at purely 2d things that the extant width doesn't make sense anymore - because there is no 3d context. 


persp(x, z, data.persp, ylab="z", zlab="y", theta=-7, phi=30, border=NA, shade=.5) 
# Note how ambiguous this figure could be - is the "wave" in z or y? is 0-pi the "near" region, or is it further away with pi-2*pi as the "near" region? It also seems to "flip" like a necker cube, at least to me.

persp(x, z, data.persp, ylab="z", zlab="y", theta=-7, phi=30, shade=.5) 
# Adding in the lines helps somewhat, but still doesn't completely resolve it.


x <- seq(0, 2*pi, length=40)
data2 <- do.call("rbind", lapply(x, function(i) data.frame(x=i, y=c(-.5, .5), z=2*sin(i) +c(-.5, .5))))

data2.persp <- acast(data2, x~y, value.var="z")
x2 <- sort(unique(data2$x))
y2 <- sort(unique(data2$y))
z2 <- sort(unique(data2$z))

# This is what pops out at me when things "flip" (the right side, at least...)
persp(x2, y2, data2.persp, theta=0, phi=0, border=NA, shade=.5) 


# Flipping 
persp(x, z, data.persp, ylab="z", zlab="y", theta=0, phi=90, border="grey20")  # A
persp(x, z, -data.persp, ylab="z", zlab="y", theta=0, phi=90, border="grey20") # B
persp(x, z, data.persp, ylab="z", zlab="y", theta=10, phi=0, border="grey20")  # C
persp(x, z, -data.persp, ylab="z", zlab="y", theta=0, phi=0, border="grey20")  # D
# I seem to flip between B and C... but clearly you need the depth cues along the nearly-linear region to accurately parse the graph, particularly for C,D.

qplot(data=createSine(40, 1, f=f, fprime=fprime, f2prime, a=0, b=2*pi), x=x, xend=x, y=ystart, yend=yend, geom="segment") + 
  geom_line(aes(y=ystart), linetype=2) + geom_line(aes(y=yend), linetype=2) + 
  geom_text(aes(y=-2.6, x=3*pi/2, label="Illusory Contour")) + 
  geom_text(aes(y=2.6, x=pi/2, label="Illusory Contour")) + 
  theme(panel.grid.major=element_blank(), panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor=element_blank(), panel.background=element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), 
        axis.text = element_blank()) + coord_equal(ratio=1)

qplot(data=createSine(40, 1, f=f, fprime=fprime, f2prime, a=0, b=2*pi), x=x, xend=x, y=ystart, yend=yend, geom="segment", colour=I("grey70")) + 
  geom_line(aes(y=ystart), linetype=2) + geom_line(aes(y=yend), linetype=2) + 
  theme(panel.grid.major=element_blank(), panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor=element_blank(), panel.background=element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), 
        axis.text = element_blank()) + coord_equal(ratio=1)

library(plyr)
dframe <- createSine(n = 150, len = 1, f=f, fprime=fprime, f2prime)
dframe$ystartcts <- dframe$ystart
dframe$yendcts <- dframe$yend
dframe[1:150,c(2, 3, 5, 6)] <- NA
dframe[(1:15)*10-5, c(2, 3)] <- dframe[(1:15)*10-5, 1] 
dframe[(1:15)*10-5, 5] <- dframe[(1:15)*10-5, 4] - .5
dframe[(1:15)*10-5, 6] <- dframe[(1:15)*10-5, 4] + .5
dframe$type <- "Data"
dframe.1 <- getSecantSegment(dframe$xstart[!is.na(dframe$xstart)], dframe, f, fprime, f2prime)
names(dframe.1) <- c("x", "y", "deriv", "xstart", "xend", "ystart", "yend", "ell", "ell.quad1", "ell.quad2", "type", "a")
dframe.1$vangle <- with(dframe.1, atan(deriv))
dframe <- rbind.fill(dframe, dframe.1)


qplot(x=x, y=y, geom="line", data=dframe, colour=I("grey50")) + theme_bw() + 
  geom_line(aes(y=ystartcts), colour="grey50", linetype=4) + 
  geom_line(aes(y=yendcts), colour="grey50", linetype=4) +
  geom_segment(data=subset(dframe, !is.na(type) & type=="Perceived Width"), aes(x=xstart, xend = xend, y=ystart, yend=yend, colour=type))  + 
  coord_equal(ratio=1) + scale_colour_manual("", values=c("black", "blue")) + theme(legend.position="bottom")  + 
  scale_x_continuous(breaks=seq(0, 2*pi, by=pi/2), 
                     labels=c("0", expression(paste(pi,"/2")), expression(pi), expression(paste("3",pi, "/2")), expression(paste("2",pi))))

p <- persp(x, z, -data.persp, ylab="z", zlab="y", theta=0, phi=90, border="grey20") # B
lines(trans3d(x, .25*sin(x)-.375, 1+0*x, p), col="red")
lines(trans3d(x, .25*sin(x)+.125, 1+0*x, p), col="red")



