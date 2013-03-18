library(ggplot2)
library(reshape2)
trans3d <- function(x,y,z, plot) {
  tmat <- t(cbind(x,y,z,1) %*% plot)
  list(x=tmat[1,]/tmat[4,], y=tmat[2,]/tmat[4,])
} 


qplot(x=x, xend=xend, y = ystart, yend=yend, geom="segment", data=createSine(40, 1, f=f, fprime=fprime, f2prime), colour=I("white"), size=I(2)) +
  geom_line(aes(x=x, y=ystart), colour="white") + 
  geom_line(aes(x=x, y=yend), colour="white") + 
  theme(panel.grid.major=element_blank(), panel.background = element_rect(fill = "black", colour = "black"),
        panel.grid.minor=element_blank(), panel.background=element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), 
        axis.text = element_blank()) + coord_equal(ratio=1)


x <- seq(-6, 6, .1)
y <- cumsum(sign(cos(x))*rep(.1, length(x)))
qplot(x,y, geom="line") + geom_line(aes(x=x, y=y+1))

x <- seq(0, 2*pi, length=40)
data <- do.call("rbind", lapply(seq(-.5, .5, .5), function(i) data.frame(x=x, y=2*sin(x), z=i)))

qplot(data=data, x, y, geom="line", colour=z, group=z) + scale_colour_continuous(guide="none", high="grey30", low="grey60") + 
  geom_segment(data=createSine(40, 1, f=f, fprime=fprime, f2prime, a=0, b=2*pi), aes(x=x, y=ystart, xend=xend, yend=yend), inherit.aes=FALSE)

data.persp <- acast(data, x~z, value.var="y")
x <- unique(data$x)
y <- unique(data$z)


persp(x, y, data.persp, ylab="z", zlab="y", theta=-10, phi=35, border=NA, shade=.75) 
# Extant width makes perfect sense here - it's exactly what you'd be judging if you were judging width for a perspective view, 
# except that you would have more information and/or be able to assume some sort of width constancy behavior. 
# It's just when we look at purely 2d things that the extant width doesn't make sense anymore - because there is no 3d context. 


persp(x, y, data.persp, ylab="z", zlab="y", theta=-10, phi=30, border=NA, shade=.5) 
# Note how ambiguous this figure could be - is the "wave" in z or y? is 0-pi the "near" region, or is it further away with pi-2*pi as the "near" region? It also seems to "flip" like a necker cube, at least to me.

persp(x, y, data.persp, ylab="z", zlab="y", theta=-10, phi=30, shade=.5) 
# Adding in the lines helps somewhat, but still doesn't completely resolve it.


qplot(data=createSine(40, 1, f=f, fprime=fprime, f2prime, a=0, b=2*pi), x=x, xend=x, y=ystart, yend=yend, geom="segment") + 
  geom_line(aes(y=ystart), linetype=2) + geom_line(aes(y=yend), linetype=2) + 
  geom_text(aes(y=-2.6, x=3*pi/2, label="Illusory Contour")) + 
  geom_text(aes(y=2.6, x=pi/2, label="Illusory Contour")) + 
  theme(panel.grid.major=element_blank(), panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.minor=element_blank(), panel.background=element_blank(),
        axis.title = element_blank(), axis.ticks = element_blank(), 
        axis.text = element_blank()) + coord_equal(ratio=1)


persp(x, y, data.persp, ylab="z", zlab="y", theta=0, phi=90, border="grey20") 
persp(x, y, -data.persp, ylab="z", zlab="y", theta=0, phi=90, border="grey20") 

