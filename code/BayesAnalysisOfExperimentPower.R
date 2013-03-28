library(ggplot2)
library(plyr)
library(reshape2)
w <- seq(0, 1.3, length=6)
theta <- .75

f <- function(y, theta, w){ # the yth largest response in w (the index, not the value of w)
  p <- function(i) 1/(.01+(w-i)^2)
  sapply(theta, function(i) p(i)[y]/sum(p(i)))*1/length(theta)
}

sum(f(5, theta=seq(0, 1.3, .01), w))


get.post.w <- function(y, theta, w){
  apply(sapply(y, f, theta=theta, w=w), 1, function(x) sum(log(x)))
}

get.post.all <- function(y, theta, w, pic.num){
  dat <- ldply(pic.num, function(i) data.frame(theta=theta, f=get.post.w(y=y[i,], theta, w=w[i,]), pic.num=i))
  dat
}
f(4, theta, w)

theta <- seq(0, 1.3, .01)
data <- ldply(1:6, function(i) data.frame(theta=theta, f=f(i, theta, w), y=i))

qplot(data=data, x=theta, y=f, colour=y, geom="line", group=y)

y <- c(4, 4, 3, 4, 4, 5, 4)

qplot(x=theta, y=get.post(y, theta, w), geom="line")
theta[which.max(get.post(y, theta, w))]


w <- rbind(seq(0, 1.3, len=6), seq(0, .5, len=6), seq(.4, .9, len=6), seq(.8, 1.3, len=6), 
           seq(0, .25, len=6), seq(.2, .45, len=6), seq(.4, .65, len=6), seq(.6, .85, len=6), seq(.8, 1.05, len=6), seq(1, 1.25, len=6))
y <- c(4, 6, 5, 1, 6, 6, 6, 4, 1, 1)


data <- NULL
for(th in seq(.1, 1.2, by=.1)){
  for(z in 1:13){
    w <- matrix(c(  0, .2,  .4,  .6,   .8,  1.0,      .1,  .3,  .5,  .7,  .9,  1.1,
                  .05, .3,  .5,  .65,  .8,  1.0,      .1,  .3,  .55, .7,  .85, 1.0,
                   .4, .6,  .7,  .8,   .9,  1.05,     .35, .65, .75, .85, .95, 1.05,
                  .35, .5,  .6,  .7,   .8,   .95,     .4,  .55, .65, .75, .85, 1.0,
                   .5, .65, .75, .8,   .9,  1.0,      .5,  .6,  .7,  .75, .85, 1.0), nrow=12, ncol=6, byrow=TRUE)
    pic.num=1:10
    y <- t(sapply(1:10, function(i) sample(1:6, 6, prob=f(1:6, th, w[i,]), replace=TRUE)))
    data1 <- ddply(get.post.all(y, theta, w, pic.num), .(theta), summarise, f=exp(sum(f)+350))
    data1$method <- 1
    data1$f <- data1$f/sum(data1$f)
    
    w <- t(sapply(1:10, function(i) sort(sample(seq(0, 1.3, .05), 6))))
    y <- t(sapply(1:10, function(i) sample(1:6, 6, prob=f(1:6, th, w[i,]), replace=TRUE)))
    data2 <- ddply(get.post.all(y, theta, w, pic.num), .(theta), summarise, f=exp(sum(f)+350))
    data2$method <- 2
    data2$f <- data2$f/sum(data2$f)
    
    data.z <- rbind(data1, data2)
    data.z$z <- z
    data.z$th <- th
    data <- rbind(data, data.z)
  }
}
qplot(data=data, x=theta, y=f, geom="line", colour=factor(method), fill=factor(method)) + scale_colour_manual(values=c("blue", "green"))  + scale_fill_manual(values=c("blue", "green")) + facet_grid(z~th) + geom_vline(aes(xintercept=th), alpha=.1)

# Conclude: Method 1 and Method 2 are approximately the same, and method 1 requires less coding for the Turk site.