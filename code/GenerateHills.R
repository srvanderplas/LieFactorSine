require(ggplot2)
require(plyr)
require(reshape2)
source("./code/themeStimuli.R")


logit <- function(x, a=1, b=0){
  1/(1+exp(a*x + b))
}

# Logistic approx. to step function I(x0[1] < x < x0[2])
logitcombo <- function(x, a=1, x0){
  logit(x+x0[1], a=a) - logit(x+x0[2], a=a)
}

# generates random start,end points for logistic approximations to step functions.
getrange <- function(){
  lower <- sample(seq(-8, 5, 1), 1)
  upper <- sample(seq(lower, 8, 1), 1)
  c(lower, upper)
}


x <- seq(-8, 8, .1)

# a can be any set of 4 positive numbers.
a <- sample(c( 1, 1, 1, 2))
data <- data.frame(
        sapply(1:4, function(i) logitcombo(x=x, a=a[i], x0=getrange())) +
        sapply(1:4, function(i) logitcombo(x=x, a=a[i], x0=getrange())) +
        sapply(1:4, function(i) logitcombo(x=x, a=a[i], x0=getrange())), 
        correct=which.max(a))
data1 <- data.frame(x=x, data)
data <- melt(data1, id.vars=c("x", "correct"), value.name="y", variable.name="type")
data$set <- as.numeric(data$type)

data.points <- ddply(data, .(x, y, set, type, correct), summarise, ypoints = rnorm(5, y, .1))

# Ribbon plot
qplot(data=data, x=x, ymin=y-.25, ymax=y+.25, geom="ribbon") + facet_wrap(~set) + theme_stimuli()

# Points
qplot(data=data.points, x=x, y=ypoints, geom="point", alpha=I(.5)) + 
  geom_line(data=data.points, aes(x=x,  y=y)) + facet_wrap(~set) + theme_stimuli()

data$correct[1]