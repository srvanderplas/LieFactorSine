require(ggplot2)
require(plyr)
require(reshape2)
source("./code/themeStimuli.R")
library(e1071)

logit <- function(x, a=1, b=0){
  data <- data.frame(y=sigmoid(a*x), deriv=dsigmoid(a*x), deriv2=d2sigmoid(a*x))
  data
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

correct <- function(data){
  fp <- data$deriv
  f2p <- data$deriv2
  ell <- .1
  lambdap <- (sqrt((fp^2+1)^2-f2p*fp^2*ell) + fp^2 + 1)^-1    
  lambdam <- -(sqrt((fp^2+1)^2+f2p*fp^2*ell) + fp^2 + 1)^-1    
  
  
  data$upper <- data$y + (4*abs(lambdap)*sqrt(1+fp^2))^-1
  data$lower <- data$y - (4*abs(lambdam)*sqrt(1+fp^2))^-1
  data
}


x <- seq(-8, 8, .1)

# a can be any set of 4 positive numbers.
a <- sample(c( 1, 1, 1, 4))
bumps <- rbind(getrange(), getrange(), getrange())
data <- do.call("rbind", lapply(1:4, function(i) 
          data.frame(x=x, set=i, correct(logitcombo(x=x, a=a[i], x0=bumps[1,]) + 
                       logitcombo(x=x, a=a[i], x0=bumps[2,]) + 
                       logitcombo(x=x, a=a[i], x0=bumps[3,])),
                     ans=which.max(a))))


data.points <- ddply(data, .(x, y, set, ans, deriv, deriv2), summarise, ypoints = runif(5, lower, upper))
data.points$resid <- data.points$ypoints-data.points$y
data.points1 <- melt(data.points, id.vars=c("x", "ans", "set", "deriv", "y"), value.vars=c("ypoints", "resid"), value.name=c("ypoints"), variable.name="type")

# Ribbon plot
qplot(data=data, x=x, ymin=-.25*(1+(deriv^2)), ymax=.25*(1+(deriv^2)), geom="ribbon") + facet_wrap(~set) + theme_stimuli()

# Points
qplot(data=data.points, x=x, y=ypoints, geom="point", alpha=I(.5)) + 
  geom_line(data=data.points, aes(x=x,  y=y)) + facet_wrap(~set) + theme_stimuli()

qplot(data=data.points1, x=x, y=ypoints, geom="point", alpha=I(.5)) + facet_grid(type~set) + theme_stimuli()
qplot(data=data.points1, x=x, y=1/(1+deriv^2), geom="line") + facet_wrap(~set) #+ theme_stimuli()

data$correct[1]