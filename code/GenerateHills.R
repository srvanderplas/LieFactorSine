require(ggplot2)
require(plyr)
require(reshape2)
source("./code/themeStimuli.R")
library(e1071)

logit <- function(x, a){
  data <- data.frame(y=sigmoid(a*x), deriv=dsigmoid(a*x), deriv2=d2sigmoid(a*x))
  data
}

# Logistic approx. to step function I(x0[1] < x < x0[2])
logitcombo <- function(x, a, x0){
  logit(x+x0[1], a=a) - logit(x+x0[2], a=a)
}

# generates random start,end points for logistic approximations to step functions.
getrange <- function(lb=-2.5, ub=2.5){
  lower <- sample(seq(from=lb, to=mean(c(lb, ub)), length=10), 1)
  upper <- sample(seq(lower, ub, length=10), 1)
  c(lower, upper)
}

correct <- function(data, ell=.5){
  fp <- data$deriv
  f2p <- data$deriv2
  lambdap <- (sqrt((fp^2+1)^2-f2p*fp^2*ell) + fp^2 + 1)^-1    
  lambdam <- -(sqrt((fp^2+1)^2+f2p*fp^2*ell) + fp^2 + 1)^-1    
  
  data$upper <- data$y + (4*abs(lambdap)*sqrt(1+fp^2))^-1
  data$lower <- data$y - (4*abs(lambdam)*sqrt(1+fp^2))^-1  
  data$upper.cor <- + (4*abs(lambdap)*sqrt(1+fp^2))^-1
  data$lower.cor <- - (4*abs(lambdam)*sqrt(1+fp^2))^-1
  data
}

x <- seq(-3, 3, .1)
coef <- c(4, -2, -3, 1)
bumps <- rbind(getrange(lb=-1, ub=1), getrange(lb=-3, ub=3), getrange(lb=0, ub=3), getrange(lb=-2, ub=1))
# Uncorrected Data
data <- cbind(set=1, x=x, coef[1]*logitcombo(x, a=1, x0=bumps[1,]) + 
                coef[2]*logitcombo(x, a=3, x0=bumps[2,]) + 
                coef[3]*logitcombo(x, a=5, x0=bumps[3,]) + 
                coef[4]*logitcombo(x, a=3, x0=bumps[4,]))
data$upper <- data$y + .5
data$lower <- data$y - .5
data$upper.cor <- .5
data$lower.cor <- -.5

# Correct using the right slope, etc.
data2 <- correct(data, ell=.5)
data2$set <- 2

# Correct using slope assuming that "steepness" (a = 3) -- both over and undercorrect
data3 <- cbind(set=3, x=x, coef[1]*logitcombo(x, a=3, x0=bumps[1,]) + 
                coef[2]*logitcombo(x, a=3, x0=bumps[2,]) + 
                coef[3]*logitcombo(x, a=3, x0=bumps[3,]) + 
                coef[4]*logitcombo(x, a=3, x0=bumps[4,]))
data3 <- correct(data3, ell=.5)
data3$y <- data$y
data3$upper <- data3$y + data3$upper.cor
data3$lower <- data3$y + data3$lower.cor

# Correct using slope assuming that "steepness" (a = 1) -- Undercorrect
data4 <- cbind(set=4, x=x, coef[1]*logitcombo(x, a=1, x0=bumps[1,]) + 
                 coef[2]*logitcombo(x, a=1, x0=bumps[2,]) + 
                 coef[3]*logitcombo(x, a=1, x0=bumps[3,]) + 
                 coef[4]*logitcombo(x, a=1, x0=bumps[4,]))
data4 <- correct(data4, ell=.5)
data4$y <- data$y
data4$upper <- data4$y + data4$upper.cor
data4$lower <- data4$y + data4$lower.cor

# Correct using slope assuming that "steepness" (a = 5) -- Overcorrect 
data5 <- cbind(set=5, x=x, coef[1]*logitcombo(x, a=5, x0=bumps[1,]) + 
                 coef[2]*logitcombo(x, a=5, x0=bumps[2,]) + 
                 coef[3]*logitcombo(x, a=5, x0=bumps[3,]) + 
                 coef[4]*logitcombo(x, a=5, x0=bumps[4,]))
data5 <- correct(data5, ell=.5)
data5$y <- data$y
data5$upper <- data5$y + data5$upper.cor
data5$lower <- data5$y + data5$lower.cor

# Correct using slope assuming that "steepness" (a = .5) -- Severely Undercorrect. 
data6 <- cbind(set=6, x=x, coef[1]*logitcombo(x, a=.5, x0=bumps[1,]) + 
                 coef[2]*logitcombo(x, a=.5, x0=bumps[2,]) + 
                 coef[3]*logitcombo(x, a=.5, x0=bumps[3,]) + 
                 coef[4]*logitcombo(x, a=.5, x0=bumps[4,]))
data6 <- correct(data6, ell=.5)
data6$y <- data$y
data6$upper <- data6$y + data6$upper.cor
data6$lower <- data6$y + data6$lower.cor

data.all <- rbind(data, data2, data3, data4, data5, data6)
data.all$correct <- which(data.all$set==2)
data.all$set <- sample(1:6, 6)[data.all$set]


qplot(data=data, x=x, y=y, geom="line")
qplot(data=data, x=x, y=deriv, geom="line")
qplot(data=data, x=x, y=deriv2, geom="line")

qplot(data=data.all, x=x, y=upper.cor, geom="line") + geom_line(aes(y=lower.cor)) + geom_line(aes(y=deriv), linetype=2) + facet_wrap(~set)

data.points <- ddply(data.all, .(set, x, y, deriv, deriv2, upper, lower, upper.cor, lower.cor, correct), summarise, ypoints=runif(7, lower, upper))

# Ribbon Plot
qplot(data=data.all, x=x, ymin=lower, ymax=upper, geom="ribbon") + coord_equal(ratio=1) + facet_wrap(~set, nrow=2) + theme_stimuli()

# Points
qplot(data=data.points, x=x, y=ypoints, geom="jitter", alpha=I(.5)) + coord_equal(ratio=1) + facet_wrap(~set, nrow=2) + theme_stimuli()

unique(data.points$set[data.points$correct])

# Size of Actual Length Change
qplot(data=data.all, x=x, y=upper-lower, geom="line") + facet_wrap(~set, nrow=2)

