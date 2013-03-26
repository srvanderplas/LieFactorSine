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
  lower <- sample(seq(from=lb, to=-1, length=10), 1)
  upper <- sample(seq(1, ub, length=10), 1)
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

getdata <- function(){
    x <- seq(-3, 3, length=40)
    coef <- c(3, 6)
    bumps <- rbind(getrange(lb=-2.5, ub=0), getrange(lb=0, ub=2.5))
    # Uncorrected Data
    data <- cbind(set=1, x=x, coef[1]*logitcombo(x, a=1, x0=bumps[1,]) + 
                    coef[2]*logitcombo(x, a=3, x0=bumps[2,]))
    data$upper <- data$y + .5
    data$lower <- data$y - .5
    data$upper.cor <- .5
    data$lower.cor <- -.5
    
    # Correct using the right slope, etc.
    data2 <- correct(data, ell=.5)
    data2$set <- 2
    
    # Correct using slope assuming that "steepness" (a = mean(coef)) -- both over and undercorrect
    data3 <- cbind(set=3, x=x, coef[1]*logitcombo(x, a=mean(coef), x0=bumps[1,]) + 
                    coef[2]*logitcombo(x, a=mean(coef), x0=bumps[2,]))
    data3 <- correct(data3, ell=.5)
    data3$y <- data$y
    data3$upper <- data3$y + data3$upper.cor
    data3$lower <- data3$y + data3$lower.cor
    
    # Correct using slope assuming that "steepness" (a = min(coef)) -- Undercorrect
    data4 <- cbind(set=4, x=x, coef[1]*logitcombo(x, a=min(coef), x0=bumps[1,]) + 
                     coef[2]*logitcombo(x, a=min(coef), x0=bumps[2,]))
    data4 <- correct(data4, ell=.5)
    data4$y <- data$y
    data4$upper <- data4$y + data4$upper.cor
    data4$lower <- data4$y + data4$lower.cor
    
    # Correct using slope assuming that "steepness" (a = max(coef)) -- Overcorrect 
    data5 <- cbind(set=5, x=x, coef[1]*logitcombo(x, a=max(coef), x0=bumps[1,]) + 
                     coef[2]*logitcombo(x, a=max(coef), x0=bumps[2,]))
    data5 <- correct(data5, ell=.5)
    data5$y <- data$y
    data5$upper <- data5$y + data5$upper.cor
    data5$lower <- data5$y + data5$lower.cor
    
    # Correct using slope assuming that "steepness" (a = .5) -- Severely Undercorrect. 
    data6 <- cbind(set=6, x=x, coef[1]*logitcombo(x, a=.5*min(coef), x0=bumps[1,]) + 
                     coef[2]*logitcombo(x, a=.5*min(coef), x0=bumps[2,]))
    data6 <- correct(data6, ell=.5)
    data6$y <- data$y
    data6$upper <- data6$y + data6$upper.cor
    data6$lower <- data6$y + data6$lower.cor
    
    data.all <- rbind(data, data2, data3, data4, data5, data6)
    
    avgcor <- mean(c(data.all$upper.cor, abs(data.all$lower.cor)))
    data.all$upper.cor[which(data.all$set==1)] <- avgcor
    data.all$lower.cor[which(data.all$set==1)] <- -avgcor
    data.all$upper[which(data.all$set==1)] <- data$y+avgcor
    data.all$lower[which(data.all$set==1)] <- data$y-avgcor
    
    data.all$type <- c("Untransformed", "Correct Transformation", 
                       "Transform for a=3", "Transform for a=1", 
                       "Transform for a=6", "Transform for a=.5")[data.all$set]
    
    data.all$correct <- data.all$set==2
    data.all$set <- sample(1:6, 6)[data.all$set]
    data.all
}

data.all <- getdata()
# 
# qplot(data=subset(data.all, type=="Untransformed"), x=x, y=y, geom="line")
# qplot(data=subset(data.all, type=="Untransformed"), x=x, y=deriv, geom="line")
# qplot(data=subset(data.all, type=="Untransformed"), x=x, y=deriv2, geom="line")
# 
# qplot(data=data.all, x=x, y=upper.cor, geom="line") + geom_line(aes(y=lower.cor)) + geom_line(aes(y=deriv), linetype=2) + facet_wrap(~set)

# Line Plot
qplot(data=data.all, x=x, xend=x, y=lower, yend=upper, geom="segment") + coord_equal(ratio=1) + facet_wrap(~set, nrow=2) + theme_stimuli()
ggsave(filename=paste("HillsLine_Pic_", i, ".png", sep=""), path="./stimuli/Hills/", width=8, height=6, units="in")

key <- ddply(data.all, .(set, correct, type), summarize, maxlength=max(upper-lower), meanlength=mean(upper-lower))

unique(data.all$set[data.all$correct])

