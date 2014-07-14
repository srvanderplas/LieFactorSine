setwd("/home/susan/Documents/R Projects/LieFactorSine/")
load("./code/BayesAnalysisW.Rdata")
require(ggplot2)
require(reshape2)
require(plyr)
require(msm)

turkdata <- read.csv("./data/turkdataclean.csv", stringsAsFactors=FALSE)

setwd("./code")
#-----------------------Distribution of W----------------------------------------------------------
logpost <- function(data, par){
  temp <- sum(dtnorm(data$answer.w, mean=par[1], sd=par[2], lower=0, upper=1.4, log=TRUE))
  temp <- temp-max(temp)*.9
}

get_posterior_density <- function(data, pars){
  temp <- sapply(1:nrow(pars), function(i) logpost(data, pars[i,]))
  temp <- exp(temp)/sum(exp(temp))
  data.frame(mean=pars[,1], sd=pars[,2], f=temp)
}

#--------------------Overall Marginals-------------------------------------------------------------------------
pars <- as.matrix(expand.grid(seq(0, 2, .01), seq(.15, .6, .01)))

overall <- ddply(turkdata, .(test_param), get_posterior_density, pars=pars)
overall.mean <- ddply(overall[,-3], .(test_param, mean), summarise, f=sum(f))
overall.mean <- ddply(overall.mean, .(test_param), transform, f=f/sum(f))

#' Posterior marginal distribution over individual and std. deviation
qplot(data=overall.mean, x=mean, y=f, geom="line", colour=test_param, group=test_param) + 
  scale_colour_discrete("Function Type") + xlim(c(0, 1.4)) +
  xlab("Mean Preferred Weighting") + ylab("Density") + theme_bw()  + theme(legend.position="bottom") 
# ggsave("figure/fig-OverallMeansW.pdf", width=4, height=4, units="in")

overall.sd <- ddply(overall[,-2], .(test_param, sd), summarise, f=sum(f))
overall.sd <- ddply(overall.sd, .(test_param), transform, f=f/sum(f))

#' Posterior marginal distribution over individual and mean
qplot(data=overall.sd, x=sd, y=f, geom="line", colour=factor(test_param), group=test_param) + theme_bw() + xlab("Posterior SD") + ylab("Density")

#' Posterior joint dist of mean, sd over individuals
#' since stat_density2d won't use weights, ... improvise!
overall.joint.sample <- sample(1:nrow(overall), size=50000, replace=TRUE, prob=overall$f)
ggplot(data=overall[overall.joint.sample,], aes(x=mean, y=sd)) + 
  stat_density2d(n=c(75, 40), geom="density2d", aes(colour=test_param)) + 
  facet_wrap(~test_param) + scale_colour_discrete(guide="none") + theme_bw() + 
  xlab("Mean Preferred Weighting") + ylab("Std. Deviation")
# ggsave("figure/fig-Joint2dDensityW.pdf", width=4, height=4, units="in")

#--------------------Individual Distribution of Theta ------------------------------------------------------------------

# test <- ddply(turkdata, .(ip.id, test_param), get_posterior_density, pars=pars)
# 
# test.mean <- ddply(test, .(ip.id, test_param, mean), summarise, f=sum(f))
# test.mean <- ddply(test.mean, .(ip.id, test_param), transform, f=f/sum(f))
# 
# participants <- dcast(ddply(turkdata, .(ip.id, test_param), summarise, n=length(test_param)), ip.id~test_param, value.var="n")
# ipsubset <- subset(participants, rowSums(is.na(participants))==0 & rowSums(participants[,2:4]>6, na.rm=TRUE)==3)$ip.id
# 
# par_labeller <- function(var, value){
#   n <- sapply(value, function(i) sum(subset(participants, ip.id%in%i)[,2:4]))
#   value <- paste("Participant ", as.character(value), "\n(n = ", n, ")", sep="")
#   return(value)
# }

#' Plot 4 individuals who did at least 6 figures of each trial 
qplot(data=subset(test.mean, ip.id%in%ipsubset), x=mean, y=f, group=test_param, colour=test_param, geom="line") + 
  facet_grid(.~ip.id, labeller=par_labeller) + scale_colour_discrete("Function Type") + theme_bw() + 
  theme(legend.position="bottom") + xlab("Mean Preferred Weighting") + ylab("Density")
# ggsave("figure/fig-IndivMeanAllFcnsW.pdf", width=7, height=3.5)      



#' Posterior mean estimates, including CI information for the individual MEAN 
#' (i.e. not for any individual observation)
# test.post.indiv<- ddply(test.mean, .(ip.id, test_param), 
#                         function(x){
#                           ex=sum(x$mean*x$f)
#                           n=nrow(subset(turkdata, turkdata$ip.id==ip.id[1] & turkdata$test_param==test_param[1]))
#                           samp <- matrix(sample(x$mean, n*11, prob=x$f, replace=TRUE), ncol=11)
#                           z <- as.numeric(quantile(rowMeans(samp), c(.025, .5, .975)))
#                           data.frame(ip.id=unique(x$ip.id), test_param=unique(x$test_param), lb=z[1], mean = ex, median=z[2], ub=z[3], n=n)
#                         })
# 
# overall.mean.f <- ddply(test.mean, .(test_param, mean), summarise, f=sum(f))
# overall.mean.f <- ddply(overall.mean.f, .(test_param), transform, f=f/sum(f))
# 
# overall.mean.bounds <- ddply(overall.mean.f, .(test_param), function(x){
#   ex=sum(x$mean*x$f)
#   n=length(unique(subset(turkdata, turkdata$test_param==test_param)$ip.id))
#   samp <- matrix(sample(x$mean, n*11, prob=x$f, replace=TRUE), ncol=11)
#   sample.mean = mean(samp)                          
#   sdev = sd(rowMeans(samp))
#   lb = as.numeric(quantile(rowMeans(samp), .025))
#   med = as.numeric(quantile(rowMeans(samp), .5))
#   ub = as.numeric(quantile(rowMeans(samp), .975))
#   data.frame(lb=lb, mean=sample.mean, median=med, ub=ub)
# })
# 
# test.post.indiv$functions <- c("Exponential", "Inverse", "Sine")[as.numeric(test.post.indiv$test_param)]
# test.post.indiv$functions <- factor(test.post.indiv$functions, levels=c("Sine", "Exponential", "Inverse"))
# overall.mean.bounds$functions <- c("Exponential", "Inverse", "Sine")[as.numeric(factor(overall.mean.bounds$test_param))]

qplot(data=test.post.indiv,  x=lb, xend=ub, y=ip.id, yend=ip.id, geom="segment", colour=test_param) + 
  facet_wrap(~functions) + geom_point(aes(x=median), colour="black") + 
  geom_vline(data=overall.mean.bounds, aes(xintercept=lb), linetype=3) + 
  geom_vline(data=overall.mean.bounds, aes(xintercept=median)) + 
  geom_vline(data=overall.mean.bounds, aes(xintercept=ub), linetype=3) + 
  ylab("Participant ID") + xlab("Mean Preferred Weighting") + theme_bw() + theme(legend.position="none") + 
  scale_colour_discrete("Function Type")
# ggsave("figure/fig-CIindivMeanW.pdf", width=6, height=6, units="in")

#' Posterior estimates, including CI information for the individual observations 
#' (i.e. not for any individual observation)
# indiv.value.bounds <- ddply(test.mean, .(ip.id, test_param), function(x){
#   lb=x$mean[which.min(abs(cumsum(x$f)-.025))]
#   med=x$mean[which.min(abs(cumsum(x$f)-.5))]
#   ub=x$mean[which.min(abs(cumsum(x$f)-.975))]
#   data.frame(lb=lb, median=med, ub=ub)
# })
# 
# overall.value.bounds <- ddply(overall.mean.f, .(test_param), function(x){
#   xnew <- sample(x$mean, length(x$mean), prob=x$f, replace=TRUE)
#   z <- as.numeric(quantile(xnew, c(.025, .5, .975)))
#   data.frame(lb=z[1], median=z[2], ub=z[3])
# })
# Posterior Distribution for theta without averaging over individuals
# qplot(data=overall.mean.f, x=mean, y=f, geom="line", colour=test_param) + 
#   xlab("Psychological Lie Factor\nEstimated Distribution for All Individuals") + 
#   theme_bw() + theme(legend.position="bottom") + scale_color_discrete("Function Type") + 
#   ylab("Density")
# 
# qplot(data=indiv.value.bounds,  x=lb, xend=ub, y=ip.id, yend=ip.id, geom="segment", colour=test_param) + 
#   facet_wrap(~test_param) + geom_point(aes(x=median), colour="black") + 
#   geom_vline(data=overall.value.bounds, aes(xintercept=lb), linetype=3) + 
#   geom_vline(data=overall.value.bounds, aes(xintercept=median)) + 
#   geom_vline(data=overall.value.bounds, aes(xintercept=ub), linetype=3) + 
#   ylab("Participant ID") + xlab("Lie Factor") + theme_bw() + theme(legend.position="bottom") + 
#   scale_colour_discrete("Function Type")
# 


#' Plot both individual user and group posterior estimates of preferred weightings. 
#' We may need more power/trials for each user in phase 2 if we want to do inference on w directly.
# test.mean.marginal <- ddply(test, .(ip.id, test_param, mean), summarise, f=sum(f))
# test.mean.marginal$f <- unlist(dlply(test.mean.marginal, .(ip.id, test_param), summarise, f=f/sum(f)))
# test.mean.marginal$functions <- c("Exponential", "Inverse", "Sine")[as.numeric(as.factor(test.mean.marginal$test_param))]
# test.mean.marginal$functions <- factor(test.mean.marginal$functions, levels=c("Sine", "Exponential", "Inverse"))
# overall.mean$functions <- c("Exponential", "Inverse", "Sine")[as.numeric(as.factor(overall.mean$test_param))]

ggplot(data=test.mean.marginal, aes(x=mean, y=f, group=ip.id, colour=test_param)) + geom_line(alpha=I(.175)) + 
  facet_wrap(~functions) + ylab("Density") + xlab("Lie Factor") + theme_bw() + scale_colour_discrete("Function Type") +
  theme(legend.position="none") + geom_line(data=overall.mean, aes(x=mean, y=f, group=functions), colour="black") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
# ggsave("figure/fig-spaghettiIndivDistsW.pdf", width=6, height=6, units="in")
# 
# posterior.modes <- ddply(test, .(ip.id, test_param), summarise, theta=mean[which.max(f)])
# qplot(data=posterior.modes, x=theta, geom="density",colour=test_param, fill=test_param, alpha=I(.25)) + ylab("Density") + xlab("Individual Posterior Lie Factor Mode")


#' User's estimates for Sine vs other experiment.
byuser <- dcast(test.post.indiv[, c(1, 2, 4)], ip.id~test_param)
byuser <- melt(byuser, id.vars=c("ip.id", "Sin"), variable.name="Experiment2", value.name="w")
qplot(data=byuser, x=Sin, y=w, colour=Experiment2, geom="point")  + 
  scale_colour_manual(values=c("red", "blue")) + 
  geom_smooth(method="lm") + theme_bw() + 
  ylab("Weight in Exp 2") + xlab("Weight for Sine")

setwd("../")
save.image("./code/BayesAnalysisW.Rdata")
setwd("./code/")