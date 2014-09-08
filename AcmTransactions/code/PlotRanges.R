require(ggplot2)
require(plyr)
library(msm)
# set.seed(70032608) # so that intervals don't change

turkdata.full <- read.csv("./data/turkdataclean.csv", stringsAsFactors=FALSE)
turkdata <- subset(turkdata, len>10)

library(doMC) 
registerDoMC(cores=12) 

turkdata <- ddply(turkdata, .(ip.id, test_param), transform, n=length(test_param))
turkdata <- subset(turkdata, n>4)


key <- read.csv("./data/pictureKey.csv", stringsAsFactors=FALSE)
key$pic_name <- tolower(key$pic_name)

tmp <- ddply(key, .(pic_name), transform, w.min = w[c(1, 1:5)], w.max=w[c(2:6,6)], d.min=round(lmax[c(1, 1:5)], 3), d.max=round(lmax[c(2:6, 6)], 3))

turkdata2 <- merge(turkdata, tmp, by.x=c("pic_name", "answer.w"), by.y=c("pic_name", "w"), all.x=TRUE)


ggplot(data=subset(turkdata2, difficulty!="test")) + 
  geom_segment(aes(x=ip.id, xend=ip.id, y=w.min, yend=w.max), alpha=.1) + 
  geom_jitter(aes(x=ip.id, y=answer.w), alpha=.1) + 
  facet_wrap(~test_param) + 
  theme_bw()

ggplot(data=subset(turkdata2, difficulty!="test")) + 
  geom_segment(aes(x=ip.id, xend=ip.id, y=d.min, yend=d.max), alpha=.1) + 
  geom_jitter(aes(x=ip.id, y=lmax), alpha=.1) + 
  facet_wrap(~test_param) + 
  theme_bw()

participants <- dcast(ddply(turkdata, .(ip.id, test_param), summarise, n=length(test_param)), ip.id~test_param, value.var="n")
ipsubset <- subset(participants, rowSums(is.na(participants))==0 & rowSums(participants[,2:4]>6, na.rm=TRUE)==3)$ip.id