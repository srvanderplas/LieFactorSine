library(ggplot2)
library(reshape2)
library(plyr)

turkdata <- read.csv("./data/raw_data_turk8.csv", stringsAsFactors=FALSE)
turkdata$ip.id <- as.numeric(factor(turkdata$ip_address))
turkdata$pic_num <- as.numeric(factor(turkdata$pic_name))

usertime <- ddply(turkdata, .(ip.id, ip_address, age, gender), summarise, 
                  len = length(time_taken), 
                  avg_time = mean(time_taken), 
                  min_time = mean(time_taken), 
                  max_time = max(time_taken))
qplot(x=usertime$avg_time, geom="density")
qplot(data=usertime, x=max_time, y=min_time, geom="point", alpha=I(.5))
usertime$gets.paid <- usertime$len>10
usertime$include <- usertime$len>5

turkdata <- merge(turkdata, usertime[,c(1, 2, 3, 4, 5, 9, 10)], keep=TRUE)
turkdata$user.seq <- unlist(dlply(turkdata, .(ip.id), function(i) seq(1, nrow(i), by=1)))

# qplot(data=turkdata, x=user.seq, y=response_no, group=ip.id, geom="line", alpha=.1) + xlim(c(0,11))


key <- read.csv("./data/pictureKey.csv", stringsAsFactors=FALSE)
key2 <- dcast(key, pic_name~response_no, value.var="w")

turkdata$answer.w <- unlist(lapply(1:nrow(turkdata), function(i) subset(key2, pic_name==turkdata$pic_name[i])[1+turkdata$response_no[i]]))



dont.pay.reason <-  function(data){
  ids <- unique(data$id)
  min_id <- data$id==min(ids)
  onesubmit <- length(unique(data$id))==1
  if(prod(data$gets.paid)==1) 
    data$reason <- paste(c("Two submissions - ", "")[1+onesubmit], 
                         c("Already paid","Pay")[1+min_id], sep="") else 
  data$reason <- "Not enough responses"
  data$pay <- data$gets.paid*min_id
  unique(data[,c("id", "ip.id", "ip_address", "reason", "pay")])
}

turkers.pay <- ddply(turkdata, .(ip_address, ip.id, gets.paid), dont.pay.reason)

