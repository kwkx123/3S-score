#utf-8 encoding
setwd("D:\\胶质瘤\\GBMmanuscript\\revision1\\3S-score-main\\3S-score-main\\genetic_algorithm")###please change this path to the "genetic_algorithm" file
library(MASS)
library(genalg)
library(GA)
library(timeROC)
library(survivalsvm)
library(survminer)
library(survival)
library(survivalROC)
library("gbm")
library(randomForestSRC)
library(ipflasso)
library(party)
data=read.table(".\\693\\693riskfile.csv",sep = ",",row.names = 1,header = T)
###
gatest=read.table(".\\TCGA\\TCGAriskfile.csv",sep = ",",row.names = 1,header = T)
gaevaluate <- function(indices) {
  result = 1
  if (sum(indices) > 0) {
    newdata=cbind(data[,1:2],data[,(which(indices[]==1)+2)])
    newgatest=cbind(gatest[,1:2],gatest[,(which(indices[]==1)+2)])
    colnames(newgatest)=colnames(newdata)
    
    rsfmodel <- rfsrc(Surv(futime, fustat) ~ ., newgatest, ntree = 50, 
                      nodedepth = 10,importance = TRUE)
    
    pred = predict(object = rsfmodel, 
                   newdata = newgatest)$predicted
    predicresult=Hmisc::rcorr.cens(-pred, Surv(newgatest$futime, newgatest$fustat))
    
    result=1-predicresult[1]
    
  }
  result
}
######30-time GA
bestout=data.frame()
averageout=data.frame()
for (i in 1:30) {
  print(i)
  monitor <- function(obj) {
    minEval = min(obj$evaluations);
  }
  woppa <- rbga.bin(size=47, mutationChance=0.21, 
                    evalFunc=gaevaluate, verbose=TRUE, monitorFunc=monitor,popSize=50, iters=50)
  bestout=rbind(bestout,woppa$best)
  averageout=rbind(averageout,woppa$mean)
}
fix(averageout)

write.table(bestout,
            file="bestout.csv",
            sep=",",
            quote = F,
            row.names=T)


write.table(averageout,
            file="averageout.csv",
            sep=",",
            quote = F,
            row.names=T)

####to draw curve of fitness value
bestout=read.table("averageout.csv",sep = ",",header = T)
fix(bestout)
colnames(bestout)[1]="C1"
colnames(bestout)[41]="C41"
c1name=rep("C1",30)
c41name=rep("C41",30)
newdata=cbind(c(c1name,c41name),c(bestout[,c("C1")],bestout[,c("C41")]))
fix(newdata)
trait="group"
bioCol=c("#0066FF","#FF0000")


rt2=newdata
colnames(rt2)=c("trait", "Score")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
rt2[,1]=as.character(rt2[,1])
rt2[,2]=as.numeric(rt2[,2])
fix(rt2)
data.class(rt2)
rt2=as.data.frame(rt2)
rt2[,1]=as.character(rt2[,1])
rt2[,2]=as.numeric(rt2[,2])
boxplot=ggboxplot(rt2, x="trait", y="Score", fill="trait",
                  xlab=trait,
                  ylab="Fitness value",
                  legend.title=trait,
                  palette=bioCol
)+  geom_jitter(width =0.2)+ 
  stat_compare_means(comparisons=my_comparisons,method= "t.test")
pdf(file=".\\50_average_boxplot.pdf",width=4,height=4.5)
print(boxplot)
dev.off()



###plot
data=read.table(".\\forbarplot.csv",sep = ",",header = T)
fix(data)
library(reshape2)
middledata=melt(data,id.vars = c("group","ed"))


library(foreign)
library(ggplot2)
library(tidyverse) 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
datass<- summarySE(middledata, measurevar="value", groupvars=c("group","ed"))
fix(datass)


pdf(file="50turn_barplot.pdf", width=8, height=5)
ggplot(datass, aes(x=ed, y=value, colour=group)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), colour="black", width=.3) +
  geom_line() +
  geom_point() + # 21 is filled circle
  xlab("Rounds") +
  ylab("Fitness value") +
  ggtitle("30-turn genetic algorithm result") +
  theme_bw()

dev.off()



#####generation=80 results###############################################################
bestout=data.frame()
averageout=data.frame()
for (i in 1:30) {
  print(i)
  monitor <- function(obj) {
    minEval = min(obj$evaluations);
    plot(obj, type="hist");
  }
  woppa <- rbga.bin(size=47,zeroToOneRatio=30,
                    evalFunc=gaevaluate, verbose=TRUE, monitorFunc=monitor,popSize=50, iters=80)
  bestout=rbind(bestout,woppa$best)
  averageout=rbind(averageout,woppa$mean)
}


write.table(bestout,
            file="80bestout.csv",
            sep=",",
            quote = F,
            row.names=T)


write.table(averageout,
            file="80averageout.csv",
            sep=",",
            quote = F,
            row.names=T)

data=read.table(".\\80forbarplot.csv",sep = ",",header = T)
library(reshape2)
middledata=melt(data,id.vars = c("group","ed"))
library(foreign)
library(ggplot2)
library(tidyverse) 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
datass<- summarySE(middledata, measurevar="value", groupvars=c("group","ed"))


pdf(file="80turn_barplot.pdf", width=8, height=5)
ggplot(datass, aes(x=ed, y=value, colour=group)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), colour="black", width=.3) +
  geom_line() +
  geom_point() + # 21 is filled circle
  xlab("Rounds") +
  ylab("Fitness value") +
  ggtitle("30-turn genetic algorithm result") +
  theme_bw()

dev.off()


####boxplot
bestout=read.table("80averageout.csv",sep = ",",header = T)
fix(bestout)
colnames(bestout)[1]="C1"
colnames(bestout)[57]="C57"
c1name=rep("C1",30)
c41name=rep("C57",30)
newdata=cbind(c(c1name,c41name),c(bestout[,c("C1")],bestout[,c("C57")]))
fix(newdata)
trait="group"
bioCol=c("#0066FF","#FF0000")


rt2=newdata
colnames(rt2)=c("trait", "Score")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
rt2[,1]=as.character(rt2[,1])
rt2[,2]=as.numeric(rt2[,2])
rt2=as.data.frame(rt2)
rt2[,1]=as.character(rt2[,1])
rt2[,2]=as.numeric(rt2[,2])
boxplot=ggboxplot(rt2, x="trait", y="Score", fill="trait",
                  xlab=trait,
                  ylab="Fitness value",
                  legend.title=trait,
                  palette=bioCol
)+  geom_jitter(width =0.2)+expand_limits(y=0.085) +
  stat_compare_means(comparisons=my_comparisons,method= "t.test")
pdf(file=".\\80_average_boxplot.pdf",width=4,height=4.5)
print(boxplot)
dev.off()


