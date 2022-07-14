#utf-8 encoding
setwd("D:\\胶质瘤\\上传的代码\\nomogram\\nomogram")#change the path to "nomogram" file. 
library(survival)
library(survminer)
library(timeROC)
library(rms)
library(regplot)
#firstly we want to establish the nomogram in training set(CGGA325).
riskFile="trainRisk.txt"     #
cliFile="clinical.txt"       #

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]

samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
regplot(res.cox,
              plots = c("density", "boxes"),
              clickable=F,
              title="",
              points=T,
              droplines=TRUE,
              #observation=rt[2,],
              rank="sd",
              failtime = c(5),
              prfail = T,
             showP = F)


save(res.cox,file=".\\nomogrammodel.RData") #save the nomogram
#draw the ROC curve in the training set

nomoRisk=predict(res.cox,newdata = rt, type="risk")
rt$nomoRisk=nomoRisk
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
	           marker=rt$nomoRisk, cause=1,
	           weighting='aalen',
	           times=c(1,3,5), ROC=TRUE)
pdf(file="TrainingROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
	    c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	      paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	      paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))
	      ),
	      col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()
#############draw the nomogram in the CGGA301(test set)
load("nomogrammodel.RData")
riskFile="301riskTest.txt"     #
cliFile="test301clinical.txt"       #

#
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)

#
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

nomoRisk=predict(res.cox,newdata = rt, type="risk")
rt$nomoRisk=nomoRisk
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$nomoRisk, cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)
pdf(file="301riskTestROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))
       ),
       col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()





#############draw the nomogram in the GSE4271(test set)
load("nomogrammodel.RData")
riskFile="4271riskTest.txt"     #
cliFile="test4271clinical.txt"       #

#
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)

#
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

nomoRisk=predict(res.cox,newdata = rt, type="risk")
rt$nomoRisk=nomoRisk
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$nomoRisk, cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)
pdf(file="4271riskTestROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))
       ),
       col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()


