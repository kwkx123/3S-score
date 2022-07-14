#UTF-8 encoding
setwd("D:\\胶质瘤\\上传的代码\\3Sscore\\ROC_survival_curve")#change the path to "ROC_survival_curve" curve.
library(survminer)
library(survival)
library(survivalROC)
library(timeROC)
library(MASS)
library(genalg)
library(timeROC)
library(survivalsvm)
library(survminer)
library(survival)
library(survivalROC)
library("gbm")
library(randomForestSRC)
library(ipflasso)
library(party)

cindex=c()
pKM=c()
load(".\\3Smodel.RData")
setname=c("TCGA","693","325","301","475","4271","4412","43378","7696","74187","83300","combined_test")
for (ID in setname) {
  print(ID)
  if(ID=="combined_test")
  {
    expdata=read.table(paste0(".\\",ID,"\\",ID,"riskfile.csv"), header=T, sep=";", check.names=F,row.names = 1)
  }
  else
  {
    expdata=read.table(paste0(".\\",ID,"\\",ID,"riskfile.csv"), header=T, sep=",", check.names=F,row.names = 1)
  }
  tcgaOut=expdata
  riskScoreTest= predict(object = rsfmodel,newdata=tcgaOut)$predicted
  
  #we obtained the RGP risk sore of this test set
  riskScoreTest=riskScoreTest
  table(riskScoreTest)
  cindexresult=Hmisc::rcorr.cens(-riskScoreTest, Surv(tcgaOut$futime, tcgaOut$fustat))
  cindex=c(cindex,cindexresult[1])
  medianTrainRisk=median(riskScoreTest)
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  riskfile=cbind(id=rownames(cbind(tcgaOut,riskScoreTest,riskTest)),cbind(tcgaOut,riskScore=riskScoreTest,risk=riskTest))
  write.table(riskfile,
              file=paste0(".\\score\\",ID,"riskTest.txt"),
              sep="\t",
              quote=F,
              row.names=F)
  
  
  
  #to draw K-M curves using below codes
  rt=riskfile
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     title=ID,
                     pval=pValue,
                     pval.size=6,
                     conf.int=F,
                     legend.title="risk group",
                     legend.labs=c("high","low"),
                     font.legend=12,
                     fontsize=4,
                     xlab="Time(years)",
                     ylab="Overall survival",
                     break.time.by = 2,
                     palette=c("red","blue"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.height=.25)
  pdf(file=paste0(".\\KM\\",ID,"KM.pdf"),width=5.5,height=5)#we can get K-M curve of this test set in "survival" file.
  print(surPlot)
  dev.off()
  pKM=c(pKM,pValue)
  
  
  
  
  ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
                 marker=rt$riskScore, cause=1,
                 weighting='aalen',
                 times=c(1,3,5), ROC=TRUE)
  pdf(file=paste0(".\\ROC\\",ID,"multiROC.pdf"), width=5, height=5)
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
  
  
  print(ID)
  
}
