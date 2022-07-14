###below is the nomogram construction process.
library(survival)
library(survminer)
library(timeROC)
library(survivalROC)
library(rms)
library(regplot)
riskFile="trainRisk.txt"     #input risk file (training set)
cliFile="clinical.txt"       #clinical data
setwd("D:\\胶质瘤\\上传的代码\\3Sscore\\permutationtest")     #change to the "permutationtest" file

#input risk file (training set)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "risk")]

#input clinical data
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]

#combine two tables
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)##have a look at the table.
#draw our nomogram using below code.the model is "res.cox".
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
regplot(res.cox,
        plots = c("density", "boxes"),
        clickable=F,
        title="",
        points=T,
        droplines=TRUE,
        #observation=rt[2,],
        rank="sd",
        failtime = c(1,3,5),
        prfail = T,
        showP = F)


nomoRisk=predict(res.cox,newdata = rt, type="risk")
rt$nomoRisk=nomoRisk
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$nomoRisk, cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)
pdf(file="trainingROC.pdf", width=5, height=5)
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




#####then we use the model "res.cox" to predict test set and draw the ROC curve.
riskFile="301riskTest.txt"     #input risk file (test set)
cliFile="test301clinical.txt"       #input clinical data (test set)

#input risk file (test set)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "risk")]

#input clinical data (test set)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
#cli$age=as.numeric(cli$age)

#combine two tables
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
fix(rt)#have a look at the test set



nomoRisk=predict(res.cox,newdata = rt, type="risk")
###nomoRisk is the final score, and we can use it to draw our ROC curve by below code.

rt$nomoRisk=nomoRisk
ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt$nomoRisk, cause=1,
               weighting='aalen',
               times=c(1,2,3), ROC=TRUE)
pdf(file="testROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))
       ),
       col=c("green","blue",'red'),lwd=2,bty = 'n')
dev.off()



##below is the permutation test process of ROC curve
permutation_test=function(survivetime,survivestate,markervalue)
{
        roc1=survivalROC(Stime=survivetime, survivestate, marker = markervalue, predict.time =1, method="KM")
        
        if(roc1$AUC>=0.5)
        {
                #roc1=survivalROC(Stime=survivetime, status=survivestate, marker = markervalue, predict.time =1, method="KM")
                result=c()
                for(i in 1:300)
                {
                        roc=survivalROC(Stime=survivetime, status=survivestate,marker = sample(markervalue,length(markervalue)), predict.time =1, method="KM")
                        result=c(result,roc$AUC)
                        print(paste0("first_year",i,"turn"))
                }
                
                p1=length(which(result[]>roc1$AUC))/300
        }
        if(roc1$AUC<0.5)
        {
                roc1=survivalROC(Stime=survivetime, status=survivestate, marker = -markervalue, predict.time =1, method="KM")
                result=c()
                for(i in 1:300)
                {
                        roc=survivalROC(Stime=survivetime, status=survivestate,marker = sample(-markervalue,length(-markervalue)), predict.time =1, method="KM")
                        result=c(result,roc$AUC)
                        print(paste0("first_year",i,"turn"))
                }
                
                p1=length(which(result[]>roc1$AUC))/300
        }
        
        
        
        
        
        roc3=survivalROC(Stime=survivetime, survivestate, marker = markervalue, predict.time =2, method="KM")
        if(roc3$AUC>=0.5)
        {
                #roc3=survivalROC(Stime=survivetime, status=survivestate, marker = markervalue, predict.time =3, method="KM")
                result=c()
                for(i in 1:300)
                {
                        roc=survivalROC(Stime=survivetime, status=survivestate,marker = sample(markervalue,length(markervalue)), predict.time =2, method="KM")
                        result=c(result,roc$AUC)
                        print(paste0("secound_year",i,"turn"))
                }
                
                p3=length(which(result[]>roc3$AUC))/300
        }
        if(roc3$AUC<0.5)
        {
                roc3=survivalROC(Stime=survivetime, status=survivestate, marker = -markervalue, predict.time =2, method="KM")
                result=c()
                for(i in 1:300)
                {
                        roc=survivalROC(Stime=survivetime, status=survivestate,marker = sample(-markervalue,length(-markervalue)), predict.time =2, method="KM")
                        result=c(result,roc$AUC)
                        print(paste0("secound_year",i,"turn"))
                }
                
                p3=length(which(result[]>roc3$AUC))/300
        }
        
        
        
        
        
        
        roc5=survivalROC(Stime=survivetime, survivestate, marker = markervalue, predict.time =3, method="KM")
        
        if(roc5$AUC>=0.5)
        {
                #roc3=survivalROC(Stime=survivetime, status=survivestate, marker = markervalue, predict.time =3, method="KM")
                result=c()
                for(i in 1:300)
                {
                        roc=survivalROC(Stime=survivetime, status=survivestate,marker = sample(markervalue,length(markervalue)), predict.time =3, method="KM")
                        result=c(result,roc$AUC)
                        print(paste0("third_year",i,"turn"))
                }
                
                p5=length(which(result[]>roc5$AUC))/300
        }
        if(roc5$AUC<0.5)
        {
                roc5=survivalROC(Stime=survivetime, status=survivestate, marker = -markervalue, predict.time =3, method="KM")
                result=c()
                for(i in 1:300)
                {
                        roc=survivalROC(Stime=survivetime, status=survivestate,marker = sample(-markervalue,length(-markervalue)), predict.time =3, method="KM")
                        result=c(result,roc$AUC)
                        print(paste0("third_year",i,"turn"))
                }
                
                p5=length(which(result[]>roc5$AUC))/300
        }
        
        
        c(p1,p3,p5)
        
        
}


permutation_test(rt$futime,rt$fustat,rt$nomoRisk)###finally we will get 3 p value for the 1st, 2nd and 3rd years, and I added the * label in Acrobat Illustrator software.



