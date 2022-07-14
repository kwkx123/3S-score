#utf-8 encoding
setwd("D:\\胶质瘤\\上传的代码\\3Sscore\\pan-cancer")#change this path to the "pan-cancer" file

#引用包
library(survival)
library(survminer)
library(forestplot)

inputFile="expTime.txt" 
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    
fix(rt)
rt$futime=rt$futime/365
gene=colnames(rt)[3]

#对肿瘤类型进行循环
outTab=data.frame()
for(i in levels(factor(rt[,"CancerType"]))){
	#获取单个肿瘤的数据
	rt1=rt[(rt[,"CancerType"]==i),]
	#COX分析
	cox=coxph(Surv(futime, fustat) ~ rt1[,gene], data = rt1)
	coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	outTab=rbind(outTab,
	             cbind(cancer=i,
	                   HR=coxSummary$conf.int[,"exp(coef)"],
	                   HR.95L=coxSummary$conf.int[,"lower .95"],
	                   HR.95H=coxSummary$conf.int[,"upper .95"],
			           pvalue=coxP) )
	
	#KM分析
	group=ifelse(rt1[,gene]>median(rt1[,"riskScore"]), "high", "low")
	diff=survdiff(Surv(futime, fustat) ~ group, data=rt1)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.05){
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		
		#绘制生存曲线
		fit=survfit(Surv(futime, fustat) ~ group, data = rt1)
		surPlot=ggsurvplot(fit, 
				    data=rt1,
				    title=paste0("Cancer: ",i),
				    pval=pValue,
				    pval.size=6,
				    conf.int=F,
				    legend.title=paste0(gene," levels"),
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
		#输出图形
		pdf(file=paste0("survival.",i,".pdf"), width=6, height=5, onefile=FALSE)
		print(surPlot)
		dev.off()
	}
	print(i)
}
#输出COX分析结果
write.table(outTab, file="cox.result.txt", sep="\t", row.names=F, quote=F)


############forest figure############
    #读取输入文件
  coxFile="cox.result.txt"
  forestFile="forest.pdf"
  forestCol="red"
	rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
	data=as.matrix(rt)
	fix(data)
	HR=data[,1:3]
	hr=sprintf("%.3f",HR[,"HR"])
	hrLow=sprintf("%.3f",HR[,"HR.95L"])
	hrHigh=sprintf("%.3f",HR[,"HR.95H"])
	pVal=data[,"pvalue"]
	pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
	#定义颜色
	clrs=fpColors(box=forestCol, line="darkblue", summary="royalblue")
	#定义图片文字
	tabletext <- 
		list(c(NA, rownames(HR)),
		    append("pvalue", pVal),
		    append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )
	#绘制森林图
	pdf(file=forestFile, width=9, height=6, onefile=FALSE)
	forestplot(tabletext, 
	           rbind(rep(NA, 3), HR),
	           col=clrs,
	           graphwidth=unit(50, "mm"),
	           xlog=T,
	           lwd.ci=4,
	           boxsize=0.6,
	           title="Overall survival",
	           xlab="Hazard ratio",
	           txt_gp=fpTxtGp(ticks=gpar(cex=1.1), xlab=gpar(cex = 1.25))
	           )
	dev.off()


################ROC curves
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    
rt$futime=rt$futime/365
fix(rt)
library(survivalROC)
library(survminer)
library(survival)
library(timeROC)
for(i in names(table(rt[,5])))
{
  print(i)
  cancername=i
  data=rt[which(rt[,5]==cancername),]
  
  
  ROC_rt=timeROC(T=data$futime, delta=data$fustat,
                 marker=data$riskScore, cause=1,
                 weighting='aalen',
                 times=c(1,3,5), ROC=TRUE)
  pdf(file=paste0(i,"multiROC.pdf"), width=5, height=5)
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('1st year AUC=',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('3rd year AUC=',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('5th year AUC=',sprintf("%.03f",ROC_rt$AUC[3]))
         ),
         col=c("green","blue",'red'),lwd=2,bty = 'n')
  dev.off()
  
  
}

