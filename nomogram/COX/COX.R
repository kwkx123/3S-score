#utf-8 encoding
setwd("D:\\胶质瘤\\上传的代码\\nomogram\\COX")     #change path to "COX" file.

library(survival)     


############forest plot function############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  
  
  
	rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrLow[hrLow<0.001]=0.001
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	pdf(file=forestFile, width=6.5, height=4.8)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
		
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	LOGindex=2 
	hrLow = log(as.numeric(hrLow),LOGindex)
	hrHigh = log(as.numeric(hrHigh),LOGindex)
	hr = log(as.numeric(hr),LOGindex)
	xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol, forestCol)
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
	a1 = axis(1,labels=F,tick=F)
	axis(1,a1,LOGindex^a1)
	dev.off()
}


  riskFile="trainRisk.txt"
  cliFile="clinical.txt"
  uniOutFile="uniCox.txt"
  multiOutFile="multiCox.txt"
  uniForest="uniForest.pdf"
  multiForest="multiForest.pdf"
  
  
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)  
	cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)  
	
	#combine the data
	sameSample=intersect(row.names(cli),row.names(risk))
	risk=risk[sameSample,]
	cli=cli[sameSample,]
	rt=cbind(futime=risk[,2], fustat=risk[,1], cli, riskScore=risk[,(ncol(risk)-1)])
	#unicox
	uniTab=data.frame()
	for(i in colnames(rt[,3:ncol(rt)])){
		 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
		 coxSummary = summary(cox)
		 uniTab=rbind(uniTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
	}
	write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")

	#multicox
	
	
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)  
	cli=read.table("clinicalmulti.txt", header=T, sep="\t", check.names=F, row.names=1)  
	
	#combine the data
	sameSample=intersect(row.names(cli),row.names(risk))
	risk=risk[sameSample,]
	cli=cli[sameSample,]
	rt=cbind(futime=risk[,2], fustat=risk[,1], cli, riskScore=risk[,(ncol(risk)-1)])
	uniTab=data.frame()
	for(i in colnames(rt[,3:ncol(rt)])){
	  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	  coxSummary = summary(cox)
	  uniTab=rbind(uniTab,
	               cbind(id=i,
	                     HR=coxSummary$conf.int[,"exp(coef)"],
	                     HR.95L=coxSummary$conf.int[,"lower .95"],
	                     HR.95H=coxSummary$conf.int[,"upper .95"],
	                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
	  )
	}
	fix(rt)
	
	uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
	multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
	multiCoxSum=summary(multiCox)
	multiTab=data.frame()
	multiTab=cbind(
	             HR=multiCoxSum$conf.int[,"exp(coef)"],
	             HR.95L=multiCoxSum$conf.int[,"lower .95"],
	             HR.95H=multiCoxSum$conf.int[,"upper .95"],
	             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	multiTab=cbind(id=row.names(multiTab),multiTab)
	write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="green")


