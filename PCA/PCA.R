#utf-8 encoding
setwd("D:\\胶质瘤\\GBM\\配对和sva的比较")#change the path to PCA
library(sva)
library(bladderbatch)
library("FactoMineR")
library("factoextra")
pca.plot = function(dat,col){
  
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               axes = c(1, 2),
               geom.ind = "point",
               col.ind = col ,
               addEllipses = TRUE,
               legend.title = "Groups"
  )
}
FPKM=read.table("allexp.csv",sep = ",",row.names = 1,header = F)
colnames(FPKM)=FPKM[1,]
FPKM=FPKM[-1,]
fix(FPKM)
for (i in 1:ncol(FPKM)) {
  FPKM[,i]=as.numeric(FPKM[,i])
}

group=read.table("group.csv",sep = ",",row.names = 1,header = T)
intername=intersect(row.names(group),colnames(FPKM))
length(intername)
group=group[intername,]
FPKM=FPKM[,intername]
sample_infor=group

pca.plot(FPKM,factor(sample_infor$group1))

combat_FPKM <- ComBat(dat = as.matrix(FPKM), batch = sample_infor$group1)
pca.plot(combat_FPKM,factor(sample_infor$group1))
####下面用配对的文件进行绘制
fix(FPKM)
rt=FPKM
tcgaPair=data.frame()
sampleNum=ncol(rt)
for(i in 1:(nrow(rt)-1)){
  for(j in (i+1):nrow(rt)){
    pair=ifelse(rt[i,]>rt[j,], 1, 0)
    pairRatio=sum(pair)/sampleNum
    if((pairRatio>0.3) & (pairRatio<0.7)){
      rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
      tcgaPair=rbind(tcgaPair, pair)
    }
  }
}

pca.plot(tcgaPair,factor(sample_infor$group1))

write.table(combat_FPKM,"allexp_aftersva.csv",sep = ",",col.names = T,row.names = T)










