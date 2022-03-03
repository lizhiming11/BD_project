data1 = read.csv("../BD2_SGB.profile",sep = "\t",check.names = F,
                 stringsAsFactors = F,header = T,row.names = 1)
mapping = read.table("mapping_file",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T)
data1 = data1[,mapping[,1]]
data1 = data1[apply(data1,1,sum)!=0,]
data2 = matrix(rep(0,nrow(data1)*5),ncol = 5)
colnames(data2) = c("species","BD_mean","HC_mean","p","roc")
data2[,1] = row.names(data1)
library(pROC)
for(i in 1:nrow(data2)){
  data2[i,2] = mean(as.double(data1[i,][mapping$GP=="BD"]))
  data2[i,3] = mean(as.double(data1[i,][mapping$GP=="HC"]))
  a = wilcox.test(as.double(data1[i,][mapping$GP=="BD"]),as.double(data1[i,][mapping$GP=="HC"]))
  data2[i,4] = a$p.value  
  ROC_v = roc(mapping$GP,as.double(data1[i,]))
  data2[i,5] = ROC_v$auc
}


library(fdrtool)
data2 = data.frame(data2,check.names = F,stringsAsFactors = F)
data2$p = as.double(data2$p)
a = fdrtool(data2$p,statistic = "pvalue")
data2$q = a$qval
#data2 = data2[data2[,6]<0.05,]

data2$enriched = "HC"
data2$BD_mean = as.double(data2$BD_mean)
data2$HC_mean = as.double(data2$HC_mean)
data2$enriched[data2$BD_mean>data2$HC_mean] = "BD"
SGB_ann = read.table("SGB_ann.txt",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T)
row.names(SGB_ann) = SGB_ann[,1]
SGB_ann = SGB_ann[data2[,1],]
data2$ann = SGB_ann$`Taxonomy lineage (NCBI)`
#write.table(data2,"diff_sep.txt",quote = F,sep ="\t",row.names = F)
