data1 = read.table("BD2_contig.profile",sep = "\t",check.names = F,
                   stringsAsFactors = F,header = T,row.names = 1)
data1_ann = read.csv("ko_list",sep = "\t",check.names = F,stringsAsFactors = F,
                       header = T)
mapping = read.table("./mapping_file",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T,row.names = 1)
library(pROC)
data1 = data1[,row.names(mapping)]
data1_ann$SGB = gsub("_[0-9].*","",data1_ann$contig_name)
name = unique(data1_ann$SGB)
  
  
for(i in 1:length(name)){
  a_ko = data1_ann[data1_ann$SGB==name[i],]
  a= data1[a_ko$contig_name,]
  cp_num =as.data.frame(table(a_ko$KO))
  a$ko = a_ko$KO
  b = aggregate(.~ko,a,sum)
  write.table(cp_num,paste("./KO/",name[i],"_KO.cope",sep =""),quote = F,sep = "\t",row.names = F)
  write.table(b,paste("./KO/",name[i],"_KO.profile",sep =""),quote = F,sep = "\t",row.names = F)
  row.names(b) = b[,1]
  b = b[,-1]
  data2 = matrix(rep(0,nrow(b)*5),ncol = 5)
  colnames(data2) = c("KO","BD_mean","CON_mean","P","ROC")
  data2[,1] =row.names(b)
  for(j in 1:nrow(b)){
      data2[j,2] = mean(as.double(b[j,mapping$GP=="BD"]))
      data2[j,3] = mean(as.double(b[j,mapping$GP=="HC"]))
      test_1 = wilcox.test(as.double(b[j,mapping$GP=="BD"]),
                           as.double(b[j,mapping$GP=="HC"]))
      data2[j,4] = test_1$p.value
      ROC_v = roc(mapping$GP,as.double(b[j,]),
                  levels=c("BD", "HC"),direction= ">")
      data2[j,5] = ROC_v$auc
    }
  write.table(data2,paste("./KO/",name[i],"_KO.diff",sep =""),quote = F,sep = "\t",row.names = F)
}
