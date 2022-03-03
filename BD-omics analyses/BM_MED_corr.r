#data1 = read.table("cluster.profile",sep = "\t",check.names = F,stringsAsFactors = F,
#                   header = T,row.names = 1)
data1 = read.csv("../../数据/BM.txt",sep = "\t",check.names = F,
                 stringsAsFactors = F,header = T,row.names = 1)
mapping = read.table("../mapping_file",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T)
row.names(mapping) = mapping[,1]
mapping = mapping[colnames(data1),]
data1 = data1[,mapping[,1]]
data1 = data1[apply(data1,1,sum)!=0,]
filter_data = read.table("filter_variant.txt",sep = "\t",check.names = F,
                         stringsAsFactors = F,header = T,row.names = 1)
#filter_BD = filter_data[filter_data$enriched=="HC",]
data1 = data1[row.names(filter_data),]
mapping = read.table("comm.txt",sep ="\t",check.names = F,stringsAsFactors = F,header = T
                     ,row.names = 1)
data1 = data.frame(t(data1),check.names = F,stringsAsFactors = F)
data1 = data1[row.names(data1)%in%row.names(mapping),]
library(gplots)
mapping = mapping[row.names(data1),]
corr_1 = function(a,b){
  cor_r = matrix(rep(0,ncol(a)*ncol(b)),ncol(a),ncol(b))
  row.names(cor_r) = colnames(a)
  colnames(cor_r) = colnames(b)
  cor_p = matrix(rep(0,ncol(a)*ncol(b)),ncol(a),ncol(b))
  row.names(cor_p) = colnames(a)
  colnames(cor_p) = colnames(b)
  for(i in 1:ncol(a)){
    for(j in 1:ncol(b)){
      cor_1 = cor.test(as.double(a[,i]),as.double(b[,j]),method = "spearman")
      cor_r[i,j] = cor_1$estimate
      cor_p[i,j] = cor_1$p.value
    }
    print(i)
  }
  k = list(p = cor_p,r = cor_r)
  return(k) 
}
mapping = mapping[,c(3,5,26,27,28,29,30,31,32,33,34,35,36:46)]
mapping = mapping[,c(-1,-2)]
data1 = data1[,-13]
cor_1 = corr_1(data1,mapping)
cor_1_r = cor_1$r
cor_1_p = cor_1$p

heat_map <- function(x,y){
  x = x[,apply(y,2,min)<=0.05]
  y = y[,apply(y,2,min)<=0.05]
  
  x = x[apply(y,1,min)<=0.05,]
  y = y[apply(y,1,min)<=0.05,]
  y[y<0.001] = "**"
  y[y>0.001&y<0.01] = "*"
  y[y>0.01&y<0.05] = "+"
  y[y>0.05] = ""
  
  p<-heatmap.2(x,col = colorRampPalette(c("#EF3F3A", "white", "#05ADAD"))(20), 
               #split = mtcars$cyl,
               key = TRUE, symkey = FALSE, density.info = "none", 
               trace = "none", cexRow = 0.5,
               main = "Heatmap",cellnote = y,notecol = "black"#,行不聚类,Rowv = F列不聚类,Colv = FALSE,
  )
  return(p)
}
heat_map(cor_1_r,cor_1_p)
