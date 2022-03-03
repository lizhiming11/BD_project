BM = read.table("BM.txt",sep= "\t",check.names = F,
                stringsAsFactors = F,header = T,row.names = 1)
SGB = read.table("BD2_SGB.profile",sep ='\t',check.names = F,
                 stringsAsFactors = F,header = T,row.names = 1)
mapping = read.table("mapping_file",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T,row.names = 1)
BM = BM[,colnames(BM)%in%row.names(mapping)]
mapping = mapping[colnames(BM),]

SGB = SGB[,colnames(BM)]
SGB = SGB[apply(SGB, 1, sum)!=0,]

library('made4')
library('ggplot2')
library('DMwR')
df2 = t(BM)
df3 = t(SGB)
dudi.chem <- dudi.pca(df2, scale=TRUE, scan=FALSE, nf=3)
dudi.topo <- dudi.pca(df3, scale=TRUE, scan=FALSE, nf=2)
dudi.chem$eig/sum(dudi.chem$eig)
all.equal(dudi.chem$lw,dudi.topo$lw)
coia.chem.topo <- coinertia(dudi.chem,dudi.topo, scan=FALSE, nf=2)
coia.chem.topo$eig[1]/sum(coia.chem.topo$eig) 
summary(coia.chem.topo)
randtest(coia.chem.topo, nrepet=99)  

eig = coia.chem.topo$eig
keqimx = coia.chem.topo$mX
colnames(keqimx)<-sapply(colnames(keqimx), FUN=function(x) paste("mx",x,sep = "_"))
keqimy = coia.chem.topo$mY
colnames(keqimy)<-sapply(colnames(keqimy), FUN=function(x) paste("my",x,sep = "_"))
keqimxy = cbind(keqimx,keqimy)
keqimxy$group = sapply(rownames(keqimxy),FUN = function(x) substr(x,1,2))
keqimxy$name = row.names(keqimxy)
#keqimxy$group = factor(keqimxy$group,levels = c("CO","KD"))
ggplot(keqimxy)+geom_point(aes(x=mx_NorS1,y=mx_NorS2,color = "#5b9bd5"),size=3)+
  geom_point(aes(x=my_NorS1,y=my_NorS2,color ="#ed7d31"),size=3,shape=15)+
  #geom_text(aes(x=my_NorS1,y=my_NorS2,label=name))+
  geom_segment(aes (x = mx_NorS1, y = mx_NorS2, xend = my_NorS1,yend = my_NorS2,alpha=1),color = "#ed7d31")+
  theme(panel.background = element_rect(fill = "transparent", color = "gray"),axis.text = element_text(color = "red"))+
  scale_color_manual(values=c("#00bfc4","#f8766d" ))+
  xlab(paste("Axis1","(",round(eig[1]/sum(eig)*100,2),"%",")",sep = ""))+
  ylab(paste("Axis2","(",round(eig[2]/sum(eig)*100,2),"%",")",sep = ""))

