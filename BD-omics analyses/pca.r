data1 = read.csv("../image_profile.txt",sep = "\t",check.names = F,
                   stringsAsFactors = F,header = T,row.names = 1)
#rm_data = read.table("deleat_sample.txt",sep = "\t",check.names = F,stringsAsFactors = F)
#data1 = data1[,!colnames(data1)%in%rm_data[,1]]
mapping = read.table("mapping_file",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T)
row.names(mapping) = mapping[,1]
mapping = mapping[colnames(data1),]
data1 = data1[,mapping[,1]]

suppressPackageStartupMessages(library(ade4))
#suppressPackageStartupMessages(library(fpc))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(vegan))

incw = 1.5
inlength = 0.08
inlift = "Axis1"
inright = "Axis2"
dat = data1
mapping$name = mapping$Sample
#mapping = mapping[,c(3,1)]
data1 = mapping
data1 = data1[,c(1,2)]
colnames(data1)[2] ="Group"

data = dat
data = t(data)
#data <- t(sqrt(data))
data.dudi <- dudi.pca(data, center=TRUE, scale=F, scan=F, nf=10)
data2 <- data.dudi$li
s.class(data2, factor(mapping$GP), cpoi = 2)
incw = as.double(incw)
inlength = as.double(inlength)
classified_c = as.character(unique(data1[,2]))
#data3 = merge(data2,data1,by = "row.names")
#row.names(data3) = data3[,1]
#data3 = data3[,-1]

adonis1<-adonis(t(dat) ~ data1[,2],permutations = 999,method = "euclidean")
color = brewer.pal(8,"Set2")
phenotype <- data1[,2]
f = classified_c
Type <- factor(phenotype,levels=f)
m = data.dudi$li
n = data.dudi$c1

lift_data = m[as.character(inlift)]
row.names(lift_data) = gsub(pattern = "[.]",replacement = "-",row.names(lift_data))
right_data = m[as.character(inright)]
row.names(right_data) = gsub(pattern = "[.]",replacement = "-",row.names(right_data))
data.dudi$li = cbind(lift_data,right_data)
num1 = substring(as.character(inlift),5,6)
num2 = substring(as.character(inright),5,6)
num1_data = n[paste("CS",num1,sep = '')]
num2_data = n[paste("CS",num2,sep = '')]
data.dudi$c1 = cbind(num1_data,num2_data)

right_data_box= cbind(data1[,2],right_data)
colnames(right_data_box)[1] = "Group" 
lift_data_box = cbind(data1[,2],lift_data)
colnames(lift_data_box)[1] = "Group" 

#png(paste(infile,"PC",num1,"-","PC",num2,".png",sep = "_"), width = 1000, height = 768, res = 100)
#png("2017.Jun11.Throat.Stool.16S.otu.Tonsil.prof.Disease.PCA.png", width = 768, height = 768, res = 72)
x1 <- min(data.dudi$li[,1]) - 0.3
y1 <- min(data.dudi$li[,2]) - 0.3
x2 <- max(data.dudi$li[,1]) + 0.3
y2 <- max(data.dudi$li[,2]) + 0.3
bb <- head(data.dudi$c1[order(sqrt((data.dudi$c1[,1])^2+(data.dudi$c1[,2])^2),decreasing=T),],n=6L)[1:6,]
rownames(bb) <- gsub("^X", "", rownames(bb))
rownames(bb) <- gsub("\\S+o__", "o__", rownames(bb))
cutoff <- (x2-0.3) / abs(bb[1,1]) * inlength
d2 <- data.frame(X=bb[1:dim(bb)[1],1]*cutoff, Y=bb[1:dim(bb)[1],2]*cutoff, LAB=rownames(bb)[1:dim(bb)[1]])


d <- data.dudi$li
eig <- ( data.dudi$eig / sum( data.dudi$eig ) )

ggdata <- as.data.frame(d)
#Type = factor(Type,levels = c("Critical","Death","Severe"))
p<-ggplot(ggdata) +
  xlab("") +
  ylab("") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  #p+geom_point(aes(x=d[,1], y=d[,2], color=Type), size=2, shape=19) +
  geom_point(aes(x=d[,1], y=d[,2], color=Type), size=3) +
  #stat_ellipse(aes(x=d[,1], y=d[,2], fill=Type), size=1, geom="polygon", level=0.8, alpha=0.3) +
  geom_text_repel(data=d2, aes(x=X, y=Y, label=LAB),
            family="Helvetica", fontface="italic", size=3, check_overlap=TRUE) +
  geom_segment(data=d2, aes(x=0, y=0, xend=X, yend=Y),
               arrow = arrow(length = unit(0.3, "cm")), size=0.8, alpha=0.5)+
  #geom_segment(data=track, aes(x=x1, y=y1, xend=x2, yend=y2),
  #             arrow = arrow(length = unit(0.4, "cm")), size=1, alpha=0.8) +
  #geom_point(data=points, aes(x=x1, y=y1, color=factor(Type)),
  #           size=6, shape=19, alpha=0.7) +
  #geom_text_repel(data=ggdata, aes(x=d[,1], y=d[,2], label=row.names(ggdata)),
  #          family="Helvetica", size=3, check_overlap=TRUE) +
scale_color_manual(values=color)+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position=c(0.9,0.9),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=9
        ),
        legend.title=element_text(
          family="Helvetica",
          colour="black",
          size=9
        ),
        legend.text=element_text(
          family="Helvetica",
          colour="black",
          size=16
        )
  )+xlim(min(d[,1], 0)*incw, max(d[,1])*incw)+ylim(min(d[,2], 0)*incw, max(d[,2])*incw)
#xlim(min(d[,1], 0)*3, max(d[,1])*3)+ylim(min(d[,2], 0)*3, max(d[,2])*3)
p <- ggplotGrob(p)
right_data_box$Group = factor(right_data_box$Group, levels=f)
d <- ggplot(right_data_box)+geom_boxplot(aes(x = Group,y = right_data_box[,2],fill = Group),width = 0.5)+
  theme_bw()+theme(panel.grid =element_blank())+
  #scale_fill_wsj("colors6", "Group")
  scale_fill_manual(values=color,breaks =f)+
  guides(fill=FALSE)+theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank(),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=12
        ))+
  ylim(min(right_data_box[,2], 0)*incw, max(right_data_box[,2])*incw)+
  # ylim(min(right_data_box[,2], 0)*3, max(right_data_box[,2])*3)+
  xlab("")+ylab(paste("PC",num2," (",round(eig[as.numeric(num2)]*100,2),"%)",sep=""))
lift_data_box$Group = factor(lift_data_box$Group, levels=f)
b<- ggplot(lift_data_box)+geom_boxplot(aes(x = Group,y = lift_data_box[,2],fill = Group),width = 0.5)+
  theme_bw()+theme(panel.grid =element_blank())+coord_flip()+
  guides(fill=FALSE)+theme(axis.text.y = element_blank())+
  theme(axis.ticks.y = element_blank(),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=12
        ))+
  scale_fill_manual(values=color)+
  ylim(min(lift_data_box[,2], 0)*incw, max(lift_data_box[,2])*incw)+
  xlab("")+ylab(paste("PC",num1," (",round(eig[as.numeric(num1)]*100,2),"%)",sep=""))
a<-ggplot()+theme_bw()+theme(panel.border = element_blank(),panel.grid =element_blank(),
                             axis.text = element_blank(),axis.title = element_blank(),
                             axis.ticks = element_blank())+
  annotate("text", x=1, y=40, label=paste("P.value =",round(adonis1[[1]][6][[1]][1],4),'\n',
                                          "R2      =",round(adonis1[[1]][5][[1]][1],4)), size=3.5)

a <- ggplotGrob(a)
d <- ggplotGrob(d)
b <- ggplotGrob(b)


grid.arrange(d,p,a,b,ncol=2,widths=c(1,4),heights = c(4,1))




