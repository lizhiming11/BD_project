 
library(vegan) 
library(ape)
library(ggplot2)
library(RColorBrewer)
library(ade4)
library(ggrepel)


inf="./KO_profile"
inp='mapping_file'

prof <- read.csv(file=inf,sep="\t",row.names=1,check.names = F,stringsAsFactors = F)  
map <- read.csv(file=inp,sep="\t",check.names = F,stringsAsFactors = F)
map = map[,c(1,2)]
prof = prof[,map[,1]]
prof <- t(prof)
prof <- prof[,colSums(prof)!=0]
prof <- sqrt(prof)

#dbRDA========================================================================================
ord <- capscale(prof~GP,map,distance = 'bray') 
sp = scores(ord,choices = 1:2,display = "sites")
s.class(sp,factor(map[,2]),col =brewer.pal(4,"Set1"))


p<-plot(ord)
site.scores <- as.data.frame(scores(p, "site")) 
species.scores <- as.data.frame(scores(p, "species"))
site.scores$lab <- row.names(site.scores)
site.scores <- cbind(site.scores,map)[-4]
site.scores$z <- NA 

x1 <- min(site.scores[,1]) - 0.3
y1 <- min(site.scores[,2]) - 0.3
x2 <- max(site.scores[,1]) + 0.3
y2 <- max(site.scores[,2]) + 0.3
inlength =0.2
bb <- head(species.scores[order(sqrt((species.scores[,1])^2+(species.scores[,2])^2),decreasing=T),],n=8L)[1:8,]


cutoff <- (x2-0.3) / abs(bb[1,1]) * inlength
d2 <- data.frame(X=bb[1:dim(bb)[1],1]*cutoff, Y=bb[1:dim(bb)[1],2]*cutoff, LAB=rownames(bb)[1:dim(bb)[1]])
d2$z = NA
color = c(brewer.pal(2,"Set1"))
site.scores$color = color[1]
site.scores$color[site.scores$GP=="BD"] = color[2]
#site.scores$color[site.scores$group=="Elderly"] = color[3]
#site.scores$color[site.scores$group=="Longevity"] = color[4]
ggplot() +
  scale_color_continuous(low='grey80',high = '#272EAF')+
  geom_point(data = site.scores,size=2,
             aes(x=CAP1,y=MDS1),colour = site.scores$color)+
  geom_text_repel(data=d2, aes(x=X, y=Y, label=LAB),
                 family="Helvetica", fontface="italic", size=3, check_overlap=TRUE) +
  stat_ellipse(data = site.scores,
               aes(x=CAP1,y=MDS1, fill=color), size=1, geom="polygon", level=0.8, alpha=0.3) +
  geom_segment(data=d2, aes(x=0, y=0, xend=X, yend=Y),
               arrow = arrow(length = unit(0.3, "cm")), size=0.8, alpha=0.5)+
  theme_bw()+
  coord_cartesian(ylim = c(-2, 2))+
  labs(title="", x="CAP1", y="MDS1")+
  theme(
    panel.grid = element_blank()
  )


