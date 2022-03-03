data1 = read.table("diff_sep.txt",sep ="\t",check.names = F,
                   stringsAsFactors = F,header = T,row.names = 1)
data1 = data1[data1$q<0.05,]
data1 = data1[data1$BD_mean/data1$HC_mean>=2|data1$BD_mean/data1$HC_mean<=0.5,]
data1$families = data1$ann
data1$families = gsub(";g__.*","",data1$families)

data1$families = gsub(".*;f__","f__",data1$families)
data1$phylum = data1$ann
data1$phylum = gsub(";c__.*","",data1$phylum)

HC = data1[data1$enriched=="HC",]
BD = data1[data1$enriched=="BD",]
HC_table = data.frame(table(HC$families))
BD_table = data.frame(table(BD$families))
data2 = merge(HC_table,BD_table,by = "Var1",all = T)
colnames(data2) = c("families","HC_enriched","BD_enriched")
data2$BD_enriched = -data2$BD_enriched

data3 = data1[,c(8,9)]
data3 = data3[!duplicated(data3$families),]
row.names(data3) = data3$families
row.names(data2) = data2$families
data3 = data3[order(data3$phylum),]
data2 = data2[row.names(data3),]
data2$families = factor(data2$families,levels = data2$families)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
ggplot(data2) +
  geom_segment( aes(x=families, xend=families, y=1, yend=HC_enriched), color="grey") +
  geom_segment( aes(x=families, xend=families, y=1, yend=BD_enriched), color="grey") +
  geom_point( aes(x=families, y=HC_enriched), color=rgb(0.2,0.7,0.1,0.5), size=5 ) +
  geom_point( aes(x=families, y=BD_enriched), color=rgb(0.7,0.2,0.1,0.5), size=5 ) +
  geom_text(aes(x=families, y=BD_enriched,label=BD_enriched))+
  geom_text(aes(x=families, y=HC_enriched,label=HC_enriched))+
  coord_flip()+
  theme_classic()+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  )  +
  xlab("") +
  ylab("Value of Y")
