library(vegan)
library(ggplot2)


species = read.csv("./BD2_SGB.profile",sep = "\t",check.names = F,
                   stringsAsFactors = F,header = T,row.names = 1)
#species = species[apply(species, 1, sum)!=0]
mapping = read.csv("./comm.txt",sep = "\t",check.names = F,
                   stringsAsFactors = F,header = T,row.names = 1)
mapping = mapping[row.names(mapping)%in%colnames(species),]
species = species[,row.names(mapping)]
species = species[apply(species, 1, sum)!=0,]
data1 = matrix(rep(0,ncol(mapping)*3),ncol = 3)
data1[,1] = colnames(mapping)
species = data.frame(t(species),check.names = F,stringsAsFactors = F)
for(i in 1:nrow(data1)){
  m = species[!is.na(mapping[,i]),]
  if(sum(is.na(mapping[,i]))==nrow(mapping)){next}
  n = mapping[,i][!is.na(mapping[,i])]
  if(length(unique(n))==1){next}
  a = adonis(m~n,permutations = 999,method = "bray")
  data1[i,2] = a$aov.tab[1,5]
  data1[i,3] = a$aov.tab[1,6]
}

colnames(data1) = c("Med","R","P")
data1 = data.frame(data1,check.names = F,stringsAsFactors = F)
data1$R = as.double(data1$R)
data1 = data1[data1$P<0.05,]
data1 = data1[order(data1$R,decreasing = T),]
data1$Med = factor(data1$Med,levels = rev(data1$Med))
ggplot(data1,aes(x=Med,y = R))+geom_bar(stat='identity')+coord_flip()+
  theme_classic()
