library (xlsx) ### Saving to spreadsheet
library (data.table) ### Fast read of large files into R
library (GO.db)
library(preprocessCore)
library (WGCNA) ### -   Clustering software. Previously reported work done using v1.34
library (flashClust) ### Clustering software
library (ppcor) ### Partial Spearman correlations, for confounder analysis. Previously reported work done using v1.0
library (gplots) ### Plotting
library (cowplot) ### Plotting; to arrange several plots on the same page
library (ggplot2) ### Plotting
library (plyr) ### Data transformations

data1 = read.csv("../Êý¾Ý/BM.txt",sep = "\t",check.names = F,
                 stringsAsFactors = F,header = T,row.names = 1)
mapping = read.table("mapping_file",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T)
row.names(mapping) = mapping[,1]
mapping = mapping[colnames(data1),]
data1 = data1[,mapping[,1]]
data1 = data1[apply(data1,1,sum)!=0,]
filter_data = read.table("diff_cluser_BM.txt",sep = "\t",check.names = F,
                         stringsAsFactors = F,header = T,row.names = 1)
#filter_BD = filter_data[filter_data$enriched=="HC",]
data1 = data1[row.names(filter_data),]
dat = t(data1) 
dat_tmp = log2 (dat)

### Settings for WGCNA generally
cor_method          = "spearman" ### for association with clinical parameters
corFun_tmp          = "bicor"
cluster_method      = "average"
corOptions_list     = list (use = 'pairwise.complete.obs') 
corOptions_str      = "use = 'pairwise.complete.obs'"
BH_pval_asso_cutoff = 0.05
NetworkType         = "signed" ### Signed-network (as the PC1 and median profile does not make sense as a summary measure of a cluster with anticorrelated metabolit


### Settings for WGCNA on polar metabolite measurements
### Once these are established the steps below can be run
RsquareCut_val      = 0.90 ### usually ranges 0.80-0.95 but requires inspecting curves
mergingThresh       = 0.10 ### Maximum dissimilarity of module eigengenes (i.e. 1-correlation) for merging modules.
minModuleSize       = 2 ### minimum number of metabolites constituting a cluster
SoftPower           = 13 ### beta-value, main parameter to optimize
A = adjacency (dat_tmp, power = SoftPower, type = NetworkType, corFnc = corFun_tmp, corOptions = corOptions_str)
colnames (A) = rownames (A) = colnames (dat_tmp)
### Define dissimilarity based on topological overlap
dissTOM = TOMdist (A, TOMType = NetworkType)
colnames (dissTOM) = rownames (dissTOM) = colnames (dat_tmp)
### Hierarchical clustering
metaTree = flashClust (as.dist (dissTOM), method = cluster_method)
### Define modules by cutting branches
moduleLabels1 = cutreeDynamic (dendro = metaTree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize)
moduleLabels1 = labels2colors (moduleLabels1)
### Automatically merge highly correlated modules
merge = mergeCloseModules (dat_tmp, moduleLabels1, corFnc = corFun_tmp, corOptions = corOptions_list, cutHeight = mergingThresh)
### Determine resulting merged module colors
moduleLabels2 = merge$colors
### Establish eigengenes of the newly merged modules, used for cluster overall abundances
MEs = merge$newMEs
### Choose final module assignments
moduleColorsMeta = moduleLabels2
names (moduleColorsMeta) = colnames (dat_tmp)
MEsMeta = orderMEs (MEs)
rownames (MEsMeta) = rownames (dat_tmp)

### Determine relevant descriptive statistics of established clusters
### kIN: within-module connectivity, determined by summing connectivity with all
###      other metabolites in the given cluster.
### kME: bicor-correlation between the metabolite profile and module eigenvector; 
### both measures of intramodular hub-metabolite status.
kIN <-      vector (length = ncol (dat_tmp)); names (kIN) = colnames (dat_tmp)
kME <-      vector (length = ncol (dat_tmp)); names (kME) = colnames (dat_tmp)
modules <-  vector (length = ncol (dat_tmp)); names (modules) = colnames (dat_tmp)


for (module in names (table (moduleColorsMeta))) {   
  
  all.metabolites = names (dat_tmp)
  inModule = (moduleColorsMeta == module)
  module.metabolites = names (moduleColorsMeta [inModule])
  modules [module.metabolites] = module 
  kIN [module.metabolites] = sapply (module.metabolites, function (x) sum (A [x, module.metabolites]) - 1)
  datKME = signedKME (dat_tmp, MEsMeta, corFnc = corFun_tmp, corOptions = corOptions_str)
  rownames (datKME) = colnames (dat_tmp)
  kME [module.metabolites] = datKME [module.metabolites, paste ("kME", module, sep = "")]   
  
}
output = data.frame ("module" = modules, "kME" = kME, "kIN" = kIN )
write.table(output,"cluster_info.txt",sep = "\t",quote = F)
write.table(MEsMeta,"cluster.profile",sep = "\t",quote = F)
