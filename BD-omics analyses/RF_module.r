data1 = read.table('../肠道微生物个代谢组的关联/KO_number.txt',sep = "\t",check.names = F,
                   stringsAsFactors = F,header = T,row.names = 1)
data2 = data1[grep("K01590",row.names(data1)),]
data2 = data2[,apply(data2, 2, sum)!=0]
data2 = data.frame(t(data2),check.names = F,stringsAsFactors = F)

write.table(data2,"Histamine_SGB.txt",sep = "\t",quote = F)

#data2 = read.table("../肠道微生物个代谢组的关联/M00125_SGB.txt",sep = "\t",check.names = F,
#                   stringsAsFactors = F)
SGB_ann = read.table('../肠道微生物个代谢组的关联/SGB_ann.txt',sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T,row.names = 1)
SGB_ann = SGB_ann[row.names(data2),]
data2$ann = SGB_ann
mapping = read.table("../mapping_file",sep = "\t",check.names = F,stringsAsFactors = F,
                     header = T,row.names = 1)
BM = read.table("../代谢组/数据/BM.txt",sep ="\t",check.names = F,stringsAsFactors = F,
                header = T,row.names = 1)
mapping = mapping[colnames(BM),]
SGB_prof = read.table("../BD2_SGB.profile",sep= "\t",check.names = F,
                      stringsAsFactors = F,header = T,row.names = 1)
SGB_prof = SGB_prof[,row.names(mapping)]
SGB_prof = SGB_prof[data2[,1],]
a1_name = "Riboflavin"
BM = BM[a1_name,]
BM = data.frame(t(BM),check.names = F,stringsAsFactors = F)
SGB_prof = data.frame(t(SGB_prof),check.names = F,stringsAsFactors = F)
library(randomForest)
Species = SGB_prof
BM = BM[,1]
fit <- randomForest(BM~ ., data=Species, importance=TRUE, proximity=TRUE, ntree=1000)
imp <- importance(fit)
impvar <- imp[order(imp[,1],decreasing = TRUE),]
Species = Species[,row.names(impvar)]
RMSE = matrix(rep(0,nrow(Species)*50),nrow = nrow(Species),ncol = 50)
RMSE = data.frame(RMSE)
colnames(RMSE)[1:2] = colnames(Species)[1:2]
#colnames(RMSE)[1] = "GUT_GENOME281138,GUT_GENOME143211"
row.names(RMSE) = row.names(Species)
result = matrix(rep(0,250),nrow = 50,ncol = 5)
colnames(result) = c("RMSE","DRMSE",'COR','Q2',"name")
result[1:2,5] = colnames(Species)[1:2]
#result[1,5] = "GUT_GENOME281138,GUT_GENOME143211"
fun_q <- function(x){
  rmse = sqrt(mean((as.double(BM)-as.double(x))^2))
  drmse = rmse/(max(as.double(BM))-min(as.double(x)))
  Q2 = 1-sum((BM-as.double(x))**2)/sum((BM-mean(as.double(x)))**2)
  a = cor(as.double(BM),as.double(x))
  b = rmse 
  m = drmse
  d = Q2
  k = c(b,m,a,d)
  return(k)
}
random_test = function(x){
  fit <- randomForest(BM~ ., data=x, importance=TRUE, 
                      proximity=TRUE, ntree=1000)
  a <- as.double(fit$predicted)
  return(a)
}
leave_one = function(x){
  a = c()
  for(i in 1:nrow(x)){
    traindata = x[-i,]
    train_BM = BM[-i]
    testdata = x[i,]
    test_BM = BM[-i]
    fit <- randomForest(train_BM~ ., data=traindata, importance=TRUE, 
                        proximity=TRUE, ntree=1000)
    preds <- predict(fit, testdata)
    m <- as.double(preds)
    a = c(a,m)
  }
  return(a)
}
calculate_q = function(x){
  name1 = colnames(Species2)[x]
  Species3 = cbind(Species1[,-ncol(Species1)],Species2[,x])
  colnames(Species3)[ncol(Species3)] = name1
  colnames(Species3)[1] = name_one
  n = random_test(Species3)
  Q_value = fun_q(n)[4]
  return(c(name1,Q_value))
}
set.seed(123)
Species1 = Species[,c(1,2)]
for(j in 1:(ncol(Species)-1)){
  m = random_test(Species1)
  a = fun_q(m)[4]
  Species2 = Species[,!colnames(Species)%in%colnames(Species1)]
  name_one = colnames(Species1)[1]
  #k = 1:ncol(Species2)
  #cl <- makeCluster(2) #CPU
  #results <- parLapply(cl,k,calculate_q) 
  #res.df <- do.call('rbind',results) 
  #stopCluster(cl) 
  for(i in 1:ncol(Species2)){
    name1 = colnames(Species2)[i]
    Species3 = cbind(Species1[,-ncol(Species1)],Species2[,i])
    colnames(Species3)[ncol(Species3)] = name1
    colnames(Species3)[1] = name_one
    n = random_test(Species3)
    Q_value = fun_q(n)[4]
    if(Q_value>a){
      Species1 = Species3
      a = Q_value
      m = n
    }
  }
  RMSE[,(j+1)] =m
  name1 = colnames(Species1)[ncol(Species1)]
  colnames(RMSE)[j+1] = name1
  result[(j+1),5] = name1
  result[(j+1),1:4] = fun_q(RMSE[,(j+1)])
  Species4 = Species[,!colnames(Species)%in%colnames(Species1)]
  Species1 = cbind(Species1,Species4[,1])
  name1 = colnames(Species4)[1]
  colnames(Species1)[ncol(Species1)] = name1
  if(result[50,5]!=0){
    break
  }
  print(j)
}

write.table(result,paste(a1_name,"result.txt",sep = "_"),quote = F,sep = "\t",row.names = F)
write.table(RMSE,paste(a1_name,"RMSE.txt",sep = "_"),quote = F,sep = "\t",row.names = F)

result = data.frame(result,stringsAsFactors = F)
result$name = as.character(result$name)

b = result[which.max(result[,4]),]
name3 = c(unlist(strsplit(colnames(RMSE)[1],",")),result[2:which.max(result[,4]),5])
data_new = Species[,name3]
random_value = leave_one(data_new)
random_q = fun_q(random_value)

name = b[,5]
row.names(b) = a1
b[1,5] = paste(colnames(data_new),collapse =";")
b[1,1:4] = random_q
write.table(b,paste(a1_name,"result_filter.txt",sep = "_"),quote = F,sep = "\t")

fit1 <- randomForest(BM~ ., data=data_new, importance=TRUE, 
                     proximity=TRUE, ntree=1000)
imp <- importance(fit1)
impvar <- imp[order(imp[,1],decreasing = TRUE),]
fit1$name = row.names(impvar)
write.table(impvar,paste(a1_name,"import.txt",sep = "_"),quote = F,sep = "\t")
save(fit1,file = paste(a1_name,"_module.RData",sep = ""))

