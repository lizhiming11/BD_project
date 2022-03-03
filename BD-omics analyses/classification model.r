data1 = read.table("../BD2_SGB.profile",sep = "\t",check.names = F,stringsAsFactors = F,
                   header = T,row.names = 1)
data1_filter = read.table("filter_SGB.txt",sep ="\t",check.names = F,stringsAsFactors = F,
                          header = T,row.names = 1)
data1 = data1[row.names(data1_filter),]
BM = read.table("mapping_file",sep = "\t",check.names = F,
                stringsAsFactors = F,header = T,row.names = 1)
BM$type = 0
BM$type[BM$GP=="BD"] = 1
data1 = data.frame(t(data1),check.names = F,stringsAsFactors = F)
library(randomForest)
Species = data1
BM = BM[,2]
fit <- randomForest(BM~ ., data=Species, importance=TRUE, proximity=TRUE, ntree=1000)
imp <- importance(fit)
impvar <- imp[order(imp[,1],decreasing = TRUE),]
Species = Species[,row.names(impvar)]
RMSE = matrix(rep(0,nrow(Species)*50),nrow = nrow(Species),ncol = 50)
RMSE = data.frame(RMSE)
colnames(RMSE)[1:2] = colnames(Species)[1:2]
row.names(RMSE) = row.names(Species)
result = matrix(rep(0,250),nrow = 50,ncol = 5)
colnames(result) = c("RMSE","DRMSE",'COR','Q2',"name")
result[1:2,5] = colnames(Species)[1:2]
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

write.table(result,"result.txt",quote = F,sep = "\t",row.names = F)
write.table(RMSE,"RMSE.txt",quote = F,sep = "\t",row.names = F)

result = data.frame(result,stringsAsFactors = F)
result$name = as.character(result$name)

b = result[which.max(result[,4]),]
name3 = c(unlist(strsplit(colnames(RMSE)[1],",")),result[2:which.max(result[,4]),5])
data_new = Species[,name3]
random_value = leave_one(data_new)
random_q = fun_q(random_value)

name = b[,5]
#row.names(b) = a1
b[1,5] = paste(colnames(data_new),collapse =";")
b[1,1:4] = random_q
write.table(b,"result_filter.txt",quote = F,sep = "\t")

fit1 <- randomForest(BM~ ., data=data_new, importance=TRUE, 
                     proximity=TRUE, ntree=1000)
imp <- importance(fit1)
impvar <- imp[order(imp[,1],decreasing = TRUE),]
data1_filter = data1_filter[row.names(imp),]
imp = data.frame(imp,check.names = F,stringsAsFactors = F)
imp$name = paste(data1_filter$`Taxonomy assignment`," (",gsub("GUT_","",row.names(data1_filter)),")",sep = "")
library(ggplot2)
imp <- imp[order(imp[,1],decreasing = TRUE),]
imp$name = factor(imp$name,levels = rev(imp$name))
ggplot(imp,aes(x= name,y = imp$`%IncMSE`))+
  geom_bar(stat = 'identity')+coord_flip()+theme_classic()

fit1$name = row.names(impvar)
write.table(impvar,"import.txt",quote = F,sep = "\t")
save(fit1,"module.RData")


library(pROC)
roc(BM,random_value)
data2 = data.frame(type = as.character(BM),value = random_value)
ggplot(data2,aes(x = type,y = value))+geom_boxplot()+theme_classic()
###################
roc1 <- roc(BM,random_value,
            percent=TRUE,
            partial.auc.correct=TRUE,
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            plot=F, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
#####################
roc1 <- roc(BM,random_value,
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            plot=TRUE, percent=roc1$percent,col=2)
sens.ci <- ci.se(roc1, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col=rgb(0,1,0,alpha=0.2))
plot(sens.ci, type="bars")
plot(roc1,col=2,add=T)
legend("bottomright",c(paste("AUC=",round(roc1$ci[2],2),"%"),
                       paste("95% CI:",round(roc1$ci[1],2),"%-",round(roc1$ci[3],2),"%")))
data1 = data.frame(name = row.names(Species2),predicted = random_value,experiment = BM)
write.table(data1,paste(a_name[1],"_",a_name[2],"Predicted.txt",sep=""),quote = F,sep = "\t",row.names = F)


