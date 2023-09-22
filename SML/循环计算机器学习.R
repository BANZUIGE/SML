# 循环机器学习
# 1. 包的载入 -----------------------------------------------------------------
library(randomForest)
library(caret)
library(tidyverse)
library(irr)
library(e1071)
library(ggplot2)
library(Metrics)

# 2. 数据载入 -----------------------------------------------------------------
asv <- read.csv("otu_CSS_table_84sample_combine9396.csv",header = TRUE, row.names = 1,sep = ",",check.names = F)  %>% as.data.frame()
IBI <- read.csv("IBI.csv",header = TRUE, row.names = 1,sep = ",",check.names = F)  %>% as.data.frame()

core <- read.csv("importance_otu_100seed.csv",header = TRUE,sep = ",",check.names = F)  %>% as.data.frame()
a <- core$Core[1:21]

asv_ibi <- cbind(t(asv),IBI) %>% as.data.frame() 
tail(names(asv_ibi))

# 3. 划分训练、测试集 -------------------------------------------------------------
dim(asv_ibi)  #查看数据多少行，多少列
set.seed(2)
trainlist <- createDataPartition(asv_ibi$IBI,p=2/3,list=FALSE) #将数据集划分为2部分
trainset <- asv_ibi[trainlist ,]
testset <- asv_ibi[-trainlist ,]


# 4.循环计算变量筛选后RF结果 ---------------------------------------------------------

rmse <- c()
explain <- c()
explain.select <- c()

R2 <- c()
R2.select <-  c()

kappa <-  c()
kappa.select <-  c()

seed <-  c()
seed.select <-  c()

kappa4_4 <- c()
kappa4_5 <- c()
kappa4.1_4 <- c()
kappa4.1_5 <- c()
kappa4.3_4 <- c()
kappa4.3_5 <- c()




for (i in 1:10){
  
  # 4. 由训练集建立模型 
  set.seed(i)
  rf.train <-  readRDS(paste("RDS_modle/rf.train.seed_",i,".rda",sep = ""))
#   rf.train <- randomForest(IBI ~.,   #as.factor做分类分析，否则做回归分析
#                            data = trainset,importance = TRUE)
#   explain <- c(explain, mean(unlist(rf.train[5])))
#   
# # 5. 根据建立的模型查看模型性能
#   set.seed(123)
#   # 使用测试集，评估预测性能， 或
#   predict.set <- testset
#   predict <- predict(rf.train, predict.set) %>% as.data.frame()
#   names(predict) = "predict"
#   test <- cbind(subset(predict.set,select = IBI), subset(predict,select = c("predict"))) 
#   
#   #根据IBI值评估EQ
#   test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
#   test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
#   
#   #lm线性拟合
#   fit1 <- lm(IBI ~ predict, data = test)
#   summary(fit1)
#   R2 <- c(R2,summary(fit1)$adj.r.squared)
#   #kappa评估
#   m <- kappa2(test[3:4])
#   kappa <- c(kappa,m$value)
  
  importance_otu <- data.frame(importance(rf.train), check.names = FALSE)
  importance_otu <- importance_otu[order(importance_otu$IncNodePurity, decreasing = TRUE), ]
  #然后取出排名靠前的 OTU，例如 top10 最重要的 OTU
  importance_otu.select <- importance_otu[1:100, ]
  otu_id.select <- rownames(importance_otu.select)
  otu.select <- asv_ibi[ ,c(otu_id.select, 'IBI')]
  
  trainset.select <- otu.select[trainlist,]
  testset.select <- otu.select[-trainlist,]
  seed <- c(seed,i)
  
  for (k in 1:10){
    set.seed(k)
    
    rf.train.select <- randomForest(IBI ~.,   #as.factor做分类分析，否则做回归分析
                                    data = trainset.select,importance = TRUE)
    # explain.select <- c(explain.select, mean(unlist(rf.train.select[5])))
    
    # 5. 根据建立的模型查看模型性能 
    set.seed(123)
    # 使用测试集，评估预测性能， 或
    predict.set <- testset.select
    predict <- predict(rf.train.select, predict.set) %>% as.data.frame()
    names(predict) = "predict"
    test <- cbind(subset(predict.set,select = IBI), subset(predict,select = c("predict"))) 
    
    # # #根据IBI值评估EQ
    #  test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
    #  test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
 
     test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
     test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)  
     m <- kappa2(test[3:4])
     kappa4_4 <- c(kappa4_4,m$value)
     test$IBI_EQ <- cut(test$IBI, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
     test$predict_EQ <- cut(test$predict, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
     m <- kappa2(test[3:4])
     kappa4_5 <- c(kappa4_5,m$value)
     max <- max(test$IBI)
     test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max * 0.75, 4.3), labels = c("1","2","3","4"), right=FALSE)
     test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max * 0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
     m <- kappa2(test[3:4])
     kappa4.1_4 <- c(kappa4.1_4,m$value)  
     test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
     test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
     m <- kappa2(test[3:4])
     kappa4.1_5 <- c(kappa4.1_5,m$value)
     max <- 4.388892073
     test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
     test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
     m <- kappa2(test[3:4])
     kappa4.3_4 <- c(kappa4.3_4,m$value)  
     test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
     test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
     m <- kappa2(test[3:4])
     kappa4.3_5 <- c(kappa4.3_5,m$value) 
     
    # 
    # #lm线性拟合
    # fit1 <- lm(IBI ~ predict, data = test)
    # R2.select <- c(R2.select,summary(fit1)$adj.r.squared)
    # #kappa评估
     # m <- kappa2(test[3:4])
     # kappa.select <- c(kappa.select,m$value)
     # 
     
     
    seed.select <- c(seed.select,paste(i,k,sep = "."))
  }
  
}

expression <- data.frame(seed.select,kappa4_4,kappa4_5,kappa4.1_4,kappa4.1_5,kappa4.3_4,kappa4.3_5)
expression
summary(expression)
expression.select <- data.frame(seed.select,rmse)
expression.select
summary(expression.select)



# 6. 循环计算变量筛选前RF结果与变量筛选后支持向量机结果 -----------------------------------------------------

seed <-  c()

explain <- c()
R2 <- c()
kappa <-  c()

R2.linear <- c()
kappa.linear <-  c()

R2.polynomial <- c()
kappa.polynomial <-  c()

R2.sigmoid <- c()
kappa.sigmoid <-  c()

R2.radial <- c()
kappa.radial <-  c()

rmse <- c()
rmse.linear <- c()
rmse.polynomial <- c()
rmse.sigmoid <- c()
rmse.radial <- c()
?svm

kappa4_4 <- c()
kappa4_5 <- c()
kappa4.1_4 <- c()
kappa4.1_5 <- c()
kappa4.3_4 <- c()
kappa4.3_5 <- c()


kappa.linear4_4 <- c()
kappa.linear4_5 <- c()
kappa.linear4.1_4 <- c()
kappa.linear4.1_5 <- c()
kappa.linear4.3_4 <- c()
kappa.linear4.3_5 <- c()

kappa.polynomial4_4 <- c()
kappa.polynomial4_5 <- c()
kappa.polynomial4.1_4 <- c()
kappa.polynomial4.1_5 <- c()
kappa.polynomial4.3_4 <- c()
kappa.polynomial4.3_5 <- c()

kappa.sigmoid4_4 <- c()
kappa.sigmoid4_5 <- c()
kappa.sigmoid4.1_4 <- c()
kappa.sigmoid4.1_5 <- c()
kappa.sigmoid4.3_4 <- c()
kappa.sigmoid4.3_5 <- c()

kappa.radial4_4 <- c()
kappa.radial4_5 <- c()
kappa.radial4.1_4 <- c()
kappa.radial4.1_5 <- c()
kappa.radial4.3_4 <- c()
kappa.radial4.3_5 <- c()





for (i in 1:100){
  
  # 4. 由训练集建立模型
  set.seed(i)
    rf.train <- randomForest(IBI ~.,   #as.factor做分类分析，否则做回归分析
                             data = trainset,importance = TRUE)
  # explain <- c(explain, mean(unlist(rf.train[5])))
   # rf.train <-  readRDS(paste("RDS_modle/rf.train.seed_",i,".rda",sep = ""))
  # # 5. 根据建立的模型查看模型性能
  # set.seed(123)
  # # 使用测试集，评估预测性能， 或
  predict.set <- testset
  predict <- predict(rf.train, predict.set) %>% as.data.frame()
  names(predict) = "predict"
  test <- cbind(subset(predict.set,select = IBI), subset(predict,select = c("predict")))

   test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
   test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)  
   m <- kappa2(test[3:4])
   kappa4_4 <- c(kappa4_4,m$value)
   test$IBI_EQ <- cut(test$IBI, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
   test$predict_EQ <- cut(test$predict, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
   m <- kappa2(test[3:4])
   kappa4_5 <- c(kappa4_5,m$value)
   max <- max(test$IBI)
   test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max * 0.75, 4.3), labels = c("1","2","3","4"), right=FALSE)
   test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max * 0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
   m <- kappa2(test[3:4])
   kappa4.1_4 <- c(kappa4.1_4,m$value)  
   test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
   test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
   m <- kappa2(test[3:4])
   kappa4.1_5 <- c(kappa4.1_5,m$value)
   max <- 4.388892073
   test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
   test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
   m <- kappa2(test[3:4])
   kappa4.3_4 <- c(kappa4.3_4,m$value)  
   test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
   test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
   m <- kappa2(test[3:4])
   kappa4.3_5 <- c(kappa4.3_5,m$value) 
   
  
  
  # #根据IBI值评估EQ
  # test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # 
  # 
  # 
  # #lm线性拟合
  # fit1 <- lm(IBI ~ predict, data = test)
  # summary(fit1)
  # R2 <- c(R2,summary(fit1)$adj.r.squared)
  # #kappa评估
  # m <- kappa2(test[3:4])
  # kappa <- c(kappa,m$value)
  # # 
  #  importance_otu <- data.frame(importance(rf.train), check.names = FALSE)
  #  importance_otu <- importance_otu[order(importance_otu$IncNodePurity, decreasing = TRUE), ]
  # # #然后取出排名靠前的 OTU，例如 top10 最重要的 OTU
  # importance_otu.select <- importance_otu[1:100, ]
  #  otu_id.select <- rownames(importance_otu.select)
  #  otu.select <- asv_ibi[ ,c(otu_id.select, 'IBI')]
  # # 
  #  trainset.select <- otu.select[trainlist,]
  # testset.select <- otu.select[-trainlist,]
 
   trainset.select <- trainset 
  testset.select <-testset
  seed <- c(seed,i)
  
  #常见参数：kernel-核函数:linear
  svm.train <- svm(IBI ~.,data = trainset.select,kernel = "linear",scale = F)
  
  predict.set <- testset.select
  predict <- predict(svm.train, predict.set) %>% as.data.frame()
  names(predict) = "predict"
  test <- cbind(subset(predict.set,select = IBI), subset(predict,select = c("predict"))) 
  
  test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)  
  m <- kappa2(test[3:4])
  kappa.linear4_4 <- c(kappa.linear4_4,m$value)
  test$IBI_EQ <- cut(test$IBI, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.linear4_5 <- c(kappa.linear4_5,m$value)
  max <- max(test$IBI)
  test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max * 0.75, 4.3), labels = c("1","2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max * 0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.linear4.1_4 <- c(kappa.linear4.1_4,m$value)  
  test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.linear4.1_5 <- c(kappa.linear4.1_5,m$value)
  max <- 4.388892073
  test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.linear4.3_4 <- c(kappa.linear4.3_4,m$value)  
  test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.linear4.3_5 <- c(kappa.linear4.3_5,m$value) 
  
  # #根据IBI值评估EQ
  # test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # 
  # #lm线性拟合
  # fit1 <- lm(IBI ~ predict, data = test)
  # R2.linear <- c(R2.linear,summary(fit1)$adj.r.squared)
  # #kappa评估
  # m <- kappa2(test[3:4])
  # kappa.linear <- c(kappa.linear,m$value)
  
  
  
  #常见参数：kernel-核函数:polynomial
  svm.train <- svm(IBI ~.,data = trainset.select,kernel = "polynomial",scale = F)
  
  predict.set <- testset.select
  predict <- predict(svm.train, predict.set) %>% as.data.frame()
  names(predict) = "predict"
  test <- cbind(subset(predict.set,select = IBI), subset(predict,select = c("predict"))) 
  
  test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)  
  m <- kappa2(test[3:4])
  kappa.polynomial4_4 <- c(kappa.polynomial4_4,m$value)
  test$IBI_EQ <- cut(test$IBI, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.polynomial4_5 <- c(kappa.polynomial4_5,m$value)
  max <- max(test$IBI)
  test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max * 0.75, 4.3), labels = c("1","2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max * 0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.polynomial4.1_4 <- c(kappa.polynomial4.1_4,m$value)  
  test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.polynomial4.1_5 <- c(kappa.polynomial4.1_5,m$value)
  max <- 4.388892073
  test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.polynomial4.3_4 <- c(kappa.polynomial4.3_4,m$value)  
  test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.polynomial4.3_5 <- c(kappa.polynomial4.3_5,m$value) 
  
  
  # #根据IBI值评估EQ
  # test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # 
  # #lm线性拟合
  # fit1 <- lm(IBI ~ predict, data = test)
  # R2.polynomial <- c(R2.polynomial,summary(fit1)$adj.r.squared)
  # #kappa评估
  # m <- kappa2(test[3:4])
  # kappa.polynomial <- c(kappa.polynomial,m$value)  
  
  
  
  #常见参数：kernel-核函数:sigmoid
  svm.train <- svm(IBI ~.,data = trainset.select,kernel = "sigmoid",scale = F)
  
  predict.set <- testset.select
  predict <- predict(svm.train, predict.set) %>% as.data.frame()
  names(predict) = "predict"
  test <- cbind(subset(predict.set,select = IBI), subset(predict,select = c("predict"))) 
  
  test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)  
  m <- kappa2(test[3:4])
  kappa.sigmoid4_4 <- c(kappa.sigmoid4_4,m$value)
  test$IBI_EQ <- cut(test$IBI, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.sigmoid4_5 <- c(kappa.sigmoid4_5,m$value)
  max <- max(test$IBI)
  test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max * 0.75, 4.3), labels = c("1","2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max * 0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.sigmoid4.1_4 <- c(kappa.sigmoid4.1_4,m$value)  
  test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.sigmoid4.1_5 <- c(kappa.sigmoid4.1_5,m$value)
  max <- 4.388892073
  test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.sigmoid4.3_4 <- c(kappa.sigmoid4.3_4,m$value)  
  test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.sigmoid4.3_5 <- c(kappa.sigmoid4.3_5,m$value) 
  
  # #根据IBI值评估EQ
  # test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # 
  # #lm线性拟合
  # fit1 <- lm(IBI ~ predict, data = test)
  # R2.sigmoid <- c(R2.sigmoid,summary(fit1)$adj.r.squared)
  # #kappa评估
  # m <- kappa2(test[3:4])
  # kappa.sigmoid <- c(kappa.sigmoid,m$value)  
  
  #常见参数：kernel-核函数:radial
  svm.train <- svm(IBI ~.,data = trainset.select,kernel = "radial",scale = F)
  
  predict.set <- testset.select
  predict <- predict(svm.train, predict.set) %>% as.data.frame()
  names(predict) = "predict"
  test <- cbind(subset(predict.set,select = IBI), subset(predict,select = c("predict"))) 
  
  test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)  
  m <- kappa2(test[3:4])
  kappa.radial4_4 <- c(kappa.radial4_4,m$value)
  test$IBI_EQ <- cut(test$IBI, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0.8,1.6,2.4,3.2,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.radial4_5 <- c(kappa.radial4_5,m$value)
  max <- max(test$IBI)
  test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max * 0.75, 4.3), labels = c("1","2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max * 0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.radial4.1_4 <- c(kappa.radial4.1_4,m$value)  
  test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.radial4.1_5 <- c(kappa.radial4.1_5,m$value)
  max <- 4.388892073
  test$IBI_EQ <- cut(test$IBI, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(0,max *0.25,max *0.5,max *0.75,4.3), labels = c("1","2","3","4"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.radial4.3_4 <- c(kappa.radial4.3_4,m$value)  
  test$IBI_EQ <- cut(test$IBI, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  test$predict_EQ <- cut(test$predict, breaks = c(max *0.2,max *0.4,max *0.6,max *0.8,4.3), labels = c("2","3","4","5"), right=FALSE)
  m <- kappa2(test[3:4])
  kappa.radial4.3_5 <- c(kappa.radial4.3_5,m$value) 
  
  
  
  # #根据IBI值评估EQ
  # test$IBI_EQ <- cut(test$IBI, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # test$predict_EQ <- cut(test$predict, breaks = c(1,2,3,4.3), labels = c("2","3","4"), right=FALSE)
  # 
  # #lm线性拟合
  # fit1 <- lm(IBI ~ predict, data = test)
  # R2.radial <- c(R2.radial,summary(fit1)$adj.r.squared)
  # #kappa评估
  # m <- kappa2(test[3:4])
  # kappa.radial <- c(kappa.radial,m$value)
  # 
  
  
    

}

seed=c(1:100)

expression <- data.frame(seed,kappa4_4,kappa4_5,kappa4.1_4,kappa4.1_5,kappa4.3_4,kappa4.3_5)
expression
summary(expression)

expression.svm <- data.frame(seed,kappa.linear4_4,kappa.polynomial4_4,kappa.sigmoid4_4,kappa.radial4_4,
                                  kappa.linear4_5,kappa.polynomial4_5,kappa.sigmoid4_5,kappa.radial4_5,
                                  kappa.linear4.1_4,kappa.polynomial4.1_4,kappa.sigmoid4.1_4,kappa.radial4.1_4,
                                  kappa.linear4.1_5,kappa.polynomial4.1_5,kappa.sigmoid4.1_5,kappa.radial4.1_5,
                                  kappa.linear4.3_4,kappa.polynomial4.3_4,kappa.sigmoid4.3_4,kappa.radial4.3_4,
                                  kappa.linear4.3_5,kappa.polynomial4.3_5,kappa.sigmoid4.3_5,kappa.radial4.3_5)
expression.svm
summary(expression.svm)

write.csv(expression,"变量选择后RF.kappa结果.csv")





# 6. 循环统计100次筛选的重要otu -----------------------------------------------------

for (i in 1:100){
  
  # 4. 由训练集建立模型 -------------------------------------------------------------
  # set.seed(i)
  # rf.train <- randomForest(IBI ~.,   #as.factor做分类分析，否则做回归分析
  #                          data = trainset,importance = TRUE)
  rf.train <- readRDS(paste("RDS/replicate.seed_",i,".rda",sep = "."))
   importance_otu <- data.frame(importance(rf.train), check.names = FALSE)
  

  importance_otu <- importance_otu[order(importance_otu$IncNodePurity, decreasing = TRUE), ]
  #然后取出排名靠前的 OTU，例如 top10 最重要的 OTU
  importance_otu.select <- importance_otu[1:100, ]
  assign(paste("replicate",i,sep = "."),rownames(importance_otu.select))
  
}
importance = as.data.frame(mget(paste("replicate",1:3,sep = ".")))
importance


# 7. 循环保存模型 ---------------------------------------------------------------

for (i in 1:100){
  
  # 4. 由训练集建立模型 -------------------------------------------------------------
  set.seed(i)
  rf.train <- randomForest(IBI ~.,   #as.factor做分类分析，否则做回归分析
                           data = trainset,importance = TRUE)
  
  saveRDS(rf.train, paste("RDS_modle/rf.train.seed_",i,".rda",sep = ""))
  
}
















