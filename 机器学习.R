# 随机森林模型，既可做分类预测，也可做回归预测 --------------------------------------------------

# 1. 随机森林模型包的载入 -----------------------------------------------------------------
library(randomForest)
library(caret)
library(pROC)
library(tidyverse)
library(irr)
library(Metrics)
# 2. 数据载入 -----------------------------------------------------------------
# data("iris")
# summary(iris)
asv <- read.csv("../all_data/otu_CSS_table_84sample_combine9396.csv",header = TRUE, row.names = 1,sep = ",",check.names = F)  %>% as.data.frame()
IBI <- read.csv("../all_data/IBI.csv",header = TRUE, row.names = 1,sep = ",",check.names = F)  %>% as.data.frame()

core <- read.csv("importance_otu_100seed.csv",header = TRUE,sep = ",",check.names = F)  %>% as.data.frame()
a <- core$Core[1:21]
asv_ibi <- cbind(t(asv[a,]),IBI) %>% as.data.frame() 
tail(names(asv_ibi))
# names(asv_ibi)[9397] <- "y"

# 3. 划分训练、测试集 -------------------------------------------------------------
dim(asv_ibi)  #查看数据多少行，多少列
set.seed(2)
trainlist <- createDataPartition(asv_ibi$IBI,p=2/3,list=FALSE) #将数据集划分为2部分
trainset <- asv_ibi[trainlist,]
testset <- asv_ibi[-trainlist,]

testset[,"IBI"]


# 4. 由训练集建立模型 -------------------------------------------------------------
set.seed(70)
rf.train <- randomForest(IBI ~.,   #as.factor做分类分析，否则做回归分析
                         data = trainset,importance = TRUE)
rf.train
summary(rf.train)
rf.train[5] %>% unlist() %>% mean()
rf.train$tr(rf.train$rsq)
names(rf.train)

plot(rf.train,main="randomforest origin")
which.min(rf.train$mse)

# 5. 根据建立的模型查看模型性能 ----------------------------------------------------------
set.seed(123)
# 使用训练集，查看预测精度，或
# 使用测试集，评估预测性能， 或
# 使用全集，查看整体情况


# predict.set <- testset
# predict.set <- asv_ibi
# rf.train <- fit1

predict <- predict(rf.train, predict.set) %>% as.data.frame()
names(predict) = "predict"


rownames(predict)
test <- cbind(subset(predict.set,select = IBI), subset(predict,select = c("predict"))) 

#根据IBI值评估EQ
max = 4.179795479
  
test$IBI_EQ <- cut(test$IBI, breaks = c(0, max *0.25,max *0.5,max * 0.75, max), 
                   labels = c("1","2","3","4"), right=FALSE)
test$predict_EQ <- cut(test$predict,breaks = c(0, max *0.25,max *0.5,max * 0.75, max), 
                       labels = c("1","2","3","4"), right=FALSE)
head(test)
str(test)

rmse(test$IBI,test$predict)
#lm线性拟合
fit2 <- lm(IBI ~ predict, data = test)
summary(fit2)
summary(fit2)$adj.r.squared
plot(test$IBI,test$predict, main = '测试集', 
     xlab = 'IBI', ylab = 'Predict')
abline(0,1)
abline(fit2)
# ?abline

#kappa评估
kappa2(test[3:4])

write.csv(test,"SVM.21.csv")



# 6. OTU 的重要性评估 --------------------------------------------------------------
##OTU 的重要性评估
#查看表示每个预测变量（细菌 OTU）重要性的得分
# summary(rf.train)
# importance_otu <- rf.train$importance
# head(importance_otu)
# 
# #或者使用函数 importance()
 importance_otu <- data.frame(importance(rf.train), check.names = FALSE)
 head(importance_otu)

#作图展示 top30 重要的 OTUs
varImpPlot(rf.train, n.var = min(30, nrow(rf.train$importance)), 
           main = 'Top 30 - variable importance')

#可以根据某种重要性的高低排个序，例如根据“IncNodePurity”指标
importance_otu <- importance_otu[order(importance_otu$IncNodePurity, decreasing = TRUE), ]
head(importance_otu,10)

#输出表格
#write.table(importance_otu, 'importance_Core.feature.csv', sep = ',', col.names = NA, quote = FALSE)

# 7. 交叉验证辅助评估选择特定数量的 OTU --------------------------------------------------
#5 次重复十折交叉验证
set.seed(123)
otu_train.cv <- replicate(5, rfcv(trainset[-ncol(trainset)], trainset$IBI, cv.fold = 10, step = 1.2), simplify = FALSE)
otu_train.cv[1] %>% head()

#提取验证结果绘图
otu_train.cv <- data.frame(sapply(otu_train.cv, '[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
# otu_train.cv.mean <- read.csv("otu_train.cv.mean.csv",row.names = 1)




#拟合线图
library(ggplot2)
library(ggbreak)
p <- ggplot(otu_train.cv.mean, aes(Group.1, x)) +
  geom_line() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +  
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')
p+scale_x_break(breaks = c(0,1))+scale_x_break(breaks = c(400,9000))




#提示保留109个 ? 重要的 OTU，可以使随机森林回归的精度最大化
#首先根据某种重要性的高低排个序，例如根据“IncNodePurity”指标
importance_otu <- importance_otu[order(importance_otu$IncNodePurity, decreasing = TRUE), ]

#然后取出排名靠前的 OTU，例如 top10 最重要的 OTU
importance_otu.select <- importance_otu[1:100, ]
importance_otu.select

#输出表格
#write.table(importance_otu.select, 'importance_otu.select.txt', sep = '\t', col.names = NA, quote = FALSE)

# #简单查看重要的OTU 丰度与响应变量的关系
#可以看到趋势非常明显，包括根际富集或排斥等都有涉及
otu_id.select <- rownames(importance_otu.select)
otu.select <- asv_ibi[ ,c(otu_id.select, 'IBI')]


otu.select <- reshape2::melt(otu.select, id = 'IBI')
ggplot(otu.select, aes(x = IBI, y = value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~variable, ncol = 3, scale = 'free_y') +
  labs(title = '',x = 'IBI', y = 'CSS_normalized')





??predict
rf.test <- predict(rf.train,newdata = testset,proximity = T)

rf.test

#统计/对比预测结果
rf.cf <- caret::confusionMatrix(as.factor(rf.test),as.factor(testset$Species)) 
rf.cf

rf.resid <- rf.test - testset$y #calculate residual
rf.resid
mean(rf.resid^2)

test <- cbind(rf.test,testset$IBI) %>% as.data.frame()

plot(test$rf.test,test$V2)

fit1 <- lm(V2 ~ rf.test, data = test)
summary(fit1)

#plot画图
plot(test$rf.test,test$V2)
abline(fit1)



# 4. 循环预测 --------------------------------------------------------------------


for (i in 1){
  set.seed(i)
  rf.train <- randomForest(IBI ~.,   #as.factor做分类分析，否则做回归分析
                           data = trainset
                           # ,scale = F
  )
  summary(rf.train)
  assign(paste("rf.test",i,sep = "."),
         predict(rf.train,newdata = testset))
  
}

predict.test = as.data.frame(mget(paste("rf.test",1,sep = "."))) %>%  as.data.frame()
predict.test <- test
predict.test$mean <- apply(predict.test[1], 1, mean)
predict.test$sd <- apply(predict.test[1], 1, sd)
predict.test$sd <- rep(0,42)
predict.test$group <- rep(c("A","B","C"),14)
head(predict.test)

test <- cbind(subset(testset,select = IBI), subset(predict.test,select = c("mean","sd","group"))) 
head(test)
str(test)

plot(test$IBI,test$mean)


fit1 <- lm(IBI ~ mean, data = test)
summary(fit1)



# 5. ggplot2画图  -------------------------------------------------------------


errorbar_ORG_adj <- predict.test
names(errorbar_ORG_adj) =c("AMBI_g","adj_ambi","sd_ORG","Farm")


plot_AMBI_ORG_adj<-ggplot(data = errorbar_ORG_adj,aes(x = AMBI_g,y = adj_ambi, group=Farm)) +
  scale_x_continuous(name="Predict", limits=c(0, 4), breaks = c(0,1,2,3,4)) +
  scale_y_continuous(name=" Observe", limits=c(0, 4), breaks = c(0,1,2,3,4)) +
  geom_point(aes(x=AMBI_g, y=adj_ambi, shape=Farm, color=Farm), size=1)+
  scale_shape_manual(values=c(16,16,16,16,16,16,16))+
  scale_color_manual(values=c( "green3","cyan", "magenta", "yellow"))+
  geom_smooth(data=errorbar_ORG_adj,aes(x=AMBI_g, y=adj_ambi),method="lm",inherit.aes=F) +
  geom_errorbar(data=errorbar_ORG_adj,aes(ymin = adj_ambi-sd_ORG, ymax = adj_ambi+sd_ORG, color=Farm), width=0.1)+
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, alpha = .2)+
  annotate("rect", xmin = 1, xmax = 2, ymin = 1, ymax = 2,  alpha = .2)+
  annotate("rect", xmin = 2, xmax = 3, ymin = 2, ymax = 3, alpha = .2)+ 
  annotate("rect", xmin = 3, xmax = 4, ymin = 3, ymax = 4, alpha = .2)+ 
  theme_classic()


plot_AMBI_ORG_adj

?randomForest


# 支持向量机模型,主要做分类预测 ---------------------------------------------------------

# 1. 支持向量机（SVM）包的载入 ------------------------------------------------------------
library(e1071)

# 2. 数据载入 -----------------------------------------------------------------
#红酒质量数据
redwine <- read.csv("winequality-red.csv",header = T,sep = ";")
summary(redwine)

redwine$quality <- as.factor(redwine$quality) #将待预测指标转换为因子
summary(redwine$quality)

# 3. 划分训练、测试集 -------------------------------------------------------------
dim(redwine)  #查看数据多少行，多少列
s <- sample(nrow(redwine),nrow(redwine)*0.7,replace = F)
length(s)

trainset <- redwine[s,]
testset <- redwine[-s,]
dim(testset)

# 4. 由训练集建立模型 -------------------------------------------------------------
# set.seed(123)
#常见参数：kernel-核函数??svm
fit1 <- svm(IBI ~.,data = trainset,kernel = "radial",scale = F)
summary(fit1) %>% names() 
summary(fit1)$kernel

fit2 <- svm(quality ~.,data = trainset,kernel = "sigmoid")
summary(fit2)

# 5. plot简洁看结果 ------------------------------------------------------------
plot(fit1,trainset,sulphates~density)
plot(fit2,trainset,sulphates~density)

# 6. 根据建立的模型进行预测 ----------------------------------------------------------
# set.seed(123)
p1 <- predict(fit1,testset)
p1
p2 <- predict(fit2,testset)
p2


# 7. 预测结果分析 ---------------------------------------------------------------

#列个表，比较结果
t1 <- table(p1,testset$quality)
t1

#acc计算准确率
acc1 <- sum(diag(t1))/nrow(testset)
acc1

t2 <- table(p2,testset$quality)
acc2 <- sum(diag(t2))/nrow(testset)
acc2





