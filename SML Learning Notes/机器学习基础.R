# 随机森林模型，既可做分类预测，也可做回归预测 --------------------------------------------------

# 1. 随机森林模型包的载入 -----------------------------------------------------------------
library(randomForest)
library(caret)
library(pROC)

# 2. 数据载入 -----------------------------------------------------------------
data("iris")
summary(iris)

# 3. 划分训练、测试集 -------------------------------------------------------------
dim(iris)  #查看数据多少行，多少列
trainlist <- createDataPartition(iris$Species,p=3/4,list=FALSE) #将数据集划分为2部分
trainset <- iris[trainlist,]
testset <- iris[-trainlist,]

# 4. 由训练集建立模型 -------------------------------------------------------------
set.seed(123)
rf.train <- randomForest(as.factor(Species) ~.,   #as.factor做分类分析，否则做回归分析
                         data = trainset,
                         importance = TRUE,na.action = na.pass)
rf.train
plot(rf.train,main="randomforest origin")

#变量重要性排序
importance(rf.train)
varImpPlot(rf.train)

# 5. 根据建立的模型进行预测 ----------------------------------------------------------
set.seed(123)
rf.test <- predict(rf.train,newdata = testset,type = "class")
rf.test

#统计/对比预测结果
rf.cf <- caret::confusionMatrix(as.factor(rf.test),as.factor(testset$Species)) 
rf.cf


# 6. ROC,AUC曲线绘制 ----------------------------------------------------------
rf.test2 <- predict(rf.train,newdata = testset,type = "prob") #ROC需要概念而不是种类
head(rf.test2)

#总的ROC,不能绘制曲线
roc.rf <- multiclass.roc(testset$Species,rf.test2) #多种类的ROC，原始标签+预测数据
roc.rf

#必须选择一类变量，才能绘制曲线
roc.rf_1 <- multiclass.roc(testset$Species,rf.test2[,1]) #选择第一类

plot(roc.rf_1$rocs[[1]],col='blue')
plot(roc.rf_1$rocs[[3]],add= TRUE,col='red')

#计算auc
auc(roc.rf_1)


#method 1.ggplot2
library(ggplot2)
library(ROCR) 
data(ROCR.simple) 
ROCR.simple$predictions
ROCR.simple$labels
summary(ROCR.simple)

pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels)  
perf <- performance(pred,"tpr","fpr") 
perf

x <- unlist(perf@x.values)
y <- unlist(perf@y.values)
plotdata <- data.frame(x,y) 
names(plotdata) <- c("x", "y")

g <- ggplot(plotdata) + 
  geom_path(aes(x = x, y = y, colour = x), size=1) + 
  labs(x = "False positive rate", y = "True positive rate", title ="ROC Curves") +
  scale_colour_gradient(name = 'False positive rate', low = 'blue', high = 'red') +
  theme(plot.title = element_text(face = 'bold',size=15))
g

#method 2. pROC
library(pROC)##roc
data(aSAH)

roc1<-roc(aSAH$outcome,aSAH$age)
roc2<-roc(aSAH$outcome,aSAH$s100b)
plot(roc1,col='blue')
plot.roc(roc2,add=TRUE,col='red')
auc(roc1, partial.auc = c(1, .9))
smooth(roc1)



# 决策树模型,主要做分类预测 -----------------------------------------------------------


# 1. 决策树模型包的载入 ------------------------------------------------------------
library(rpart)
library(rpart.plot)


# 2. 数据载入 -----------------------------------------------------------------
data("iris")
summary(iris)

# 3. 划分训练、测试集 -------------------------------------------------------------
dim(iris)  #查看数据多少行，多少列
 #从总行数里随机抽取一定比例的行为训练集，replace=F表示不重复抽
s=sample(c(1:nrow(iris)),120,replace = FALSE)
s
trainset <- iris[s,]
testset <- iris[-s,]

# 4. 由训练集建立模型 -------------------------------------------------------------
# set.seed(123)

fit1 <- rpart(Species ~ .,trainset)
summary(fit1)  #结果很难读


# 5. plot简洁看结果 ------------------------------------------------------------

rpart.plot(fit1,type = 2)

# 6. 根据建立的模型进行预测 ----------------------------------------------------------
# set.seed(123)
pre <- predict(fit1,newdata = testset,type = "class")
pre <- predict(fit1,newdata = testset)
pre

t <- table(pre,testset$Species) #列个表，比较结果
t
#acc计算准确率
acc=sum(diag(t))/nrow(testset)*100
acc




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
#常见参数：kernel-核函数
fit1 <- svm(quality ~.,data = trainset,kernel = "linear")
summary(fit1)
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



# 分类模型构建和多分类ROC曲线绘制,以决策树模型为例 -------------------------------------------------------
# 1. 包的载入 ------------------------------------------------------------
library(ggplot2)
library(ROCR)
library(rpart)
library(rpart.plot)
library(caret)

# 2. 数据载入 -----------------------------------------------------------------
data("iris")
attach(iris)
summary(iris)

# 3. 划分训练、测试集 -------------------------------------------------------------
dim(iris)  #查看数据多少行，多少列

#方法1：sample
#从总行数里随机抽取一定比例的行为训练集，replace=F表示不重复抽
s=sample(c(1:nrow(iris)),100,replace = FALSE)
s
trainset <- iris[s,]
testset <- iris[-s,]

#方法2：caret包createDataPartition,只需指定因变量
trainlist <- createDataPartition(iris$Species,times=1, p=0.5,list=FALSE) #将数据集划分为2部分
trainset <- iris[trainlist,]
testset <- iris[-trainlist,]


# 4. 由训练集建立模型 -------------------------------------------------------------
# set.seed(123)

fit1 <- rpart(Species ~ .,trainset)
summary(fit1)  #结果很难读

# 5. plot简洁看结果 ------------------------------------------------------------

rpart.plot(fit1,type = 2)

# 6. 根据建立的模型进行预测 ----------------------------------------------------------
# set.seed(123)
pre <- predict(fit1,newdata = testset,type = "class")
pre <- predict(fit1,newdata = testset)
pre

t <- table(pre,testset$Species) #列个表，比较结果
t
#acc计算准确率
acc=sum(diag(t))/nrow(testset)*100
acc

# 7. ROC曲线绘制 --------------------------------------------------------------
prep2 <- predict(fit1,testset,type = "prob") #注意绘制曲线一定是概率而不是分类标签

roc1 <- multiclass.roc(testset$Species,prep2[,1])
roc1
plot(roc1$rocs[[1]], col = 'blue')
plot.roc(roc1$rocs[[3]],add = TRUE,col = 'red')

#计算auc值
auc(roc1)


# 逻辑回归（二分类问题）的十折交叉验证 -------------------------------------------------------------
# 1. 包的载入 --------------------------------------------------------------------
library(caret)

# 2. 数据载入 -----------------------------------------------------------------
data("iris")
summary(iris)
#二分类逻辑回归问题
iris = iris[iris$Species != "virginica",]

iris$Species <- factor(iris$Species,labels=c('setosa','versicolor'))
summary(iris$Species)

# 3. 划分训练、测试集 -------------------------------------------------------------
dim(iris)  #查看数据多少行，多少列
#从总行数里随机抽取一定比例的行为训练集，replace=F表示不重复抽
s=sample(c(1:nrow(iris)),70,replace = FALSE)
s
trainset <- iris[s,]
testset <- iris[-s,]

# 4. 逻辑回归训练 ------------------------------------------------------------------
modle_logit <- glm(Species ~ .,data=trainset,family = binomial(link = 'logit'))
summary(modle_logit)


# 5. 交叉验证 -----------------------------------------------------------------
set.seed(123) #固定folds函数的分组

#方法1：for循环
folds <- createFolds(y=iris$Species,k=10)
folds 

max=0  #最高准确率
num=0  #最佳一折

for (i in 1:10) {
  print(i)
}

#方法2：10-folds cv
set.seed(123)
fit= train(Species ~ . , data=iris, method = "glm",
           trControl =trainControl(method = "cv", number = 10))
fit


# 线性回归 --------------------------------------------------------------------

# 1. 数据载入 -----------------------------------------------------------------
attach(cars)
summary(cars)


# 2. 简单线性回归 ---------------------------------------------------------------
#lm拟合
fit1 <- lm(dist ~ speed, data = cars)
summary(fit1)

#plot画图
plot(cars$speed,cars$dist)
abline(fit1)

#查看预测效果
fitted(fit1)  #预测值
residuals(fit1)  #残差


# 3. 多项式回归 ----------------------------------------------------------------
fit2 <- lm(dist ~ speed +I(speed ^ 2), data = cars)
summary(fit2)

#plot画图
plot(cars$speed,cars$dist)
abline(fit2)

#car包scatterplot画图
library(car)
scatterplot(dist ~ speed, data = cars,spread= FALSE,
            main="cars",xlab = "speed",ylab = "dist")


# 4. 多元线性回归 ---------------------------------------------------------------
attach(ChickWeight)
summary(ChickWeight)

#先看变量之间的相关性，如果不相关，就没必要做回归
cor(ChickWeight$weight,ChickWeight$Time)
scatterplotMatrix(ChickWeight,spread=F,main="chik")

fit3 <- lm(weight ~ Time + Chick + Diet,data = ChickWeight) #diet分类变量
fit3 <- lm(weight ~ Time + Chick,data = ChickWeight)
summary(fit3)



# 主成分分析（三种方法） -------------------------------------------------------------------
# 1. 数据载入 -----------------------------------------------------------------
data("USArrests")
str(USArrests)
summary(USArrests)

# 2方法1：prcomp --------------------------------------------------------------
#prcomp(formula, data = NULL, subset, na.action, ...)
prcomp(USArrests)  #没有scale

prcomp(USArrests, scale = TRUE) #ֱscale标准化

p <- prcomp(~ Murder + Assault + Rape, data = USArrests, scale = TRUE) #直接方程
p
plot(p)
plot(prcomp(USArrests))

summary(prcomp(USArrests, scale = TRUE))

biplot(prcomp(USArrests, scale = TRUE))



# 3. 方法2：princomp -----------------------------------------------------------
#princomp(formula, data = NULL, subset, na.action, ...) #??????ȫһ??
princomp(USArrests, cor = TRUE) 
pc.cr <- princomp(USArrests, cor = TRUE)
summary(pc.cr)

loadings(pc.cr)  #rotation各个主成分得分


plot(pc.cr) # screeplot

biplot(pc.cr)


# 方法3：psych包 --------------------------------------------------------------
library(psych) 

df <- USArrests
df.cor <- cor(df) # 相关性分析
df.cor  

#自动判断主成分个数 + 碎石图
fa.parallel(df, fa = "pc", n.iter = 100,
            show.legend = F, main = "Scree plot with parallel analysis")

#根据上一步判断，选2个主成分
pc<-principal(df, nfactors = 2, score = T, rotate = "varimax") 

summary(pc)

pc$loadings

round(unclass(pc$weights),2)


