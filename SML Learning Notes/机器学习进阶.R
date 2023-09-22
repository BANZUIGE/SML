# ##随机森林模型 ------------------------------------------------------------------
# 1. 加载程序包 ----------------------------------------------------------------
library(rpart) #classification and regression trees
library(partykit) #treeplots
library(MASS) #breast and pima indian data
library(ElemStatLearn) #prostate data
library(randomForest) #random forests
library(xgboost) #gradient boosting 
library(caret) #tune hyper-parameters



################RF

# 2. 随机森林回归模型构建 ----------------------------------------------------------------
#  划分数据集
data(prostate)
prostate$gleason <- ifelse(prostate$gleason == 6, 0, 1)
pros.train <- subset(prostate, train == TRUE)[, 1:9]
pros.test = subset(prostate, train == FALSE)[, 1:9]
set.seed(123)
rf.pros <- randomForest(lpsa ~ ., data = pros.train)
rf.pros

#查看误差
plot(rf.pros)
which.min(rf.pros$mse)  #80个树的时候误差最小

#根据上面指定树的个数，即80个树的时候，均方误差稍微降低了
set.seed(123)
rf.pros.2 <- randomForest(lpsa ~ ., data = pros.train, ntree = 75)
rf.pros.2

#查看影响因变量的自变量重要性
varImpPlot(rf.pros.2, scale = TRUE,
           main = "Variable Importance Plot - PSA Score")
importance(rf.pros.2)


# 4. 查看模型在验证集上的表现 ---------------------------------------------------------
rf.pros.test <- predict(rf.pros.2, newdata = pros.test)

#查看均方误差
#plot(rf.pros.test, pros.test$lpsa)
rf.resid <- rf.pros.test - pros.test$lpsa #calculate residual
rf.resid
mean(rf.resid^2)


# 5. 随机森林分类模型构建 -----------------------------------------------------------
#划分数据集
data(biopsy)
biopsy <- biopsy[, -1]
names(biopsy) <- c("thick", "u.size", "u.shape", "adhsn", "s.size", "nucl", "chrom", "n.nuc", "mit", "class")
biopsy.v2 <- na.omit(biopsy)
set.seed(123) #random number generator
ind <- sample(2, nrow(biopsy.v2), replace = TRUE, prob = c(0.7, 0.3))
biop.train <- biopsy.v2[ind == 1, ] #the training data set
biop.test <- biopsy.v2[ind == 2, ] #the test data set
str(biop.test)

#建立模型
set.seed(123)
rf.biop <- randomForest(class ~ ., data = biop.train)
rf.biop

#画图
plot(rf.biop)
which.min(rf.biop$err.rate[, 1])   #查看均方误差最小时树的个数

set.seed(123)
rf.biop.2 <- randomForest(class ~ ., data = biop.train, ntree = 125)
#getTree(rf.biop,1)
rf.biop.2


# 6. 验证集验证 ----------------------------------------------------------------
rf.biop.test <- predict(rf.biop.2, 
                        newdata = biop.test, 
                        type = "response")   #response反映率

table(rf.biop.test, biop.test$class)
(138 + 67) / 209    #计算误差率，对角线之和/总数

#查看变量重要性
varImpPlot(rf.biop.2)
importance(rf.biop.2)


#其他数据集
data(Pima.tr)
data(Pima.te)
pima <- rbind(Pima.tr, Pima.te)

#划分训练集，验证集
set.seed(502)
ind <- sample(2, nrow(pima), replace = TRUE, prob = c(0.7, 0.3))
pima.train <- pima[ind == 1, ]
pima.test <- pima[ind == 2, ]

#建模
set.seed(321)
rf.pima <- randomForest(type ~ ., data = pima.train)
rf.pima

# plot(rf.pima)
which.min(rf.pima$err.rate[,1])

set.seed(321)
rf.pima.2 <- randomForest(type ~ ., data = pima.train, ntree = 80)
rf.pima.2
rf.pima.test <- predict(rf.pima.2, 
                        newdata = pima.test, 
                        type = "response")
table(rf.pima.test, pima.test$type)
(76+33)/147
#varImpPlot(rf.pima.2)


# ##支持向量机模型 ---------------------------------------------------------------
##
library(class)
library(kknn)
library(e1071)
library(kernlab)
library(caret)
library(MASS)
library(reshape2)
library(ggplot2)
library(pROC)

#划分数据集
data(Pima.tr)
str(Pima.tr)
data(Pima.te)
str(Pima.te)
pima <- rbind(Pima.tr, Pima.te)
str(pima)

pima.melt <- melt(pima, id.var = "type")
pima.scale <- data.frame(scale(pima[, -8]))
#scale.pima = as.data.frame(scale(pima[,1:7], byrow=FALSE)) #do not create own function
str(pima.scale)
pima.scale$type <- pima$type

pima.scale.melt <- melt(pima.scale, id.var = "type")
cor(pima.scale[-8])
table(pima.scale$type)
set.seed(502)
ind <- sample(2, nrow(pima.scale), replace = TRUE, prob = c(0.7, 0.3))
train <- pima.scale[ind == 1, ]
test <- pima.scale[ind == 2, ]
str(train)
str(test)


#线性SVM
set.seed(123)
#turn.svm() 使用交叉验证使调优参数达到最优
linear.tune <- tune.svm(type ~ ., data = train, 
                        kernel = "linear",    #线性核函数
                        cost = c(0.001, 0.01, 0.1, 1, 5, 10))
summary(linear.tune)  #10folds交叉验证

best.linear <- linear.tune$best.model  #提取最好的模型

#在测试集上进行验证
tune.test <- predict(best.linear, newdata = test)
tune.test 
table(tune.test, test$type)
(80+32)/147   #准确率计算



#多项核函数
set.seed(123)
poly.tune <- tune.svm(type ~ ., data = train, 
                      kernel = "polynomial",   
                      degree = c(3, 4, 5),   #多项式的阶
                      coef0 = c(0.1, 0.5, 1, 2, 3, 4))  #核系数
summary(poly.tune)

#提取最优模型
best.poly <- poly.tune$best.model
poly.test <- predict(best.poly, newdata = test)
table(poly.test, test$type)
(81 + 26) / 147

#雷达核函数
set.seed(123)
rbf.tune <- tune.svm(type ~ ., data = train, 
                     kernel = "radial", 
                     gamma = c(0.1, 0.5, 1, 2, 3, 4))
summary(rbf.tune)
best.rbf <- rbf.tune$best.model
rbf.test <- predict(best.rbf, newdata = test)
table(rbf.test, test$type)
(73+21)/147

#sigmoid核函数
set.seed(123)
sigmoid.tune <- tune.svm(type ~ ., data = train, 
                         kernel = "sigmoid", 
                         gamma = c(0.1, 0.5, 1, 2, 3, 4),
                         coef0 = c(0.1, 0.5, 1, 2, 3, 4))
summary(sigmoid.tune)

best.sigmoid <- sigmoid.tune$best.model
sigmoid.test <- predict(best.sigmoid, newdata = test)
table(sigmoid.test, test$type)
(74+31)/147


#caret包confusionMatrix混淆矩阵函数计算更多信息
confusionMatrix(sigmoid.test, test$type, positive = "Yes")
confusionMatrix(tune.test, test$type, positive = "Yes")


# 4. SVM中的特征选择 ------------------------------------------------------------

set.seed(123)
#建立模型
rfeCNTL <- rfeControl(functions = lrFuncs, method = "cv", number = 10) #定义验证方法：10折
svm.features <- rfe(train[, 1:7], train[, 8],
                    sizes = c(7, 6, 5, 4),  #指定变量数量
                    rfeControl = rfeCNTL, 
                    method = "svmLinear")
svm.features  #筛选出了5个重要特征

#利用上面筛选出来的特征建模
svm.5 <- svm(type ~ glu + ped + npreg + bmi + age, 
             data = train, 
             kernel = "linear")

#验证对应的5个变量
svm.5.predict = predict(svm.5, newdata=test[c(1,2,5,6,7)])
table(svm.5.predict, test$type)
(79 + 33)/147
#caret包confusionMatrix混淆矩阵函数计算更多信息
confusionMatrix(svm.5.predict, test$type, positive = "Yes")




# 主成分分析 -------------------------------------------------------------------
library(psych)

#碎石图，判断主成分个数
fa.parallel(USJudgeRatings[,-1], fa="pc", n.iter=100,
            show.legend=FALSE, main="Scree plot with parallel analysis")
abline(h=1,lwd=1,col="green") 

# Listing 14.1 - Principal components analysis of US Judge Ratings(原始数据集)
#library(psych)
pc <- principal(USJudgeRatings[,-1], nfactors=1)
pc

# Principal components analysis Harman23.cor data
#library(psych)
fa.parallel(Harman23.cor$cov, n.obs=302, fa="pc", n.iter=100,
            show.legend=FALSE, main="Scree plot with parallel analysis")
abline(h=1,lwd=1,col="green")

# Listing 14.2 - Principal components analysis of body measurements(相关系数矩阵)
#library(psych)
PC <- principal(Harman23.cor$cov, nfactors=2, rotate="none")
PC

# Listing 14.3 - Principal components analysis with varimax rotation
rc <- principal(Harman23.cor$cov, nfactors=2, rotate="varimax")   #旋转
rc

# Listing 14.4 - Obtaining componenet scores from raw data
library(psych)
pc <- principal(USJudgeRatings[,-1], nfactors=1, score=TRUE)
head(pc$scores)
cor(USJudgeRatings$CONT, pc$score)

# Listing 14.5 - Obtaining principal component scoring coefficients
library(psych)
rc <- principal(Harman23.cor$cov, nfactors=2, rotate="varimax")
round(unclass(rc$weights), 2)


# 主成分分析在线性回归中的应用 ----------------------------------------------------------

# 案例3
library(car)
example16_2  <- read.table ("example16_2.csv", header=TRUE, sep=",")
example16_2
fit <- lm(y~x1+x2+x3, data=example16_2)
summary(fit)

#判断变量之间的共线性
vif(fit)


library(psych)

#描述数据集
describe(example16_2)

#碎石图，判断主成分个数
fa.parallel(example16_2[-4], fa="pc", n.iter=100, show.legend=FALSE, main="Screen plot with parallel analysis")
abline(1,0)

#提取主成分并旋转
pc <- principal(example16_2[-4], nfactors=2, rotate= "varimax", score=TRUE)
pc

#提取每个观测值的得分并加到数据框
pc$weights
pc$scores
newdata <- data.frame(example16_2,  pc$scores)
newdata

#与因变量进行拟合
fit <- lm(y~ RC1+RC2, data=newdata)
summary(fit)
vif(fit)


####
library(ggplot2) #support scatterplot
# library(GPArotation) #support rotation
library(psych) #PCA package

#载入数据
train <- read.csv("NHLtrain.csv")
str(train)
names(train)

#数据标准化并计算相关性
train.scale <- scale(train[, -1:-2])
nhl.cor = cor(train.scale)
# cor.plot(nhl.cor)

#主成分分析，直接分析并画图
pca <- principal(train.scale, rotate="none")   #默认只提取一个主成分？
pca
plot(pca$values, type="b", ylab="Eigenvalues", xlab="Component")

#用旋转的方法提取5个主成分
pca.rotate <- principal(train.scale, nfactors = 5, rotate = "varimax")
pca.rotate

#提取主成分得分并添加Y轴（因变量）
pca.scores <- data.frame(pca.rotate$scores)
head(pca.scores)
pca.scores$ppg <- train$ppg

#线性拟合
nhl.lm <- lm(ppg ~ ., data = pca.scores)   #所有主成分，解释69.81%的变异
summary(nhl.lm)

nhl.lm2 <- lm(ppg ~ RC1 + RC2, data = pca.scores) #所有主成分，解释69.37%的变异
summary(nhl.lm2)

#简单散点图绘制，nhl.lm2$fitted.values预测值
plot(nhl.lm2$fitted.values, train$ppg, 
     main="Predicted versus Actual",
     xlab="Predicted",ylab="Actual")


train$pred <- round(nhl.lm2$fitted.values, digits = 2)

p <- ggplot(train, aes(x = pred,
                       y = ppg,
                       label = Team)) 
p + geom_point() + 
  geom_text(size=3.5, hjust=0.1, vjust=-0.5, angle=0) + 
  xlim(0.8, 1.4) + ylim(0.8, 1.5) +
  stat_smooth(method="lm", se=FALSE)

pca.scores$Team <- train$Team
p2 <- ggplot(pca.scores, aes(x = RC1, y = RC2, label = Team))
p2 + geom_point() +
  geom_text(size = 2.75, hjust = .2, vjust = -0.75, angle = 0) +
  xlim(-2.5, 2.5) + ylim(-3.0, 2.5)

sqrt(mean(nhl.lm2$residuals^2))

#利用主成分分析结果进行预测
test <- read.csv("NHLtest.csv")
test.scores <- data.frame(predict(pca.rotate, test[, c(-1:-2)]))
test.scores$pred <- predict(nhl.lm2, test.scores)

test.scores$ppg <- test$ppg
test.scores$Team <- test$Team

library(ggplot2)
p <- ggplot(test.scores, aes(x = pred,
                             y = ppg,
                             label = Team)) 
p + geom_point() + 
  geom_text(size=3.5, hjust=0.4, vjust = -0.9, angle = 35) + 
  xlim(0.75, 1.5) + ylim(0.5, 1.6) +
  stat_smooth(method="lm", se=FALSE)


#计算均方误差
resid <- test.scores$ppg - test.scores$pred
sqrt(mean(resid^2))
####
































