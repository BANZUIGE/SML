# 机器学习算法-随机森林初探（1） --------------------------------------------------------
# 读入数据
test <- read.table("DLBCL.expr.txt", row.names = 1, header = T, sep="\t")
expr_mat <- test[1:7070]
metadata <- test[7071]
expr_mat[1:4,1:5]
head(metadata)
group = "class"
metadata[[group]] <- as.factor(metadata[[group]])
class(metadata)
str(metadata)

library(randomForest)
?randomForest

# 直接分析一下，看到结果再调参
# 设置随机数种子，
set.seed(304)

# 直接使用默认参数
rf <- randomForest(expr_mat, metadata$class)
rf
# 增加决策树的数目到1000测试下分类率是否会降低
set.seed(304)
rf1000 <- randomForest(expr_mat, metadata$class, ntree=1000)
rf1000
library(ggplot2)

YSX::sp_lines(as.data.frame(rf1000$err.rate), manual_color_vector="Set2", 
              x_label="Number of trees", y_label="Error rate",line_size=0.6,
              width=6, height=4
)

devtools::install_github("Tong-Chen/YSX")

# 增加树的数目没有给出好的结果，这里还是用的默认的500棵树以便获得较快的运行速度
set.seed(304)
rf_mtry100 <- randomForest(expr_mat, metadata$class, mtry=100)
rf_mtry100

# 一个个测试也不是办法，tuneRF给我们提供了一个根据OOB值迭代鉴定最合适的mtry值的函数。
# (测试几次，不太好用，改一个stepfactor，结果变化很大)

# mtryStart: 从多少个变量开始尝试，可用默认值。程序会自动向更多变量或更少变量迭代。
# ntreeTry: 迭代时构建多少棵树，这里设置为500，与我们上面获得的效果比较好的树一致。
# stepFactor: 迭代步长，mtryStart向更多变量迭代时乘以此值；mtryStart向更少变量迭代时除以此值。
# improve：如果迭代不能给OOB带来给定值的改善，则停止迭代。
set.seed(304)
tuneRF(expr_mat, metadata$class, ntreeTry=500, stepFactor=1.1, improve=1e-5)


# randomForest自带了另一个函数rfcv，通过嵌套交叉验证方式评估了根据变量重要性降低预测变量后的模型的性能。

result = rfcv(expr_mat, metadata$class, cv.fold=10)
result$error.cv




# 假设有一个数据集，包含 6 个样品
m = 6
train_set <- paste0('a', 1:m)
train_set
K = 2
set.seed(1)
# 下面这行代码随机从1:K 中有放回的选取与样品数目一致的索引值
# 从属于相同索引值的样本同属于一个fold
kfold <- sample(1:K, size=m, replace=T)
kfold

table(kfold)



# 机器学习 - 随机森林手动10 折交叉验证 ---------------------------------------------------
library(randomForest)
set.seed(304)
rf1000 <- randomForest(expr_mat, metadata[[group]], ntree=1000)
rf1000

# 查看模型的类，为randomForest
class(rf1000)
?predict
?predict.randomForest

# 开始预测
preds <- predict(rf1000, expr_mat, type="response")
preds

# 计算模型效果评估矩阵（也称混淆矩阵）
caret::confusionMatrix(preds, metadata[[group]])

preds_prob <- predict(rf1000, expr_mat, type="prob")
head(preds_prob)

# predict还可以返回分类的vote值。

preds_prob <- predict(rf1000, expr_mat, type="vote")
head(preds_prob)

# 把前面的代码串起来，就构成了一个随机森林的 10 折交叉验证代码
K = 10
m = nrow(expr_mat)
set.seed(1)
kfold <- sample(rep(1:K, length.out=m), size=m, replace=F)

randomForestCV <- function(x, y, xtest, ytest, type="response", seed=1, ...){
  set.seed(seed)
  model <- randomForest(x, y, ...)
  preds <- predict(model, xtest, type=type)
  return(data.frame(preds, real=ytest))
}

CV_rf <- lapply(1:K, function(x, ...){ 
  train_set = expr_mat[kfold != x,]
  train_label = metadata[[group]][kfold!=x]
  
  validate_set = expr_mat[kfold == x,]
  validate_label = metadata[[group]][kfold==x]
  
  randomForestCV(x=train_set, y=train_label, xtest=validate_set, ytest=validate_label, ...)
})
CV_rf[1]

kfold_estimate <- do.call(rbind, CV_rf)


# 基于Caret和RandomForest包进行随机森林分析的一般步骤 （1） ----------------------------------
# 拆分数据为测试集和训练集
seed <- 1
set.seed(seed)
train_index <- createDataPartition(metadata[[group]], p=0.75, list=F)
train_data <- expr_mat[train_index,]
train_data_group <- metadata[[group]][train_index]

test_data <- expr_mat[-train_index,]
test_data_group <- metadata[[group]][-train_index]

# 借助Caret快速构建一个初步模型：默认参数

# Create model with default parameters，设置10折交叉验证
?trainControl
trControl <- trainControl(method="repeatedcv", number=10, repeats=3)
# train model<span style="color: rgb(51, 51, 51);font-weight: bold;">if(file.exists('rda/rf_default.rda')){
  # rf_default <- readRDS("rda/rf_default.rda")
# } else {  # 设置随机数种子，使得结果可重复
  seed <- 1
  set.seed(seed)
  
  rf_default <- train(x=train_data, y=train_data_group, method="rf", 
                      trControl=trControl)
  rf_default
  rf_default$bestTune
  
  
  #将结果保存为rda文件
  saveRDS(rf_default, "rf_default.rda")
  ls <- predict(rf_default,test_data)
ls

# Caret训练最终模型
# method="none" 不再抽样，用全部数据训练模型
trControl <- trainControl(method="none", classProbs = T)  # 设置随机数种子，使得结果可重复
seed <- 1
set.seed(seed)

rf_final <- train(x=train_data, y=train_data_group, method="rf", 
                  # 用最好的模型的调参值
                  tuneGrid = rf_default$bestTune,
                  trControl=trControl)
rf_final 
saveRDS(rf_final, "rf_finaltest.rda")
readRDS()

# 基于模型对测试集进行预测
# ?predict.train 
# type: raw (返回预测的类或回归的数值) or prob （返回分类到每个类的概率） 
predict(rf_final, newdata=head(test_data))

?train
svm_default <- train(x=train_data, y=train_data_group, method="svmLinear")
svm_default$bestTune




# R包randomForest的随机森林回归模型以及对重要变量的选择 ---------------------------------------
##数据预处理
#读取 OTUs 丰度表
otu <- read.delim('otu_table.txt', row.names = 1)

#过滤低丰度 OTUs 类群，它们对分类贡献度低，且影响计算效率
#例如剔除总丰度低于 0.05% 的值
otu <- otu[which(rowSums(otu) >= 0.0005), ]

#合并有关于植物生长时间的信息
plant <- read.delim('plant_age.txt', row.names = 1)

otu <- data.frame(t(otu))
otu <- otu[rownames(plant), ]
otu <- cbind(otu, plant)

#为了方便后续评估随机森林模型的性能
#将总数据集分为训练集（占 70%）和测试集（占 30%）
set.seed(123)
train <- sample(nrow(otu), nrow(otu)*0.7)
otu_train <- otu[train, ]
otu_test <- otu[-train, ]

##randomForest 包的随机森林
library(randomForest)

#随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
set.seed(123)
otu_train.forest <- randomForest(plant_age~., data = otu_train, importance = TRUE)
otu_train.forest

#使用训练集，查看预测精度
plant_predict <- predict(otu_train.forest, otu_train)

plot(otu_train$plant_age, plant_predict, main = '训练集', 
     xlab = 'Plant age (days)', ylab = 'Predict')
abline(1, 1)

#使用测试集，评估预测性能
plant_predict <- predict(otu_train.forest, otu_test)

plot(otu_test$plant_age, plant_predict, main = '测试集',
     xlab = 'Plant age (days)', ylab = 'Predict')
abline(1, 1)

##OTU 的重要性评估
#查看表示每个预测变量（细菌 OTU）重要性的得分
#summary(otu_train.forest)
importance_otu <- otu_train.forest$importance
head(importance_otu)

#或者使用函数 importance()
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)

#作图展示 top30 重要的 OTUs
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)), 
           main = 'Top 30 - variable importance')

#可以根据某种重要性的高低排个序，例如根据“IncNodePurity”指标
importance_otu <- importance_otu[order(importance_otu$IncNodePurity, decreasing = TRUE), ]
head(importance_otu)

#输出表格
#write.table(importance_otu, 'importance_otu.txt', sep = '\t', col.names = NA, quote = FALSE)

##交叉验证辅助评估选择特定数量的 OTU
#5 次重复十折交叉验证
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu_train[-ncol(otu_train)], otu_train$plant_age, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv[1]

#提取验证结果绘图
otu_train.cv <- data.frame(sapply(otu_train.cv, '[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 10)

#拟合线图
library(ggplot2)

ggplot(otu_train.cv.mean, aes(Group.1, x)) +
  geom_line() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +  
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')

#提示保留 9-13 个重要的 OTU，可以使随机森林回归的精度最大化
#首先根据某种重要性的高低排个序，例如根据“IncNodePurity”指标
importance_otu <- importance_otu[order(importance_otu$IncNodePurity, decreasing = TRUE), ]

#然后取出排名靠前的 OTU，例如 top10 最重要的 OTU
importance_otu.select <- importance_otu[1:10, ]
importance_otu.select

#输出表格
#write.table(importance_otu.select, 'importance_otu.select.txt', sep = '\t', col.names = NA, quote = FALSE)

#不妨简单查看下这些重要的 OTU 丰度与植物生长时间的关系
#可以看到趋势非常明显，包括根际富集或排斥等都有涉及
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu[ ,c(otu_id.select, 'plant_age')]
otu.select <- reshape2::melt(otu.select, id = 'plant_age')

ggplot(otu.select, aes(x = plant_age, y = value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~variable, ncol = 3, scale = 'free_y') +
  labs(title = '',x = 'Plant age (days)', y = 'Relative abundance')



##只包含 10 个重要预测变量的简约回归
otu.select <- otu[ ,c(otu_id.select, 'plant_age')]

#为了方便后续评估随机森林模型的性能，将总数据集分为训练集（占 70%）和测试集（占 30%）
set.seed(123)
train <- sample(nrow(otu.select), nrow(otu.select)*0.7)
otu_train.select <- otu.select[train, ]
otu_test.select <- otu.select[-train, ]

#随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
set.seed(123)
otu_train.select.forest <- randomForest(plant_age~., data = otu_train.select, importance = TRUE)
otu_train.select.forest

#使用训练集，查看预测精度
plant_predict <- predict(otu_train.select.forest, otu_train.select)

plot(otu_train.select$plant_age, plant_predict, main = '训练集', 
     xlab = 'Plant age (days)', ylab = 'Predict')
abline(1, 1)

#使用测试集，评估预测性能
plant_predict <- predict(otu_train.select.forest, otu_test.select)

plot(otu_test.select$plant_age, plant_predict, main = '测试集',
     xlab = 'Plant age (days)', ylab = 'Predict')
abline(1, 1)

































