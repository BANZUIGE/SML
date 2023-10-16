# 1. 包的载入 ------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(patchwork)

# 2. 数据载入 -----------------------------------------------------------------
IBI <- read.csv("Core_otu_CSS_table_84sample_21core.csv",header = TRUE, row.names = 1,sep = ",",check.names = F)  %>% as.data.frame()
group <- read.csv("../all_data/IBI指数/IBI_group.csv",header = TRUE, row.names = 1,sep = ",",check.names = F)  %>% as.data.frame()


# 3. 数据整理 -----------------------------------------------------------------------
data <- cbind(group[1:2],IBI)
names(data)

data$df <- data$SVM.9396 -data$IBI
data$group <- factor(data$group,levels= c('Test_set','Train_set'),ordered = T)
str(data$group)


max = 4.179795479
data$IBI_EQ <- cut(data$IBI, breaks = c(0, max *0.25,max *0.5,max * 0.75, 4.3), 
                   labels = c("1","2","3","4"), right=FALSE)
data$predict_EQ <- cut(data$SVM.9396,breaks = c(0, max *0.25,max *0.5,max * 0.75, 4.3), 
                       labels = c("1","2","3","4"), right=FALSE)

data$EQ_df <- as.numeric(data$predict_EQ) - as.numeric(data$IBI_EQ) 
head(data)

data2 <- subset(data, group=="Test_set")


# 4. 散点图 ------------------------------------------------------------------
# ??geom_smooth
p <- ggplot(data,aes(IBI,SVM.9396,color=group)) +
  geom_point(size=4)+
  scale_color_manual(values = c("red","gray"))+
  geom_smooth(data=data2,aes(x=IBI, y= SVM.9396),method="lm",se = F)+
  geom_segment(aes(x=0, xend=4.22, y=0, yend=4.22), linetype=2,color= "blue") +

  labs(x="Me-IBI actual value", y = "Me-IBI predicted value") +
  
  annotate("rect", xmin = 3.134847, xmax = 4.22, ymin = 3.134847, ymax = 4.22, alpha = 0,color="#00BFFF") +
  annotate("rect", xmin = 2.089898, xmax = 3.134847, ymin = 2.089898, ymax = 3.134847, alpha = 0,color="#00FA9A")+ 
  annotate("rect", xmin = 1.044949, xmax = 2.089898, ymin = 1.044949, ymax = 2.089898,  alpha = 0,color="#FFA500" )+
  annotate("rect", xmin = 0, xmax = 1.044949, ymin = 0, ymax = 1.044949, alpha = 0,color="red")+

  theme_bw() +
  theme(text=element_text(face="bold",size=20))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        legend.position = "none" 
  )
 
p

ggsave("SVM.9396.a.pdf",width = 5,height = 5,dpi = 300)


# 5. 柱形图 ------------------------------------------------------------------
barplot(table(data2$EQ_df))

adj <- group_by(data2,EQ_df) %>% summarize(count=n()) %>% mutate(perc=((count/28)*100))
adj
p2 <- ggplot(adj, aes(EQ_df, perc)) +
  geom_bar(stat="identity")+
  labs(x=NULL, y = NULL) +
  theme_classic()+
  scale_y_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100)) +
   # scale_x_continuous(limits=c(0,2), breaks=c(0,1)) +
  scale_x_continuous(breaks = seq(-1, 1, by = 1), limits = c(-1.5, 1.5))+
  theme(text=element_text(face="bold",size = 36))


p2

ggsave("SVM.9396.b.pdf",width = 4,height = 4,dpi = 300)


# 6. 箱线图 ------------------------------------------------------------------
data2$x <- rep(1,28)
summary(data2$df)
p3 <- ggplot(data2,aes(x,df))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(width=0.75) +
  xlim(0, 2)+
  labs(x=NULL, y = NULL) +
  stat_summary(fun="mean", geom="point", shape=20, size=10, color="red", fill="red",alpha=1)+
  theme_bw() +
  theme(text=element_text(face="bold",size=36))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.background=element_blank(),
        panel.grid=element_blank(),
        legend.position = "none" 
  )


p3


ggsave("SVM.9396.c.pdf",width = 4,height = 4,dpi = 300)







