# 1. 包的载入 ------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(patchwork)

# 2. 数据载入 -----------------------------------------------------------------
IBI <- read.csv("Core_otu_CSS_table_84sample_21core.csv",header = TRUE, row.names = 1,sep = ",",check.names = F)  %>% as.data.frame()
group <- read.csv("../all_data/IBI指数/IBI_group.csv",header = TRUE, row.names = 1,sep = ",",check.names = F)  %>% as.data.frame()


# 4. 数据整理 -----------------------------------------------------------------------
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

# ??geom_smooth
# 5. 散点图 ------------------------------------------------------------------
p <- ggplot(data,aes(IBI,SVM.9396,color=group)) +
  geom_point(size=4)+
  scale_color_manual(values = c("red","gray"))+
  geom_smooth(data=data2,aes(x=IBI, y= SVM.9396),method="lm",se = F)+
  geom_segment(aes(x=0, xend=4.22, y=0, yend=4.22), linetype=2,color= "blue") +

  labs(x="Mt-IBI actual value", y = "Mt-IBI predicted value") +
  
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


# 6. 柱形图 ------------------------------------------------------------------
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


# 7. 箱线图 ------------------------------------------------------------------
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


# 8. 核心Feature密度图 ---------------------------------------------------------



otu.select <- reshape2::melt(IBI, id = 'IBI')

p <- ggplot(otu.select, aes(x = IBI, y = value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~variable, ncol = 4, scale = 'free_y') +
  labs(title = '',x = 'IBI', y = 'CSS_normalized')
p 

ggsave("Core.pdf",width = 10,height = 10,dpi = 300)

 
  max = 4.179795479
  IBI$ESG <- cut(IBI$IBI, breaks = c(0, max *0.25,max *0.5,max * 0.75, 4.3), 
                   labels = c("1","2","3","4"), right=FALSE)

  
df <- as.data.frame(IBI) %>% 
rownames_to_column("sample") 

names(df)

df_filt <- filter(df, OTU1192 > 0)

p1 <- ggplot(df_filt, aes(x=IBI, y=OTU1192, col=as.factor(ESG)), inherit.aes=F) +
  geom_point(size=3) +
  theme_classic()+
  scale_x_continuous(limits=c(0,4.22), breaks=c(1,2,3,4)) +
  labs(x="Mt-IBI", y="CSS normalized abundance") +
  ggtitle("OTU1192") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  scale_color_manual(values=c("red","#FFA500","#00FA9A","#00BFFF" ))


p2 <- ggplot(df_filt, aes(IBI)) + 
  geom_density(aes(fill=OTU1192), adjust = 2 ) +
  scale_x_continuous(limits=c(0,4.22), breaks=c(1,2,3,4))+
  # scale_y_continuous(limits=c(0,4.22), breaks=c(0.25,0.5,0.75,1))+
  theme_classic()

p <- p1+p2

ggsave("OTU1192.pdf",width = 6,height = 3)








p <- p1+p2



p


ggsave("test.pdf",width = 8,height = 4)


??geom_density
p <- ggplot(ciliates_df_filt_OTU8366, aes(x=IBI, y=OTU8366), inherit.aes=F) +
  geom_point(size=3) +
  geom_density(aes(OTU8366,y= ..scaled..), 
               fill='transparent',cex=1)




library(patchwork)
plot_points_OTU8366
plot_curve_OTU8366
??inset_element
plot_points_OTU8366 + inset_element(plot_curve_OTU8366,left = 1.5,right = 2,bottom = 1.5,top = 2,clip = F) 
p3 <- ggplotGrob(p2)
p4 <- p1+annotation_custom(p3,xmin = 0.1,xmax = 3,ymin = 0.1,ymax = 3)


grid_points <- ggarrange(plot_points_OTU8366,plot_points_OTU8366,ncol=2, nrow = 1)

grid_curve <- ggarrange(plot_curve_OTU8366,plot_points_OTU8366,ncol=2, nrow = 1)

library(tidyverse)
library(reshape2)
library(patchwork)

set.seed(135)
df <- data.frame(GFP=c(rnorm(60,5,1.5),rnorm(90,12,2)),
                 MSR=c(rnorm(110,5,1.5),rnorm(40,12,2.5)))
data <- melt(df)

ggplot(data,aes(x=value,fill=variable))+
  geom_histogram(binwidth=1,
                 color='black')+ #使用aes(y=after_stat(count/sum(count)))则可以绘制频率分布
  scale_fill_manual(values = c('#03B0AB','#FB632E'))

ggplot(data,aes(x=value,fill=variable))+
  geom_histogram(data = filter(data,variable=='GFP'),
                 binwidth=1,alpha=0.7,color='black')+
  geom_histogram(data = filter(data,variable=='MSR'),
                 binwidth=1,alpha=0.7,color='black')+ 
  scale_fill_manual(values = c('#03B0AB','#FB632E'))+
  #加上密度曲线，默认是频率密度，加上y=..count..则为频数密度
  geom_density(aes(value,y=..count..,color=variable), 
               fill='transparent',cex=1)+
  scale_color_manual(values = c('#03B0AB','#FB632E'))






