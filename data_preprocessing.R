####data_preprocessing####
library(statTarget)
library(caret)
library(dplyr)
setwd()
root <- getwd()

path_samplelist <- paste0(root,"/samplelist.csv")

path_xcmsdata_pos <- paste0(root,"/afterxcms_pos.csv")

dir.create(paste0(root,"/path_statTarget_pos"))
setwd(paste0(root,"/path_statTarget_pos"))
shiftCor(path_samplelist,path_xcmsdata_pos, Frule = 0.8,MLmethod = "QCRFSC",ntree = 500, imputeM = "KNN",coCV = 30,  plot = FALSE)



path_xcmsdata_neg <- paste0(root,"./afterxcms_neg.csv")

dir.create(paste0(root,"/path_statTarget_neg"))
setwd(paste0(root,"/path_statTarget_neg"))
shiftCor(path_samplelist,path_xcmsdata_neg, Frule = 0.8,MLmethod = "QCRFSC",ntree = 500, imputeM = "KNN",coCV = 30,  plot = FALSE)
#
setwd(root)

df_pos <- read.csv(paste0(root,"/path_statTarget_pos/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv"))
df_neg <- read.csv(paste0(root,"/path_statTarget_neg/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv"))

QC_pos_data <- read.csv(paste0(root,"/path_statTarget_pos/statTarget/shiftCor/Before_shiftCor/shift_QC_raw.csv"))
QC_neg_data <- read.csv(paste0(root,"/path_statTarget_neg/statTarget/shiftCor/Before_shiftCor/shift_QC_raw.csv"))

QC_mean_pos  <- mean(colSums(QC_pos_data[,-c(1:3)]))
QC_mean_neg  <- mean(colSums(QC_neg_data[,-c(1:3)]))

df_pos_norm <- cbind(df_pos[,c(1,2)],log2(QC_mean_pos*prop.table(as.matrix(df_pos[,-c(1,2)]),margin = 1)))
df_neg_norm <- cbind(df_neg[,c(1,2)],log2(QC_mean_neg*prop.table(as.matrix(df_neg[,-c(1,2)]),margin = 1)))

#
setwd(root)
df_pos<- df_pos_norm  %>% 
  mutate(sample = toupper(sample))

df_neg <- df_neg_norm  %>%
  mutate(sample = toupper(sample))

class <- read.csv("./sample_group.csv")

#
df_pos1 <- merge(class,df_pos,sort = F,all.y = T)
df_pos1$group <- ifelse(df_pos1$group%in%NA,"QC",df_pos1$group)
df_pos2 <- subset(df_pos1,group!="QC")
df_pos2<- df_pos2 %>% 
  mutate(class = factor(class, levels = c("B", "A"))) 
#
df_neg1 <- merge(class,df_neg,all.y = T) 
df_neg1$group <- ifelse(df_neg1$group%in%NA,"QC",df_neg1$group)
df_neg2 <- subset(df_neg1,group!="QC")
df_neg2<- df_neg2 %>% 
  mutate(class = factor(class, levels = c("B", "A")))
#
set.seed(10)

train <- createDataPartition(df_pos2$group,p=2/3,list=F)

trainset_sample <- df_pos2[train,] %>% select(sample) 
testset_sample <-df_pos2[-train,]%>% select(sample) 

#
df_pos1 <- df_pos1 %>% mutate(group=case_when(group=="Benign"~"BOT",
                                              group=="Malignant"~"EOC",
                                              T~group))
df_train_pos <- df_pos1%>% subset(!sample%in%testset_sample$sample)%>% select(-"class")
df_test_pos <- df_pos1 %>% subset(sample%in%testset_sample$sample ) %>% select(-"class")
write.csv(df_train_pos ,"./df_trainset_pos.csv",row.names = F)
write.csv(df_test_pos,"./df_testset_test_pos.csv",row.names = F)
#
df_neg1 <- df_neg1%>% mutate(group=case_when(group=="Benign"~"BOT",
                                              group=="Malignant"~"EOC",
                                              T~group))
df_train_neg <- df_neg1%>% subset(!sample%in%testset_sample$sample)%>% select(-"class")
df_test_neg <- df_neg1 %>% subset(sample%in%testset_sample$sample ) %>% select(-"class")
write.csv(df_train_neg,"./df_trainset_neg.csv",row.names = F)
write.csv(df_test_neg ,"./df_testset_test_neg.csv",row.names = F)
