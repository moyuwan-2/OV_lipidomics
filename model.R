####single_lipid-based_ROC_model######
library(caret)
library(e1071)
library(pROC)
library(dplyr)
library(purrr)
#
lipid_data <-  readxl::read_xlsx("./lipid_data.xlsx") %>% subset(select!=0) %>% subset(dif_change!="D") 
sample_pheno <- read.csv("./sample_pheno.csv")


names <- c("EOC (n = 24) vs. BOT(n = 22)",
           "EOC (n = 24) vs. HC (n = 26)", 
           "BOT (n = 22) vs. HC (n = 26)")

TOP_roc <- function(groups,vsnames,change_name){
  
  f3 <- lipid_data%>% subset(get(change_name)!="NS")
  df <- f3 %>% select(c(1,24:131)) %>% t() %>% as.data.frame() %>% setNames(.[1,]) %>% .[-1,] %>% mutate(name=rownames(.)) %>% select(name,everything())
  my_colros <- c("#FF4040","#00C5CD","#4876FF","#32CD32","#EE9A00")
  
  
  
  pheno <- sample_pheno %>% subset(set%in%"train") %>% subset(group%in%groups) %>% select(all_of(c("sample","group")) )
  
  df_tr <- pheno %>% left_join(df,by=c("sample"="name"))%>% mutate_at(vars(-c(1, 2)), as.numeric)
  
  name_to_lipid <- setNames(lipid_data$`Metabolite name`,lipid_data$name)
  ROC <- function(df,var){
    roc_analysis <- roc(response = df$group,predictor = as.numeric(df[,var]))
    return(roc_analysis)
  }
  var_names <-colnames(df_tr)[-c(1:2)] 
  roc_re <-setNames(map(var_names ,~ROC(df_tr,.x)),var_names ) 
  AUC_val <- lapply(roc_re,auc)
  AUC_val1 <-Reduce(rbind,map2(AUC_val,names(AUC_val),~cbind(.x,name=.y))) %>% as.data.frame()%>% slice_max(.x,n=5) #TOP5AUC
  TOP_AUC <- AUC_val1$name
  roc_re_top <- setNames(lapply(TOP_AUC,function(x) roc_re[[x]]),TOP_AUC)
  
  proc <-  ggroc(roc_re_top,alpha=1,legacy.axes=T)+
    scale_color_manual(values =  my_colros)+
    ggtitle(paste(vsnames))+
    theme_bw()+
    #annotate(geom="segment",x=0,y=0,xend=1,yend=1,linetype=2)+
    theme(legend.title = element_blank() ,plot.background = element_blank(),legend.position ="",  panel.grid=element_blank(),plot.title = element_text(hjust = 0.5,size = 10))
  
  lipid_name <- setNames(lapply(TOP_AUC, function(x) name_to_lipid[[x]]),TOP_AUC)
  auc_values <-  setNames(lapply(TOP_AUC, function(x) AUC_val[[x]]),TOP_AUC)
  for(i in 1: 5){
    proc <- proc  +
      annotate("segment", x = 0.28, xend = 0.33, y = 0.09 * i, yend = 0.09 * i,colour = my_colros[i])+
      annotate("text", x = 0.35, y = 0.09 * i, label = lipid_name[[i]],hjust=0,size=2.2)+
      annotate("text", x = 0.75, y = 0.09 * i, label =  paste("AUC =", sprintf("%.3f", auc_values[[i]])),hjust=0,size=2.2)
  }
  
  return(proc)
}
ROC_EOC_HC <- TOP_roc(groups=c("EOC","HC"),vsnames=names[2],change_name="change_EOC_HC")
ROC_BOT_HC <- TOP_roc(groups=c("BOT","HC"),vsnames=names[3],change_name="change_BOT_HC")
ROC_EOC_BOT <- TOP_roc(groups=c("EOC","BOT"),vsnames=names[1],change_name="change_EOC_BOT")
###TOP_AUC_performence



TOP_roc_acc <- function(groups,change_name){
  
  f3 <- lipid_data%>% subset(get(change_name)!="NS")
  df <- f3 %>% select(c(1,24:131)) %>% t() %>% as.data.frame() %>% setNames(.[1,]) %>% .[-1,] %>% mutate(name=rownames(.)) %>% select(name,everything())
  
  pheno_tr <- sample_pheno %>% subset(set%in%"train") %>% subset(group%in%groups) %>% select(all_of(c("sample","group")) )
  
  pheno_te <- sample_pheno %>% subset(set%in%"test") %>% subset(group%in%groups) %>% select(all_of(c("sample","group")) )
  
  df_tr <- pheno_tr %>% left_join(df,by=c("sample"="name")) %>% mutate_at(vars(-c(1, 2)), as.numeric)
  df_te <- pheno_te %>% left_join(df,by=c("sample"="name")) %>% mutate_at(vars(-c(1, 2)), as.numeric)
  
  
  
  ROC <- function(df,var){
    
    roc_analysis <- roc(response = df$group,predictor = as.numeric(df[,var]))
    return(roc_analysis)
  }
  
  
  var_names <-colnames(df_tr)[-c(1:2)] 
  roc_re <-setNames(map(var_names ,~ROC(df_tr,.x)),var_names ) 
  
  AUC_val <- lapply(roc_re,auc)
  AUC_val1 <-Reduce(rbind,map2(AUC_val,names(AUC_val),~cbind(.x,name=.y))) %>% as.data.frame()%>% slice_max(.x,n=5) #TOP5AUC
  TOP_AUC <- AUC_val1$name
  roc_re_top <- setNames(lapply(TOP_AUC,function(x) roc_re[[x]]),TOP_AUC)
  
  threshold <- lapply(roc_re_top, function(x){roc_result <- coords(x,"best",best.method=c("youden"),ret=c("threshold"))})
  
  threshold1 <- Reduce(rbind,map2(threshold,names(threshold),~data.frame(.x,name=.y))) %>% mutate(change=map_chr(
    name,~lipid_data%>% filter(name%in%.x) %>% pull(change_name) %>% as.character()
  ))
  
  df_te_1 <- df_te %>% 
    select(all_of(c("sample","group",threshold1$name)))
  
  result <- list()
  Accuracy <- list()
  Accuracy_low <- list()
  Accuracy_up <- list()
  Accuracy_p <- list()
  Sensitivity <- list()
  Specificity <- list()
  acc_result <- list()
  for(i in  threshold1$name){
    
    dir <- threshold1 %>% filter(name%in%i) %>% pull(change) %>% as.character()
    if(dir=="Up"){
      threshold_val <- threshold1 %>% filter(name%in%i) %>% pull(threshold) %>% as.numeric()
      df <- df_te_1 %>% select(all_of(c("group",i) ))
      df$pre <- ifelse(df[[i]]>=threshold_val,groups[1],groups[2])
      
      result[[i]] <- confusionMatrix(as.factor(df$group),as.factor(df$pre),positive = groups[1])#ture,predict
      
      
      
    }else{
      threshold_val <- threshold1 %>% filter(name%in%i) %>% pull(threshold) %>% as.numeric()
      df <- df_te_1 %>% select(all_of(c("group",i) ))
      df$pre <- ifelse(df[[i]]<=threshold_val,groups[1],groups[2])
      
      result[[i]] <- confusionMatrix(as.factor(df$group),as.factor(df$pre),positive = groups[1])
      
    }
    
    
    
    
    Accuracy[[i]] <- result[[i]][["overall"]][["Accuracy"]]
    Accuracy_low[[i]]<- result[[i]][["overall"]][["AccuracyLower"]]
    Accuracy_up[[i]]<- result[[i]][["overall"]][["AccuracyUpper"]]
    Accuracy_p[[i]]<- result[[i]][["overall"]][["AccuracyPValue"]]
    Sensitivity[[i]]<- result[[i]][["byClass"]][["Sensitivity"]]
    Specificity[[i]]<- result[[i]][["byClass"]][["Specificity"]]
    acc_result[[i]] <- c(Accuracy[[i]],Accuracy_low[[i]],Accuracy_up[[i]],
                         Accuracy_p[[i]], Sensitivity[[i]],Specificity[[i]])
    
  }
  acc_result_1 <- Reduce(rbind,acc_result) %>% as.data.frame() %>% mutate(name=names(acc_result))
  colnames(acc_result_1) <- c("Accuracy","Accuracy_low",
                              "Accuracy_up"," Accuracy_p",
                              "Sensitivity","Specificity","name")
  threshold_result <- left_join(threshold1,acc_result_1,by="name")
  return(threshold_result)
}
ROC_EOC_HC_acc <- TOP_roc_acc(groups=c("EOC","HC"),change_name="change_EOC_HC") %>% mutate(class="EOC vs HC")
ROC_BOT_HC_acc <- TOP_roc_acc(groups=c("BOT","HC"),change_name="change_BOT_HC") %>% mutate(class="BOT vs HC")
ROC_EOC_BOT_acc <- TOP_roc_acc(groups=c("EOC","BOT"),change_name="change_EOC_BOT") %>% mutate(class="EOC vs BOT")

ROC_all_acc <- rbind(ROC_EOC_HC_acc,ROC_BOT_HC_acc ,ROC_EOC_BOT_acc)
name_to_lipid <- setNames(lipid_data$`Metabolite name`,lipid_data$name)
ROC_all_acc$lipid <- name_to_lipid [ROC_all_acc$name]
dir.create("./Single_lipid/")
write.csv(ROC_all_acc,"./Single_lipid/top5_lipid_performance.csv")
######multiple_lipid-based_ML_model####

lipid_data <- readxl::read_xlsx("./lipid_data.xlsx") %>% subset(select!=0) %>% subset(dif_change!="D") 

#
df_raw <-read.csv("./df_raw.csv")

trainset <- subset(df_raw,set%in%"train") %>% select(all_of(c(colnames(df_raw)[1:6],lipid_data$name))) 
testset <- subset(df_raw,set%in%"test") %>% select(all_of(c(colnames(df_raw)[1:6],lipid_data$name))) 
#
#Binary classification modeling
fitControl <- trainControl( 
  method = "repeatedcv",
  number = 3,
  repeats = 5)
#
levels <- c("EOC","BOT","HC")
exgroup <- c("HC","BOT","EOC")
#
trainset$group <- factor(trainset$group,levels = levels)
testset$group <- factor(testset$group,levels = levels)
trainset1  <- trainset [,-c(3:6)]
testset1  <- testset[,-c(3:6)]
#
ml_result <- list()
result <- list()
#knn
fit.knn <- list()
te_knn_result <- list()
#pls
fit.pls <- list()
te_pls_result <- list()
#rf
fit.rf <- list()
te_rf_result <- list()
#svm_linear
fit.svm.linear <- list()
te_svm_linear_result <- list()
trainset_snv_0 <- list()
testset_snv_0 <- list()

model <- c("knn","pls","rf","svm_linear")
#

RESULT <- function(data,i){
  setnames <-data.frame(exgroup=c("EOC","BOT","HC"),
                        name=c("BOT_HC","EOC_HC","EOC_BOT")) 
  setname <- setNames(setnames$name,setnames$exgroup)
  name <- setname[[i]]
  te<- get(paste0("te_",data,"_result"))#字符串变为变量
  statistics <- c("Accuracy","AccuracyLower","AccuracyUpper","Sensitivity",
                  "Specificity","F1")
  acc <- c("Accuracy","AccuracyLower","AccuracyUpper")
  result$exclude <- i
  result$group <- name
  
  for (j in statistics){
    
    if( j %in% acc){
      type <- "overall"
    } else{
      type <- "byClass"
    }
    result[[paste0("te_",j)]] <- te[[i]][[type]][[j]]
    result[["ML"]] <- data
    
  }
  
  
  
  result<- as.data.frame(result)
  
  return(result)
}

for (k in model){
  for(i in exgroup ){
    trainset2 <- subset(trainset1,group!=i)
    testset2 <- subset(testset1,group!=i)
    
    level <- levels(droplevels(trainset2$group))
    
    trainset2$group <- factor(trainset2$group, levels =level)
    testset2$group<- factor(testset2$group, levels = level)
    
    prePro <- preProcess(trainset2[,-c(1:2)],
                         method = c("center", "scale"))
    train_snv <- predict(prePro,trainset2[,-c(1:2)])
    test_snv <- predict(prePro,testset2[,-c(1:2)])
    
    trainset_snv_0[[i]] <- cbind(trainset2[,c(1:2)],train_snv)
    testset_snv_0[[i]] <- cbind(testset2[,c(1:2)],test_snv)
    
    trainset_snv <-  trainset_snv_0[[i]]
    testset_snv <- testset_snv_0[[i]] 
    
    ####class_
    #KNN
    set.seed(38)
    system.time(fit.knn[[i]]<- train(x = trainset_snv[,-c(1:2)], 
                                     y = trainset_snv$group,
                                     method = "knn",
                                     trControl = fitControl,
                                     tuneLength = 10))
    
    
    te_knn_result[[i]]<- confusionMatrix(predict(fit.knn[[i]],testset_snv),as.factor(testset_snv$group))
    ##pls
    set.seed(38)
    system.time(
      fit.pls[[i]]<- train(y = trainset_snv$group ,
                           x = trainset_snv[,-c(1:2)], 
                           method = "pls",
                           trControl = fitControl,
                           verbose = FALSE
      ))
    
    
    te_pls_result[[i]] <- confusionMatrix(predict(fit.pls[[i]] ,testset_snv),as.factor(testset_snv$group))
    #RF
    set.seed(38)
    system.time(
      fit.rf[[i]]<- train(y = trainset_snv$group, 
                          x = trainset_snv[,-c(1:2)], 
                          method = "rf",   
                          trControl = fitControl,
                          verbose = FALSE
      ))
    
    te_rf_result[[i]] <-confusionMatrix(predict(fit.rf[[i]],testset_snv),as.factor(testset_snv$group))
    ###svm_linear
    set.seed(38)
    
    fit.svm.linear[[i]] <- svm( y = as.factor(trainset_snv$group), 
                                x = trainset_snv[,-c(1:2)],
                                cross = 10,
                                kernel = "linear",
                                type = "C-classification",
                                scale =F,
                                probability= TRUE,
                                cost =1)#     from 0.001 to 100
    te_svm_linear_result[[i]] <- confusionMatrix(predict(fit.svm.linear[[i]],testset_snv[,-c(1:2)]),as.factor(testset_snv$group)); te_svm_linear_result[[i]] 
    result <- RESULT(data=k,i=i)
    ml_result <- rbind(ml_result,result) 
  }
}
dir.create("./multiple_lipid/")
write.csv(ml_result,"./multiple_lipid/ml_performance.csv")
#

names <- c("EOC","EOC","BOT")
vsnames <- c("EOC (n = 24) vs. BOT(n = 22)","EOC (n = 24) vs. HC (n = 26)", "BOT (n = 22) vs. HC (n = 26)")
model_names <- c("K-NN","PLS","SVM","RF")
my_colros <- c("K-NN"="#BF3EFF","PLS"="#76EE00",
               "SVM"="#1E90FF","RF"="#FF3E96")
knn <- list()
pls <- list()
svm <- list()
rf <- list()
rknn <- list()
rpls <- list()
rsvm <- list()
rrf <- list()
roc_ <- list()
auc_values <- list()
roc_result <- list()
#ML_roc
for(i in seq_along(exgroup) ){
  ex <- exgroup[i]
  name  <- names [i]
  vsname <- vsnames[i]
  
  
  knn[[ex]] <- predict(fit.knn[[ex]],trainset_snv_0[[ex]],type = "prob")
  knn[[ex]]$group<- trainset_snv_0[[ex]][["group"]]
  
  
  #
  pls[[ex]] <- predict(fit.pls[[ex]],trainset_snv_0[[ex]],type = "prob")
  pls[[ex]]$group<- trainset_snv_0[[ex]][["group"]]
  
  #
  svm[[ex]] <- predict(fit.svm.linear[[ex]],trainset_snv_0[[ex]][,-c(1:2)],probability = TRUE)
  svm[[ex]] <- as.data.frame(attr(svm[[ex]], "probabilities"))
  svm[[ex]]$group <- trainset_snv_0[[ex]]$group
  
  rf[[ex]] <- predict(fit.rf[[ex]],trainset_snv_0[[ex]],type = "prob")
  rf[[ex]]$group<- trainset_snv_0[[ex]][["class2"]]
  
  
  knn[[ex]]$group <- factor(knn[[ex]]$group,levels = exgroup)
  level <- knn[[ex]] %>% .$group %>% droplevels() %>% levels()
  
  rknn[[ex]][[name]] <- roc(trainset_snv_0[[ex]]$group,as.numeric(unlist(knn[[ex]][[name]])),levels=level,direction="<")
  rpls[[ex]][[name]] <- roc(trainset_snv_0[[ex]]$group,as.numeric(unlist(pls[[ex]][[name]])),levels=level,direction="<")
  rsvm[[ex]][[name]] <- roc(trainset_snv_0[[ex]]$group,as.numeric(unlist(svm[[ex]][[name]])),levels=level,direction="<")
  rrf[[ex]][[name]] <- roc(trainset_snv_0[[ex]]$group,as.numeric(unlist(rf[[ex]][[name]])),levels=level,direction="<")
  
  
  roc_[[ex]][[name]]<- ggroc(list("K-NN"=rknn[[ex]][[name]],"PLS"=rpls[[ex]][[name]],
                                  "SVM"=rsvm[[ex]][[name]],"RF"=rrf[[ex]][[name]] ),alpha=1,legacy.axes=T)+
    theme_bw()+
    ggtitle(paste(vsnames[i]))+
    
    theme(legend.title = element_blank() ,plot.background = element_blank(),legend.position =" ",  panel.grid=element_blank(),plot.title = element_text(size = 10,hjust =0.5 ))+
    scale_color_manual(values = my_colros )
  
  
  
  auc_values[[ex]][[name]] <- c(auc(rknn[[ex]][[name]]),
                                auc(rpls[[ex]][[name]]),
                                auc(rsvm[[ex]][[name]]),
                                auc(rrf[[ex]][[name]]))
  
  for (j in 1:length(model_names)) {
    
    roc_[[ex]][[name]] <- roc_[[ex]][[name]] + 
      annotate("segment", x = 0.28, xend = 0.33, y = 0.09 * j, yend = 0.09 * j,colour = my_colros[model_names[j]])+
      annotate("text", x = 0.35, y = 0.09 * j, label = model_names[j],hjust=0,size=2.2)+
      annotate("text", x = 0.75, y = 0.09 * j, label =  paste("AUC =", sprintf("%.3f",auc_values[[ex]][[name]][j])),hjust=0,size=2.2)
    
    
    
  }
} 


#important_lipids_in_ML#
fit_TOP10 <- function(fit){
  var <- varImp(fit,scale = T)
  var_1  <-   var$importance
  var_2<-   var_1 %>% slice_max(.[,1],n=10)
}
knn_top10 <- lapply(fit.knn, fit_TOP10)
pls_top10 <- lapply(fit.pls, fit_TOP10)
rf_top10 <- lapply(fit.rf, fit_TOP10)
#svm_RFE
SVM_rank <- function(exgroup,level){
  
  svmrfeFeatureRanking = function(x,y){
    n = ncol(x)
    survivingFeaturesIndexes = seq(1:n)
    featureRankedList = vector(length=n)
    rankedFeatureIndex = n
    while(length(survivingFeaturesIndexes)>0){
     
      svmModel <- svm(x[, survivingFeaturesIndexes], 
                      y, 
                      cost = 1, 
                      cachesize=500,  
                      scale= T, 
                      type="C-classification", 
                      kernel="linear" )
      
      w <- t(svmModel$coefs) %*% svmModel$SV 
     
      rankingCriteria <- apply(w,2,
                               function(v){
                                 sqrt(sum(v^2))
                               })
      #rank the features from minimum to maximum.
      ranking <- sort(rankingCriteria, index.return = TRUE)$ix
      #update feature ranked list 
      featureRankedList[rankedFeatureIndex] <- survivingFeaturesIndexes[ranking[1]]
      rankedFeatureIndex <- rankedFeatureIndex - 1
      #eliminate the feature with smallest ranking criterion
      (survivingFeaturesIndexes <- survivingFeaturesIndexes[-ranking[1]])
    }
    return (featureRankedList) 
  }
  
  
  
  trainset2 <- subset(trainset1,group!=exgroup)
  
  
  level <- levels(droplevels(trainset2$group))
  
  trainset2$group <- factor(trainset2$group, levels =level)
  
  
  prePro <- preProcess(trainset2[,-c(1:2)],
                       method = c("center", "scale"))
  trainset_snv0 <- predict(prePro,trainset2[,-c(1:2)])
  trainset_snv <- cbind(trainset2[,c(1:2)],trainset_snv0 )
  
  
  rankedsvm1 <- svmrfeFeatureRanking(trainset_snv[,-c(1:2)],trainset_snv$group)   
  
  col <- colnames(trainset_snv[,-c(1:2)][c(rankedsvm1[1:317])])
  svm_ranking_fearture <- data.frame(col) %>% mutate(ranking=rownames(.))
  return(svm_ranking_fearture)
}
SVM_ranking_EOC_HC<- SVM_rank("BOT",c("EOC","HC")) %>% subset(ranking%in%c(1:317))
SVM_ranking_BOT_HC <- SVM_rank("EOC",c("BOT","HC"))%>% subset(ranking%in%c(1:317))
SVM_ranking_EOC_HC <- SVM_rank("HC",level=c("EOC","BOT"))%>% subset(ranking%in%c(1:317))
