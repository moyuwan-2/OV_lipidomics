#####multivariate_statistics####
library(dplyr)
library(ropls)
library(plot3D)
library(ggplot2)
#####PCA
#pos


df_train_pos <- read.csv("./df_trainset_pos.csv")

pca <- opls(df_train_pos[,-c(1:2)],crossvalI = 88,predI = 3)
df_pca<- as.data.frame(pca@scoreMN)
df_pca$class0 <- df_train_pos$group
percentage1 <- round(pca@modelDF$R2X,3)*100
percentage2 <- paste0("PC",1:3,"(",percentage1[1:3],"%",")")

colors <- c("BOT"="#FF7F00", "EOC"="#008B8B","HC"= "#104E8B","QC"="#CD0000")
df_pca$color <- colors[df_pca$class0]
scatter3D(x = df_pca[,1], y = df_pca[,2], z = df_pca[,3],
          xlab = paste0("\nPC",1,"(",percentage1[1],"%",")"), 
          ylab = paste0("\nPC",2,"(",percentage1[2],"%",")"), 
          zlab =  paste0("\nPC",3,"(",percentage1[3],"%",")"), 
          colvar = NULL, col =df_pca$color,pch = 16, cex = 1,
          ticktype = "simple", bty = "b2", theta =350,phi =20,d = 90,cex.lab = 0.8)

#


#neg
rm(df_pca)
df_train_neg <- read.csv("./df_trainset_neg.csv") 

pca <- opls(df_train_neg[,-c(1:2)],crossvalI = 88,predI = 3)
df_pca<- as.data.frame(pca@scoreMN)
df_pca$class0 <-df_train_neg$group
percentage1 <- round(pca@modelDF$R2X,3)*100
percentage2 <- paste0("PC",1:3,"(",percentage1[1:3],"%",")")

colors <- c("BOT"="#FF7F00", "EOC"="#008B8B","HC"= "#104E8B","QC"="#CD0000")
df_pca$color <- colors[df_pca$class0]
scatter3D(x = df_pca[,1], y = df_pca[,2], z = df_pca[,3],
          xlab = paste0("\nPC",1,"(",percentage1[1],"%",")"), 
          ylab = paste0("\nPC",2,"(",percentage1[2],"%",")"), 
          zlab =  paste0("\nPC",3,"(",percentage1[3],"%",")"), 
          colvar = NULL, col =df_pca$color,pch = 16, cex = 1,
          ticktype = "simple", bty = "b2", theta =280,phi =10,d = 90,cex.lab = 0.8)
leg<-    legend("right",
                title = "", 
                legend = unique(df_pca$class0),  
                pch = 21,     
                pt.cex = 2,  
                cex = 1,
                pt.bg = c( "EOC"="#008B8B","BOT"="#FF7F00","HC"= "#104E8B","QC"="#CD0000"), 
                bg = "white",  
                bty = "n",   
                horiz = F   
)
#PLSDA

PLSDA <- function(data0,case,control){
  class <- c(case,control)
  class <- as.factor(class)
  data1<- data0[data0$group %in% c(class), ]
  data1$group <- factor(data1$group,levels = class)
  n <- length(data1$group )
  plsda  <- opls(data1[,-c(1:2)], data1$group, predI = 2, permI = 20,crossvalI = n-1)
  df_plsda <- as.data.frame(plsda@scoreMN)
  df_plsda$class0 <- data1$group
  df_perm <- as.data.frame(plsda@suppLs$permMN)
  #
  
  vip <- plsda@vipVn  
  Pvalue <- vector()
  log2FC <- vector()
  
  for(i in 3:ncol(data1)){
    t <- t.test(data1[,i] ~  data1$group, data=data1)
    
    Pvalue[i-2] <- t$p.value
    log2FC[i-2] <-t$estimate[[paste0("mean in group ", class[1])]]-t$estimate[[paste0("mean in group ", class[2])]]
  }
  
  
  dif_ion<- data.frame(vip,Pvalue,log2FC)
  dif_ion$change <- "NS"
  dif_ion$change[which(dif_ion$vip > 1 & dif_ion$Pvalue <0.05 & dif_ion$log2FC  > log2(1.5))] <- "Up"
  dif_ion$change[which(dif_ion$vip > 1 & dif_ion$Pvalue <0.05 & dif_ion$log2FC  < -log2(1.5))] <- "Down"
  
  #plsda 
  percentage3 <- round(plsda@modelDF$R2X,2)*100
  percentage4 <- paste0("PC",1:2,"(",percentage3[1:2],"%",")")
  plot_plsda0 <- ggplot(df_plsda,aes(x = p1,y = p2, color = class0)) + 
    geom_point(size=2,alpha = 0.8) +
    labs(title = paste(class[1]," vs.",class[2],seq=" ")) +
    labs(color = "Class")+
    theme_bw()+theme(panel.grid=element_blank())+
    xlab(percentage4[1]) +
    ylab(percentage4[2]) +
    ylim(-130,130)+
    scale_color_manual(values = c("BOT"="#FF7F00", "EOC"="#008B8B","HC"= "#104E8B"))+
    stat_ellipse(linetype = 2,level=0.95)+
    annotate("text",label=paste("R2Y=",plsda@summaryDF[["R2Y(cum)"]],seq="")
             ,x=30,y=-100,hjust=0,size=2)+
    annotate("text",label=paste("Q2=",plsda@summaryDF[["Q2(cum)"]],seq="")
             ,x=30,y=-100-20,hjust=0,size=2)
  #permutation test
  colnames(df_perm)[1:3] <- c("R2X","R2Y","Q2")
  df_perm2 <- df_perm[,c(2,3,7)]
  df_perm3 <- reshape2::melt(df_perm2,measure.vars = c("R2Y","Q2"))
  R2_y <- df_perm$R2Y[1]
  Q2_y <- df_perm$Q2[1]
  
  y1 <- df_perm$R2Y - R2_y
  x1 <- df_perm$sim -1
  y1 <- y1[2:length(y1)]
  x1 <- x1[2:length(x1)]
  coefficients(lm(y1 ~ x1+0))
  
  y2 <- df_perm$Q2 - Q2_y
  x2 <- df_perm$sim -1
  y2 <- y2[2:length(y2)]
  x2 <- x2[2:length(x2)]
  coefficients(lm(y2 ~ x2+0))
  plot_perm0 <- ggplot(data = df_perm3,aes(x = sim, y = value, colour = variable))+
    geom_point(size = 1,alpha = 0.6)+
    scale_color_manual(values= c("#00CD66","#1E90FF"))+
    theme_classic()+
    xlab(expression("Similarity(y,"~"y"[perm]*")"))+
    ylab(expression("R2Y/Q2"))+ 
    geom_abline(slope = coefficients(lm(y1 ~ x1+0)), intercept = R2_y -coefficients(lm(y1 ~ x1+0)), linewidth = 0.5, color = "#00CD66",linetype = "dashed")+
    geom_abline(slope = coefficients(lm(y2 ~ x2+0)), intercept = Q2_y -coefficients(lm(y2 ~ x2+0)),linewidth = 0.5, color = "#1E90FF",linetype = "dashed")+
    scale_x_continuous(limits = c(-0,1))+
    geom_hline(yintercept = 0,linewidth = 0.5  )+
    guides(colour = guide_legend(reverse = TRUE))+
    
    theme(legend.title = element_blank())
  
  
  
  
  
  result <- list(plsda_plot=plot_plsda0,permutation_plot= plot_perm0 ,dif_ion=dif_ion)
  return(result)
  
}
#pos
pos_EOC_HC <- PLSDA( df_train_pos,"EOC","HC")
pos_BOT_HC <- PLSDA( df_train_pos,"BOT","HC")
pos_EOC_BOT <- PLSDA( df_train_pos,"EOC","BOT")
#neg
neg_EOC_HC <- PLSDA( df_train_neg,"EOC","HC")
neg_BOT_HC <- PLSDA( df_train_neg,"BOT","HC")
neg_EOC_BOT <- PLSDA( df_train_neg,"EOC","BOT")
#
dir.create("./multivariate_statistics_result")
dir.create("./multivariate_statistics_result/plsda_plot/")
dir.create("./multivariate_statistics_result/plsda_plot/permutation_plot/")
dir.create("./multivariate_statistics_result/differential_table/")
save_plasd_result <- function(x){
  ggsave(paste0("./multivariate_statistics_result/plsda_plot/",x,".tiff"),
         get(x)$plsda_plot,width = 6,height=5)
  ggsave(paste0("./multivariate_statistics_result/plsda_plot/permutation_plot/",x,".tiff"),
         get(x)$permutation_plot,width = 6,height=5)
  write.csv(get(x)$dif_ion,
            paste0("./multivariate_statistics_result/differential_table/",x,"_table.csv"))
}
save_plasd_result ("neg_EOC_BOT")
#########volcano_plot#######
Volcano <- function(pospath,negpath,class){
  
  posvol <- read.csv(pospath,header = T,row.names = 1,stringsAsFactors = F)
  
  negvol <- read.csv(negpath,header = T,row.names = 1,stringsAsFactors = F)
  
  mycol <-c("#1874CD", "#D1D1D1", "#EE6363")
  dt_vol <- rbind(posvol,negvol)
  
  dt_vol$"-log10Pvalue" <- -1*log10(dt_vol$Pvalue)
  
  plot_vol0  <- ggplot(dt_vol, aes(x = log2FC, y = `-log10Pvalue`))+ 
    theme_bw()+theme(panel.grid=element_blank())+
    geom_point(data = subset(dt_vol, change =="NS"), aes(color = change), size=1,alpha = 0.8) +
    geom_point(data = subset(dt_vol, change != "NS"), aes(color = change), size=1,alpha = 0.8) +
    labs(title = paste(class[1]," vs.",class[2],seq=" "))+
    scale_color_manual(values = mycol ) +
    theme(legend.title = element_blank())+
    xlab(expression("Log"[2]~"(Fold change)"))+
    ylab(expression("-Log"[10]~"P-value"))+
    geom_vline(xintercept = log2(1.5),linetype = 6,linewidth = 0.5) + 
    geom_hline(yintercept = -log10(0.05),linetype = 6,linewidth = 0.5) + 
    geom_vline(xintercept = -log2(1.5),linetype = 6,linewidth = 0.5)+
    xlim(-7,7)+theme(legend.title = element_blank())
  return( plot_vol0)
}
vol_EOC_HC <- Volcano("./multivariate_statistics_result/differential_table/pos_EOC_HC_table.csv","./multivariate_statistics_result/differential_table/neg_EOC_HC_table.csv",c("EOC","HC"))
vol_BOT_HC <- Volcano("./multivariate_statistics_result/differential_table/pos_BOT_HC_table.csv","./multivariate_statistics_result/differential_table/neg_BOT_HC_table.csv",c("BOT","HC"))
vol_EOC_BOT <- Volcano("./multivariate_statistics_result/differential_table/pos_EOC_BOT_table.csv","./multivariate_statistics_result/differential_table/neg_EOC_BOT_table.csv",c("EOC","BOT"))
