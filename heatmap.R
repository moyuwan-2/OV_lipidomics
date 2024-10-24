########heatmap#####
library(tidyverse )
library(pheatmap)
lipid_data <- readxl::read_xlsx("./heatmap/lipid_data.xlsx") %>% subset(select!=0) %>% subset(dif_change!="D")  
sample_pheno <- read.csv("./heatmap/sample_pheno.csv",row.names = 1) %>% subset(set%in%"train")
col<- c("#EEA9B8", "#DDA0DD", "#7AC5CD", "#7CCD7C", "#CDC673", "#EE0000", "#9B30FF", "#4876FF", "#00CED1", "#00FF7F", "#00CD00", "#FFD700", "#FF8C00", "#A52A2A", "#483D8B", "#00868B", "#698B22", "#53868B","#CDB38B")
Heatmap <- function(colname ,groups){
  
  f1_1 <- f1 %>% subset(get(colname) %in%c("Up","Down")) %>% column_to_rownames(var ="Metabolite name")
  lev <-unique(f1_1$class) 
  f2 <-  f1_1 %>% .[,-c(1:22)]
  
  
  Md <- sample_pheno %>% select(all_of(c("group"))) %>% subset(group%in%groups)
  Md$group <- factor(Md$group,levels =groups)
  
  
  colnames(Md) <- "Group"
  
  cdMB <- Md %>% mutate(na=1)
  l <- c("HC","BOT","EOC")
  cdMB$Group <- factor(cdMB$Group,levels = l)
  cdMB <-cdMB[order(cdMB$Group),] %>%select(all_of("Group"))
  
  fMB<-subset(f2,select = c(rownames(cdMB)))
  
  cdlipid <- data.frame(f1$class,f1$`Metabolite name`) 
  row.names(cdlipid) <- f1$`Metabolite name`
  colnames(cdlipid) <- c("Category","name")
  
  
  
  
  cdlipid$Category <- factor(cdlipid$Category,levels=lev)
  cdlipid1<- cdlipid[order(cdlipid$Category),] %>% subset(name%in%rownames(fMB) )
  cdlipid2 <- data.frame(cdlipid1$Category)
  
  
  row.names(cdlipid2) <- cdlipid1$name
  colnames(cdlipid2)[1] <- "Category"
  
  
  
  cd<- unique(cdMB$Group) %>% as.character()
  
  cd_col <- c(BOT="#FFD700", EOC="#00EE76",HC= "#00BFFF" )
  lip_col <- c(FA=col[1],ST=col[10],SP=col[2],HexCer=col[17],
               CAR=col[3],Cer=col[16],DG=col[4],PI=col[15],
               TG=col[5],LPC=col[14],PE=col[6],LPE=col[13],
               "Ether PC"=col[7],"Ether LPC"=col[12],"Ether PE"=col[8],SM=col[11], CE=col[9],OxFA=col[18])
  
  
  coul<- colorRampPalette(c("#104E8B","white","#FF4040"))(100)
  
  anno_color <- list(Group=c(cd_col[cd[1]],cd_col[cd[2]]),
                     Category=c(lip_col[lev[1:length(lev)]]))
  fMB1 <- lapply(fMB, function(x) as.numeric(as.character(x)))
  
  fMB1 <- as.data.frame(fMB1) 
  rownames(fMB1 ) <-row.names(fMB) 
  pMB<- pheatmap(fMB1,
                 scale = "row",
                 color = coul,
                 clustering_method = "complete",
                 cluster_cols = F,treeheight_col = 0,cutree_cols =2,
                 cluster_rows = F,treeheight_row = 0,cutree_rows =2, 
                 annotation_col = cdMB,
                 annotation_row = cdlipid2,
                 annotation_legend  =T,
                 show_rownames = F,
                 show_colnames = F,
                 annotation_colors = anno_color,
                 fontsize_row = 7,
                 fontsize_col = 7 , 
                 labels_col = " ",
                 border_color = NA,
                 legend = T
  )
  
  return(pMB)
}
Heat_MH <- Heatmap (colname ="change_EOC_HC",groups=c("EOC","HC"))
Heat_BH <- Heatmap (colname ="change_BOT_HC",groups=c("BOT","HC"))
Heat_MB <- Heatmap (colname ="change_EOC_BOT",groups=c("EOC","BOT"))
