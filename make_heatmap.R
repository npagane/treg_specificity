suppressPackageStartupMessages({
  library(limma)
  library(Glimma)
  library(edgeR)
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(pheatmap)
  library(gplots)
  library(viridis)
  library(cluster)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(dendextend)
  library(vioplot)
  library(MASS)
})

# read in the data 
#setwd("~/Dropbox (MIT)/wong/projects/collabs/dave/")
df <- read.csv("new_tempmat_tconv_noperipheralization.csv")

mat <- as.matrix(df[,1:8])

mat_scale <- scale(mat)

df$cluster <- factor(df$cluster)

ylabs <- c("Ki67 on local Tregs", "PD-1 on local Tregs", "pSTAT5 on local Tregs", "Ki67 on MJ23 Tconv","PD-1 on MJ23 Tconv", "pSTAT5 on MJ23 Tconv", "local Treg density", "local MJ23 Treg number")

bk <- seq(-1.5,1.5,by=0.1)
color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(bk))

#temp_markers = ["C1", "C2", "A1", "C3", "A2"]
#temp_colors =  ["plum", "purple", "goldenrod", "k", "tan"] 
ann_colors = list(
  treatment = c("WT"="royalblue", "dTEC" = 'red'),
  cell = c("MJ23 Treg" = "black", "polyclonal Treg" = "white"),
  cluster = c("C1"= "purple", "C2"= "black", "A1"="goldenrod", "BA" = "grey", "A2"="wheat")
)

pheatmap(t(mat_scale),clustering_distance_cols = 'correlation', color=color, breaks=bk, 
         cluster_rows = F, labels_row=ylabs,  cutree_cols=6, clustering_method = "average",
         annotation_col=df[,c('treatment','cluster')], annotation_colors = ann_colors,
         treeheight_col = 100, legend_breaks = seq(-1.5,1.5,by=0.5), row_order = colnames(mat),
         fontsize=15, cellwidth=2.0, cellheight=40, row_split= c(rep("mean local Treg",3), rep("MJ23 Tconv", 3), rep("local env",2)))


#############################

alltreg <- read.csv("temp_tregs.csv")
mjtconv <- read.csv("temp_mjtconvs.csv")
otiitconv <- read.csv("temp_otiitconvs.csv")
mjtregs <- read.csv("temp_mjtregs.csv")


ggplot(data=alltreg)+
  stat_density_2d(aes(x=x, y=y,fill=stat(level)),h=c(75,75),geom="polygon")+
  geom_point(data=mjtconv,aes(x=x,y=y), colour="white",size=2)+
  geom_point(data=otiitconv,aes(x=x,y=y), colour="green",size=2)+
  geom_point(data=mjtregs,aes(x=x,y=y),colour="darkorange",size=1.75)+
  scale_fill_viridis(trans="log10",option="mako",alpha=0.9)+
  coord_fixed()+
  theme_bw()+theme(axis.ticks.length=unit(0.3,"cm"))


