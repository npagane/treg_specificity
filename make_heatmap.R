suppressPackageStartupMessages({
  library(limma)
  library(Glimma)
  library(edgeR)
  library(dplyr)
  library(Mus.musculus)
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
df <- read.csv("")

mat <- as.matrix(df[,1:(ncol(df)-2)])

mat <- mat[,c(4:6, 1:3, 7:8)]


mat_scale = scale(mat)

df$cluster <- factor(df$cluster)

ylabs <- c("Ki67 on MJ Tconv","PD-1 on MJ Tconv", "pSTAT5 on MJ Tconv", "Ki67 sum on local Tregs", "PD-1 sum on local Tregs", "pSTAT5 sum on local Tregs", "local Treg density", "local MJ Treg number")


library(circlize)

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("green", "white", "pink"), hcl_palette = 'PiYG')
col_fun(seq(-3, 3))


ha = HeatmapAnnotation(
  treatment = c("WT"="royalblue", "dTEC" = 'red'),
  cell = c("MJ Treg" = "black", "polyclonal Treg" = "white"),
  cluster = c("A1"= "chocolate", "A2"= "goldenrod", "C1"="mediumpurple", "C2"="purple",  "C3"="cornflowerblue", "BA"= "grey"),
  annotation_name_side = "left")

Heatmap(t(mat_scale),cluster_rows=FALSE,clustering_distance_columns = "pearson",
        clustering_method_columns = "average", col=col_fun,column_dend_height = 100,
        column_dend_reorder = TRUE, row_labels = ylabs, top_annotation= ha)

Heatmap(t(mat_scale),cluster_rows=FALSE,clustering_distance_columns = "pearson",
        clustering_method_columns = "average",row_labels = ylabs, cluster_column_slices = FALSE,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                            labels = c("group0", "group1", "group2", "group3", "group4", "group5", "group6"), 
                                                            labels_gp = gpar(col = "white", fontsize = 10))),
        column_km = 7,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                         ))
        row_split = c(rep("MJ Tconv",3), rep("local enviornment",2), rep("MJ Treg",3)))

bk <- seq(-1.5,1.5,by=0.1)
color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(bk))

ann_colors = list(
  treatment = c("WT"="royalblue", "dTEC" = 'red'),
  cell = c("MJ Treg" = "black", "polyclonal Treg" = "white"),
  cluster = c("A1"= "chocolate", "A2"= "goldenrod", "C1"="mediumpurple", "C2" = "purple", "C3"="cornflowerblue", "BA"= "grey")
)

pheatmap(t(mat_scale),clustering_distance_cols = 'correlation', color=color, breaks=bk, 
         
         cluster_rows = F, labels_row=ylabs,  cutree_cols=10, clustering_method = "average",
         annotation_col=df[,c('treatment','cluster')], annotation_colors = ann_colors,
         treeheight_col = 100, legend_breaks = seq(-1.5,1.5,by=0.5), row_order = colnames(mat),
         fontsize=15, cellwidth=2, cellheight=40, row_split= c(rep("MJ Tconv",3), rep("Treg", 3), rep("local env",2)))


#############################


alltreg <- read.csv("old_alltregs.csv")
tconvs <- read.csv("old_tconvs.csv")
mjtregs <- read.csv("old_mjtregs.csv")

ggplot(data=alltreg[alltreg$n=="n=2",])+
  stat_density_2d(aes(x=x, y=y,fill=stat(level)),h=c(75,75),geom="polygon")+
  geom_point(data=tconvs[tconvs$n=='n=2',],aes(x=x,y=y),width=50,height=10,colour="red",size=1.5)+
  geom_point(data=mjtregs[mjtregs$n=='n=2',],aes(x=x,y=y),width=50,height=10,colour="blue",size=1.5)+
  scale_fill_viridis(trans="log10",option="mako",alpha=0.9)+
  coord_fixed()+
  theme_bw()+theme(axis.ticks.length=unit(0.3,"cm"))


alltreg <- read.csv("transfers_alltregs.csv")
mjtconv <- read.csv("transfers_mjtconv.csv")
wttconv <- read.csv("transfers_wttconv.csv")

ggplot(data=alltreg[alltreg$n=="2-2",])+
  stat_density_2d(aes(x=x, y=y,fill=stat(level)),h=c(75,75),geom="polygon")+
  geom_point(data=mjtconv[mjtconv$n=='2-2',],aes(x=x,y=y),width=50,height=10,colour="red",size=1.5)+
  geom_point(data=wttconv[wttconv$n=='2-2',],aes(x=x,y=y),width=50,height=10,colour="darkorange",size=1.5)+
  scale_fill_viridis(trans="log10",option="mako",alpha=0.9)+
  coord_fixed()+
  theme_bw()+theme(axis.ticks.length=unit(0.3,"cm"))


alltreg <- read.csv("temp_tregs.csv")
mjtconv <- read.csv("temp_tconvs.csv")
mjtregs <- read.csv("temp_mjtregs.csv")


ggplot(data=alltreg)+
  stat_density_2d(aes(x=x, y=y,fill=stat(level)),h=c(75,75),geom="polygon")+
  geom_point(data=mjtconv[mjtconv$Ki67.SFI.normalized>=4,],aes(x=x,y=y),width=50,height=10,colour="yellow1",size=2)+
  geom_point(data=mjtconv[mjtconv$Ki67.SFI.normalized<4,],aes(x=x,y=y),colour="white",size=2)+
  geom_point(data=mjtregs,aes(x=x,y=y),colour="blue",size=1.75)+
  scale_fill_viridis(trans="log10",option="mako",alpha=0.9)+
  coord_fixed()+
  theme_bw()+theme(axis.ticks.length=unit(0.3,"cm"))

alltreg <- read.csv("temp_tregs.csv")
mjtconv <- read.csv("temp_tconvs.csv")
#mjtregs <- read.csv("temp_mjtregs.csv")

ggplot(data=alltreg)+
  stat_density_2d(aes(x=x, y=y,fill=stat(level)),h=c(75,75),geom="polygon")+
  geom_point(data=mjtconv[mjtconv$Ki67.SFI.normalized>=4,],aes(x=x,y=y),width=50,height=10,colour="gold",size=2)+
  geom_point(data=mjtconv[mjtconv$Ki67.SFI.normalized<4,],aes(x=x,y=y),colour="white",size=2)+
  scale_fill_viridis(trans="log10",option="mako",alpha=0.9)+
  coord_fixed()+
  theme_bw()+theme(axis.ticks.length=unit(0.3,"cm"))

