

# Set Library
library(ggplot2)
library(Seurat) # 4.3.0

# Calculation after normalization, dimentional reduction, clustering and UMAP


# Set the signature genes as list
exhaustion_markers<-
  list(c("PDCD1","TOX","CXCL13","TIGIT","CTLA4","TNFRSF9","HAVCR2","LAG3"))

cytokine_effector_markers<-
  list(c("CST7","GZMK","GZMA","NKG7","IFNG","PRF1","GZMB","GNLY"))
naiive_markers<-
  list(c("TCF7","LEF1","CCR7","SELL","MALL"))


# Calculate the Module score 
lung<-AddModuleScore(lung,features=naiive_markers, name = "naiive_markers")
lung<-AddModuleScore(lung,features=exhaustion_markers, name = "exhaustion_markers")
lung<-AddModuleScore(lung,features=cytokine_effector_markers, name = "cytokine_effector_markers")
lung<-AddModuleScore(lung,features=terminal_exhausted_markers,name="terminal_exhausted_markers")

# Compare the T cell status between Cluster5 and cluster11
df_immunescore<-lung@meta.data
df_immunescore$Cluster<-paste("Cluster",df_immunescore$Cluster,sep="")
df_immunescore_cls5_11<-subset(df_immunescore,Cluster=="Cluster5"|Cluster=="Cluster11")
df_immunescore_cls5_11$Cluster<-factor(df_immunescore_cls5_11$Cluster,levels=c("Cluster5","Cluster11"))

# Perform the statistical test
cls5<-subset(df_immunescore_cls5_11,Cluster=="Cluster5")
cls11<-subset(df_immunescore_cls5_11,Cluster=="Cluster11")
naiive_wilcoxon<-wilcox.test(x=cls5$Naiive,y=cls11$Naiive,paired=F)

# Visualization for Naiive Tcell Scores as an example

Naiive_plot<-ggplot(data=df_immunescore_cls5_11,aes(x=Cluster,y=Naiive,fill=Cluster))+
  geom_violin()+
  labs(title="Naiive")+
  xlab(NULL)+
  ylab(label="Module Score")+
  scale_y_continuous(limits = c(-0.6,1.3),breaks = seq(-0.5,1.0,0.5)) +
  scale_fill_manual(values = col_list[c(6,12)])+
  geom_text(x = 1.5, y = 1.25, label = "***",size=12)+
  geom_segment(x = 1, xend = 1, y = 1.2, yend = 1.1) +
  geom_segment(x = 1, xend = 2, y = 1.2, yend = 1.2) +
  geom_segment(x = 2, xend = 2, y = 1.2, yend = 1.1)+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,face="bold",size = 28),
        axis.title=element_text(size=18,face="bold"),
        axis.text.x=element_text(size=16,color = "black"),
        axis.text.y=element_text(size=16,color = "black")
  )
plot(Naiive_plot)
