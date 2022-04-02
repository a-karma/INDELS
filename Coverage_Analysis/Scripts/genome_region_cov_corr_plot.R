library(ggplot2)
library(ggrepel)

tab<-read.table(region_df, header = T, sep="\t",dec = ".")
colnames(tab)<-c("ID","Region","Length","Mean_Cov_Depth","Breadth","Genome_Depth","Type")

dat<-data.frame(tab$ID,tab$Mean_Cov_Depth,tab$Genome_Depth,tab$Type)
names(dat)<-c("ID","R_doc","G_doc","Type")
x_label="Mean genome Depth of Coverage"
y_label="Mean region Depth of Coverage"

col<-c("#1A80C4","#CC3D3D")

pdf(paste("Plots/",levels(tab$Region),"_",levels(as.factor(tab$Length)),"_bp_corr_plot.pdf", sep="", collapse=NULL))
print(
ggplot(dat,aes(G_doc,R_doc))+
  geom_point(aes(colour = Type), size=3, alpha=0.7) +
  labs(x=x_label , y=y_label)+
  scale_colour_manual(values=col) +
  geom_abline(intercept=0, slope=1)+
  geom_abline(intercept=0, slope=2, linetype="dashed")+
  geom_abline(intercept=0, slope=4, linetype="dashed")+
  geom_abline(intercept=0, slope=0.5, linetype="dashed")+
  geom_abline(intercept=0, slope=0.25, linetype="dashed")+
  geom_label_repel(aes(label = ID), box.padding = 0.9, point.padding=0.2, size=2.5)+
  ggtitle(paste("Coverage Depth correlation plot: ",region)))
dev.off()

