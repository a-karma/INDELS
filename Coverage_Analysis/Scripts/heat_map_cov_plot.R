# May need to modify this accordingly to your filesystem
# setwd("INDELS/Coverage_Analysis/Scripts/")

All_region_df<-"./Output/edited_coverage_database.csv"

library(ggplot2)
library(reshape)
library(scales)
library(wesanderson)
library(ggtext)
library(tidyverse)

tab<-read.table(All_region_df, header = T, sep="\t",dec = ".")
dat1<-data.frame(tab$Sample_ID,tab$Region,tab$E.DoC.,tab$Type)
names(dat1)<-c("ID","Region","Mean_Coverage","Type")
pal <- wes_palette("Zissou1", 21, type = "continuous")

for (i in seq(length(dat1$Mean_Coverage))){
if ( dat1$Mean_Coverage[i] > 20 ){
  dat1$Mean_Coverage[i]=20
}}

pdf(paste("./Plots/heat_map_mean_doc.pdf", sep="", collapse=NULL), width=7, height=21, pointsize=10)
print(
  ggplot(dat1, aes(x = Region, y = ID, fill = Mean_Coverage)) +
        geom_tile() +
        scale_fill_gradientn(colours = pal ,values = rescale(c(0,0.01,0.1,0.25,0.5,1))) +
        theme(axis.text = element_text(size=12), axis.title = element_blank())+
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)))
dev.off()



dat2<-data.frame(tab$Sample_ID,tab$Region,tab$BoC)
names(dat2)<-c("ID","Region","Breadth")

pdf(paste("./Plots/heat_map_boc.pdf", sep="", collapse=NULL), width=7, height=21, pointsize=10)
print(
  ggplot(dat2, aes(x = Region, y = ID, fill =Breadth)) +
    geom_tile() +
    scale_fill_gradientn(colours = pal ,values = rescale(c(0,0.01,0.1,0.2,0.4,1))) +
    theme(axis.text = element_text(size=12), axis.title = element_blank())+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)))
dev.off()
