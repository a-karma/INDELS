library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(caret)
library(factoextra)

# Importing dataframe
setwd("~/Desktop/Willis_Synd/Insert_Size")
full_data<-read.table("final_summary.tsv", header=F, col.names=c(
  "ID","Missing", "Ins_events","Del_events","mean_len_Ins","sd_len_Ins","mean_len_Del","sd_len_Del","total_net_len","Status"))


################################################# Net length distro  MaxMiss 50% #####################################################

# Subsetting the dataset
dat<-subset(full_data, Missing<=99)
nrow(dat) # There are 1667 samples 
wild<-subset(dat, Status=='W')
nrow(wild) # of which only 133 are wild canids

# Calculating weighted mean
mu <- dat %>% 
  group_by(Status) %>%
  summarise(grp.mean = mean(total_net_len))

# plotting Densities separating Wild and Domestic
a<-ggplot(dat, aes(x = total_net_len))
#Change line color by sex
a + geom_density(aes(color = Status)) +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))
# Change fill color by sex and add mean line
# Use semi-transparent fill: alpha = 0.4
a + geom_density(aes(fill = Status), alpha = 0.4) +
  geom_vline(aes(xintercept = grp.mean, color = Status),
             data = mu, linetype = "dashed") +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))+
  scale_fill_manual(values = c("#868686FF", "#EFC000FF"))

##################################################### Net length distro  MaxMiss 25% #################################################

# Subsetting the dataset
dat<-subset(full_data, Missing<=50)
nrow(dat) # There are 1112 samples 
wild<-subset(dat, Status=='W')
nrow(wild) # of which only 39 are wolves

# Calculating weighted mean
mu <- dat %>% 
  group_by(Status) %>%
  summarise(grp.mean = mean(total_net_len))

# plotting Densities separating Wild and Domestic
a<-ggplot(dat, aes(x = total_net_len))
#Change line color by sex
a + geom_density(aes(color = Status)) +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))
# Change fill color by sex and add mean line
# Use semi-transparent fill: alpha = 0.4
a + geom_density(aes(fill = Status), alpha = 0.4) +
  geom_vline(aes(xintercept = grp.mean, color = Status),
             data = mu, linetype = "dashed") +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))+
  scale_fill_manual(values = c("#868686FF", "#EFC000FF"))


####################################################### Insertion distro Maxmiss 25% ####################################################

# Density of the average insertion length
a<-ggplot(dat, aes(x = mean_len_Ins))
a + geom_density(fill = "lightgray") +
  geom_vline(aes(xintercept = mean(mean_len_Ins)), 
             linetype = "dashed", size = 0.6, color = "#FC4E07")

# Calculating weighted mean
mu <- dat %>% 
  group_by(Status) %>%
  summarise(grp.mean = mean(mean_len_Ins))
# Checking mean of distros
mu

# plotting Densities separating Wild and Domestic
a<-ggplot(dat, aes(x = mean_len_Ins))
#Change line color by sex
a + geom_density(aes(color = Status)) +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))
# Change fill color by sex and add mean line
# Use semi-transparent fill: alpha = 0.4
a + geom_density(aes(fill = Status), alpha = 0.4) +
  geom_vline(aes(xintercept = grp.mean, color = Status),
             data = mu, linetype = "dashed") +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))+
  scale_fill_manual(values = c("#868686FF", "#EFC000FF"))


###################################################### Deletion distro MaxMiss 25% ###################################################

a<-ggplot(dat, aes(x = mean_len_Del))
a + geom_density(fill = "lightgray") +
  geom_vline(aes(xintercept = mean(mean_len_Del)), 
             linetype = "dashed", size = 0.6, color = "#FC4E07")

# Calculating weighted mean
mu <- dat %>% 
  group_by(Status) %>%
  summarise(grp.mean = mean(mean_len_Del))
# Checking mean of distros
mu

# plotting Densities separating Wild and Domestic
a<-ggplot(dat, aes(x = mean_len_Del))
#Change line color by sex
a + geom_density(aes(color = Status)) +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))
# Change fill color by sex and add mean line
# Use semi-transparent fill: alpha = 0.4
a + geom_density(aes(fill = Status), alpha = 0.4) +
  geom_vline(aes(xintercept = grp.mean, color = Status),
             data = mu, linetype = "dashed") +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))+
  scale_fill_manual(values = c("#868686FF", "#EFC000FF"))


##################################################### Clustering MaxMiss 50% #########################################################

#==================================== clustering on the whole spectrum, no transformation =====================================
dat<-subset(full_data, Missing<=99 )
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
dist_matrix <- dist(numeric_dat, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")

species_color=c("#800080", "#EFC000FF")
p0a <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: MaxMiss 50%, No transformation, var=(net, mean, sd)")


dat<-subset(full_data, Missing<=99 )
numeric_dat<-select(dat,total_net_len, Ins_events,Del_events)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
dist_matrix <- dist(numeric_dat, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")

species_color=c("#800080", "#EFC000FF")
p0b <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: MaxMiss 50%, No transformation, var=(net, events)")

#==================================== clustering of the extremes, no transformations ===========================================

dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)

numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
dist_matrix <- dist(numeric_dat, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")

species_color=c("#800080", "#EFC000FF")
p1a <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") +
  ggtitle("Indels Dendrogram: No Transformations, var=(net, mean, sd), tails only")

dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len, Ins_events,Del_events)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
dist_matrix <- dist(numeric_dat, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")

species_color=c("#800080", "#EFC000FF")
p1b <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: MaxMiss 50%, No Transformations, var=(net, events), tails only")

#==================================== with Normalized values ==========================================
dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)

pre_data_norm <- preProcess(numeric_dat,method=c("range"))
data_norm<-predict(pre_data_norm,numeric_dat)

dist_matrix <- dist(data_norm, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")

species_color=c("#800080", "#EFC000FF")
p2a <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: MaxMiss 50%, Normalised var=(net, mean, sd), tails only")

dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len, Ins_events,Del_events)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)

pre_data_norm <- preProcess(numeric_dat,method=c("range"))
data_norm<-predict(pre_data_norm,numeric_dat)

dist_matrix <- dist(data_norm, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")

species_color=c("#800080", "#EFC000FF")
p2b <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: MaxMiss 50%, Normalised var=(net, events), tails only")


#==================================== with Standirdised values ====================================================

dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
pre_data_std <- preProcess(numeric_dat,method=c("center", "scale"))
data_std<-predict(pre_data_std,numeric_dat)

dist_matrix <- dist(data_std, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")

species_color=c("#800080", "#EFC000FF")
p3a <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: MaxMiss 50%, Standardised var=(net_len, mean, sd), tails only")

dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len, Ins_events,Del_events)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
pre_data_std <- preProcess(numeric_dat,method=c("center", "scale"))
data_std<-predict(pre_data_std,numeric_dat)

dist_matrix <- dist(data_std, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")

species_color=c("#800080", "#EFC000FF")
p3b <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: MaxMiss 50%, Standardised var=(net_len, events), tails only")

#==================================== Plots ========================================================

pdf("Hierarchical_Clustering_plots_MaxMiss50.pdf", width=15, height=17)
print(p0a)
print(p0b)
print(p1a)
print(p1b)
print(p2a)
print(p2b)
print(p3a)
print(p3b)
dev.off()

##################################################### Clustering MaxMiss 25% #########################################################

# Min Indel size 5 bp =============================================================================================================
full_data<-read.table("final_summary_min5.tsv", header=F, col.names=c(
  "ID","Missing","Ref_Alleles", "Discarded","Ins_events","Del_events","mean_len_Ins",
  "sd_len_Ins","mean_len_Del","sd_len_Del","total_net_len","Status"))

species_color<-c(rainbow(16))

dat<-subset(full_data, Missing<=50)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)

pre_data_norm <- preProcess(numeric_dat,method=c("range"))
data_norm<-predict(pre_data_norm,numeric_dat)

dist_matrix <- dist(data_norm, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")


P_25a <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: Min Indel size 5 bp, MaxMiss 25%, Normalised var=(net, mean, sd), All Ind")

dat<-subset(full_data, Missing<=50 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)

pre_data_norm <- preProcess(numeric_dat,method=c("range"))
data_norm<-predict(pre_data_norm,numeric_dat)

dist_matrix <- dist(data_norm, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")


P_25b <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: Min Indel size 5 bp, MaxMiss 25%, Normalised var=(net, mean, sd), tails only")

# Min Indel Size 10 bp ==========================================================================================================
full_data<-read.table("final_summary_min10.tsv", header=F, col.names=c(
  "ID","Missing","Ref_Alleles", "Discarded","Ins_events","Del_events","mean_len_Ins",
  "sd_len_Ins","mean_len_Del","sd_len_Del","total_net_len","Status"))

dat<-subset(full_data, Missing<=50)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)

pre_data_norm <- preProcess(numeric_dat,method=c("range"))
data_norm<-predict(pre_data_norm,numeric_dat)

dist_matrix <- dist(data_norm, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")


P_25c <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: Min Indel size 10 bp, MaxMiss 25%, Normalised var=(net, mean, sd), All Ind")

dat<-subset(full_data, Missing<=50 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)

pre_data_norm <- preProcess(numeric_dat,method=c("range"))
data_norm<-predict(pre_data_norm,numeric_dat)

dist_matrix <- dist(data_norm, method = "euclidean")
dendrogram <- as.dendrogram(hclust(dist_matrix, method = "complete"))
dend_data <-dendro_data(dendrogram)
dend_segments<-dend_data$segments
dend_ends<-dend_segments %>%
  filter(yend==0) %>% left_join(dend_data$labels, by="x") %>%
  rename(ID=label) %>%
  left_join(metatdata, by= "ID")


P_25d <- ggplot() +
  geom_segment(data = dend_segments, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = dend_ends, aes(x=x, y=y.x, xend=xend, yend=yend, color = Status))+
  scale_color_manual(values = species_color) +
  scale_y_reverse() +
  coord_flip() + theme_bw() + theme(legend.position = "none") + ylab("Distance") + 
  ggtitle("Indels Dendrogram: Min Indel size 10 bp, MaxMiss 25%, Normalised var=(net, mean, sd), tails only")


# Plots =====================================================================================================================
pdf("Hierarchical_Clustering_plots_MaxMiss25.pdf", width=15, height=17)
print(P_25a)
print(P_25b)
print(P_25c)
print(P_25d)
dev.off()

################################################################ PCA plots #############################################################

setwd("~/Desktop/Willis_Synd/Insert_Size")

# Min Indel size ====================================================================================================================
full_data<-read.table("final_summary_min1.tsv", header=F, col.names=c(
  "ID","Missing","Ref_Alleles", "Discarded","Ins_events","Del_events",
  "mean_len_Ins","sd_len_Ins","mean_len_Del","sd_len_Del","total_net_len","Status"))

dat<-subset(full_data, Missing<=99)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml1_mm50_all<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none", mean.point = FALSE,geom.ind = "point", 
                               pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >1bp, missing_data<50%, All ind")

dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml1_mm50_tails<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point", 
                                 pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >1bp, missing_data<50%, tails only") 
             
dat<-subset(full_data, Missing<=50)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml1_mm25_all<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                               pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >1bp, missing_data<25%, All ind")

dat<-subset(full_data, Missing<=50 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml1_mm25_tails<-fviz_pca_ind(pca_dat, col.ind = dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                                 pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >1bp, missing_data<25%, tails only") 

# Min Indel size 2bp ==================================================================================================================

full_data<-read.table("final_summary_min2.tsv", header=F, col.names=c(
  "ID","Missing","Ref_Alleles", "Discarded","Ins_events","Del_events","mean_len_Ins",
  "sd_len_Ins","mean_len_Del","sd_len_Del","total_net_len","Status"))

dat<-subset(full_data, Missing<=99)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml2_mm50_all<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                               pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >2bp, missing_data<50%, All ind")

dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml2_mm50_tails<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                                 pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >2bp, missing_data<50%, tails only") 


dat<-subset(full_data, Missing<=50)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml2_mm25_all<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                               pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >2bp, missing_data<25%, All ind")

dat<-subset(full_data, Missing<=50 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml2_mm25_tails<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                                 pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >2bp, missing_data<25%, tails only") 
# Min Indel size 5bp ==================================================================================================================

full_data<-read.table("final_summary_min5.tsv", header=F, col.names=c(
  "ID","Missing","Ref_Alleles", "Discarded","Ins_events","Del_events","mean_len_Ins",
  "sd_len_Ins","mean_len_Del","sd_len_Del","total_net_len","Status"))

dat<-subset(full_data, Missing<=99)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml5_mm50_all<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                               pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >5bp, missing_data<50%, All ind")

dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml5_mm50_tails<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                                 pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >5bp, missing_data<50%, tails only") 


dat<-subset(full_data, Missing<=50)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml5_mm25_all<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                               pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >5bp, missing_data<25%, All ind")

dat<-subset(full_data, Missing<=50 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml5_mm25_tails<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                                 pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >5bp, missing_data<25%, tails only")

# Min Indel size 10bp ================================================================================================================

full_data<-read.table("final_summary_min10.tsv", header=F, col.names=c(
  "ID","Missing","Ref_Alleles", "Discarded","Ins_events","Del_events","mean_len_Ins",
  "sd_len_Ins","mean_len_Del","sd_len_Del","total_net_len","Status"))

dat<-subset(full_data, Missing<=99)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml10_mm50_all<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                                pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >10bp, missing_data<50%, All ind")

dat<-subset(full_data, Missing<=99 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml10_mm50_tails<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                                  pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >10bp, missing_data<50%, tails only") 


dat<-subset(full_data, Missing<=50)
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metadata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml10_mm25_all<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE,geom.ind = "point",
                                pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >10bp, missing_data<25%, All ind")

dat<-subset(full_data, Missing<=50 & (total_net_len>=50 | total_net_len<=-50))
numeric_dat<-select(dat,total_net_len,mean_len_Ins,mean_len_Del,sd_len_Ins,sd_len_Del)
rownames(numeric_dat)<-dat$ID
metatdata<-select(dat,ID,Status)
pca_dat<-prcomp(numeric_dat,center = TRUE,scale. = TRUE)

PCA_ml10_mm25_tails<-fviz_pca_ind(pca_dat, col.ind=dat$Status, label = "none",  mean.point = FALSE, geom.ind = "point",
                                  pointshape = 21,fill.ind = dat$Status, alpha.ind = 0.7, legend.title="Categories") + 
  ggtitle("PCA: min_Indel_size >10bp, missing_data<25%, tails only")
fviz_pca_var(pca_dat)
pdf("PCA_plots.pdf")
print(PCA_ml1_mm50_all)
print(PCA_ml1_mm50_tails)
print(PCA_ml1_mm25_all)
print(PCA_ml1_mm25_tails)
print(PCA_ml2_mm50_all)
print(PCA_ml2_mm50_tails)
print(PCA_ml2_mm25_all)
print(PCA_ml2_mm25_tails)
print(PCA_ml5_mm50_all)
print(PCA_ml5_mm50_tails)
print(PCA_ml5_mm25_all)
print(PCA_ml5_mm25_tails)
print(PCA_ml10_mm50_all)
print(PCA_ml10_mm50_tails)
print(PCA_ml10_mm25_all)
print(PCA_ml10_mm25_tails)
dev.off()

########################################################################################################################################

