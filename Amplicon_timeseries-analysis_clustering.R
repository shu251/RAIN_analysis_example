# Import diel data
setwd("~/repos/times_series_methods_scripts/")
load("diel_df_Robj.RData",verbose=T) # load diel dataset
library(reshape2); library(vegan); library(dplyr); library(ggplot2)
#
## Pre-process dataset:
# (1) Remove unwanted taxa
unique(diel_data$Taxa)
rm<-c("Opisthokont-Fungi","Opisthokont-Metazoa","Opisthokonts-Other","Other/unknown","Unassigned")
tmp<-subset(diel_data, !(Taxa %in% rm))
#
# (2) Normalize by relative abundance
relAbun<-tmp %>% 
  group_by(Num) %>% 
  mutate(DNA_relabun=DNA/sum(DNA), RNA_relabun=RNA/sum(RNA)) %>%
  as.data.frame
names(relAbun) # output
#
# (3) Reformat data frame (to wide format) and separate DNA and RNA
names(relAbun)
samplenum<-as.character(c(1:19))
new_samplename<-c("1_6PM","2_10PM","3_2AM","4_6AM","5_10AM","6_2PM","7_6PM","8_10PM","9_2AM","10_6AM","11_10AM","12_2PM","13_6PM","14_10PM","15_2AM","16_6AM","17_10AM","18_2PM","19_6PM")
relAbun$Sample<-factor(relAbun$Num, levels = samplenum, labels = new_samplename)
names(relAbun)
dna_base<-dcast(relAbun[c(1,19,17)], Sample~OTU.ID, fill=0)
rna_base<-dcast(relAbun[c(1,19,18)], Sample~OTU.ID, fill=0)
#
# (4) Format for PCA analysis input
row.names(dna_base)<-dna_base$Sample; dna_base$Sample<-NULL
row.names(rna_base)<-rna_base$Sample; rna_base$Sample<-NULL
# 
# (5) Remove all zero rows
head(dna_base[1:4])
rowsum<-apply(dna_base,1,sum);rowsum
rowsum<-apply(rna_base,1,sum);rowsum
#
#
## NMDS - calculate for RNA and DNA
# (1) Calculate bray-curtis distance, creates a matrix
# First calcualte for RNA
NMDS_rna=metaMDS(rna_base,distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
head(NMDS_rna$points) #given points for MDS plot
NMDS_rna$stress # calculated stress
RNA_pts <- NMDS_rna$points[1:nrow(NMDS_rna$points),]
RNA_pts<-as.data.frame(RNA_pts)
plot(RNA_pts$MDS1, RNA_pts$MDS2)
RNA_pts$Sample<-row.names(RNA_pts)
#
# Repeat for DNA
NMDS_dna=metaMDS(dna_base,distance="bray",k=2,trymax=100,engine=c("monoMDS"),autotransform=FALSE)
head(NMDS_dna$points) #given points for MDS plot
NMDS_dna$stress
DNA_pts <- NMDS_dna$points[1:nrow(NMDS_dna$points),]
DNA_pts<-as.data.frame(DNA_pts)
plot(DNA_pts$MDS1, DNA_pts$MDS2)
DNA_pts$Sample<-row.names(DNA_pts)
#
# Plot:
# Import manual curation of colors for sample names.
## We need this in there
key<-read.delim("key_CCA.txt"); head(key)
key$Sample<-paste(key$Sample.Number, key$Time, sep="_")
label<-as.character(key$Sample)
col<-as.character(key$Color)
shap<-as.numeric(key$Shape)
colScale<-scale_color_manual(values=col)
# Factor:
DNA_pts$order<-factor(DNA_pts$Sample, levels = label)
RNA_pts$order<-factor(RNA_pts$Sample, levels = label)

# DNA
DNA_MDS <- ggplot(DNA_pts, aes(x = MDS1, y = MDS2, label=order, shape=order), color="black") + geom_point(size=6,aes(fill=order)) + theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +scale_shape_manual(values = shap)+scale_fill_manual(values=col)+ geom_hline(yintercept=0, linetype="dashed", color = "#252525")+ geom_vline(xintercept=0, linetype="dashed", color = "#252525")
DNA_MDS+theme(legend.position = "none")
# 
#RNA
RNA_MDS <- ggplot(RNA_pts, aes(x = MDS1, y = MDS2, label=order, shape=order), color="black") + geom_point(size=6,aes(fill=order)) + theme_bw() +theme(panel.grid.minor = element_line(colour="white"), panel.grid.major = element_line(colour = "white")) +scale_shape_manual(values = shap)+scale_fill_manual(values=col)+ geom_hline(yintercept=0, linetype="dashed", color = "#252525")+ geom_vline(xintercept=0, linetype="dashed", color = "#252525")
RNA_MDS+theme(legend.position = "none")
# Make Figure S4
library(cowplot)
figS4<-plot_grid(RNA_MDS+theme(legend.position = "none"), DNA_MDS+theme(legend.position = "none"), align="h", ncol=2, nrow=1,labels=c("A", "B"))
figS4 #svg save W:1000, H: 478
#
#
# CCA analysis:
# First thing -- data exist on a simplex so covariance matrix is singular:
# Evidence:
cov_determinant<-det(as.matrix(rna_base)%*%t(rna_base))
# Next thing, cca w/out constraining variable is just PCA so it is more clear and precise to call it that
# But before we do it we need to get the data out of the simplex. Recommendation is log ratios or
# estimating read counts using a package like dada2.
# Alternatively, you can do a PCoA using your favorite distance metric. For example, Jaccard
rna_dmat<-vegdist(rna_base,method="jaccard")
pcoa_rna<-princomp(rna_dmat)
# Then you'd need to look at the skree plot to determine appropriate # axes.
plot(pcoa_rna)
# This looks like most of the variance is explained by one axis, but the following 2 have similar
# amounts of variance explained, so you'd need to go for 3 axes.
scatterplot3d::scatterplot3d(pcoa_rna$scores[,1:3])
# There is a pretty clear hidden axis here, which seems to suggest some information is excluded from the
# distance matrix provided: recommendation is to use Unifrac distance w/trees of relatedness to capture
# missing phylogenetic information.
# Alternatively, you can use pca/cca if you convert to log ratios.
log_rats<-compositions::ilr(rna_base)
# These are complicated to interpret however because you move from D simplex to D-1 euclidean
# But using these we can see we are now in invertible covariance regime
new_covdet<-det(log_rats%*%t(log_rats))
new_pca<-rda(log_rats)
# Look at this
plot(1:length(new_pca$CA$eig),new_pca$CA$eig,main='Skree Plot',ylab='Eigenvalue',xlab='Axis')
cca_rna<-cca(rna_base, scaling=TRUE)
# In this case eigenvalues are still one big one and then a bunch of small ones, so high dimensional PCA is 
# a preferable visualization to 2D. 
# CCA ploting set up:
# par(xpd=FALSE,oma=c(0,0,0,8)) 
eig<-eigenvals(cca_rna); x<-signif((eig[1]*100),3); y<-signif((eig[2]*100),3) #extract variances for CCA1 and CCA2 axes
plot(cca_rna,type=c("points"),display=c("wa", "sites"),main="",xlab=bquote(.(x)*"%"), ylab=(bquote(.(y)*"%")))
points(cca_rna,pch=key$Shape,bg=key$Color,cex=1.5)
