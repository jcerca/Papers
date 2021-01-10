A principal components analysis should be one of your first analysis! Here, we will be using the .vcf file generated in step 02.

We will be doing this in R.

```
# Few notes. The native format for adegenet, a library to generate PCAs, is "genlight".
# Useful tutorial from Filip Kolar which goes from vcf to genlight directly.
# Filip's tutorial  https://botany.natur.cuni.cz/hodnocenidat/Lesson_05_tutorial.pdf

#Clean up!
rm(list=ls())

#Set Working directory
setwd("C:/Users/Cerca/Desktop/ongoing_stuff_ambientedetrabalho/project3/02_PCA/")

#Load up those Libraries
library(vcfR)       #Data loading
library(adegenet)   #Data analysis
library(tidyverse)  #PCA
library(gridExtra)  #PCA
library(ape)        #Phylogenetic analyses
library(lemon)      #for the g_legend

# Loading up the VCF
vcf_file<-read.vcfR("./Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.indswith90missingnessRemoved.randomSNP.vcf")
str(vcf_file)

# Converting from vcf to genlight
genlight_file<-vcfR2genlight(vcf_file)


#Loci missingness
glPlot(genlight_file, posi="topleft")


# PCA
pca1 <- glPca(genlight_file)
4 #axis

# Now we need to convert the samples into a data_frame to plot.

temp<-rownames_to_column(as.data.frame(pca1$scores),"rownames") #passing the row names to the first column
temp$LAB_ID<-temp$rownames #copying the first column to last
temp$pop<-word(temp$LAB_ID,1,sep="_")
temp$clade<-word(temp$LAB_ID,-2,sep="_")

# calculating the % of variance explained
pve <- (pca1$eig/sum(pca1$eig))*100

a <- ggplot(temp, aes(PC1, PC2, colour = pop, pch=clade)) + geom_point() +
  xlab(paste0("PC1 (", signif(pve[1], 2), "%)")) + ylab(paste0("PC2 (", signif(pve[2], 2), "%)"))
a
b <- ggplot(temp, aes(PC1, PC3, colour = pop, pch=clade)) + geom_point() +
  xlab(paste0("PC1 (", signif(pve[1], 2), "%)")) + ylab(paste0("PC3 (", signif(pve[3], 2), "%)"))
b

c <- ggplot(temp, aes(PC2, PC3, colour = pop, pch=clade)) + geom_point() +
  xlab(paste0("PC2 (", signif(pve[2], 2), "%)")) + ylab(paste0("PC3 (", signif(pve[3], 2), "%)"))
c

d <- ggplot(temp, aes(PC1, PC4, colour = pop, pch=clade)) + geom_point() +
  xlab(paste0("PC1 (", signif(pve[1], 2), "%)")) + ylab(paste0("PC4 (", signif(pve[4], 2), "%)"))
d

pdf("pca.pdf", width = 11, height = 8.5, pointsize = 10, family = "Helvetica")
a
b
c
d
dev.off()
```
