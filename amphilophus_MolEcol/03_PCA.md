It's PC time!

We start with a vcf file, and do the PCA on R directly in this case.

```
#Setwd
setwd("C:/Users/Cerca/Desktop/ongoing_stuff_ambientedetrabalho/AMPHILOPHUS/PCA")

#libraries
library(vcfR)       #Data loading
library(adegenet)   #Data analysis
library(tidyverse)  #PCA
library(gridExtra)  #PCA
library(ape)        #Phylogenetic analyses
library(lemon)      #for the g_legend
#loading_vcf
vcf_file<-read.vcfR("amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.LDprunned.vcf.IndsWith0.5MissingnessRemoved.vcf")
str(vcf_file)

#converting_vcf_to_genlight
genlight_file<-vcfR2genlight(vcf_file)


# Plotting missing data
glPlot(genlight_file, posi="topleft")

###### Neighbourjoining TREE
tre <- nj(dist(as.matrix(genlight_file)))
tre
plot(tre, typ="fan", cex=0.7)
title("NJ tree Amphilophus")


#  PC time
pca1 <- glPca(genlight_file)
4 #axis
# Looks fishy

scatter(pca1, posi="bottomright")
title("PC_Amphilophus\n axes 1-2")


# Making it pretty, as suggested by Mark Ravinet
library(plyr); library(dplyr)

#CONVERTING THE PCA'S RESULTS TO A DATA_FRAME AND MANUALLY MAKING THE PLOT.

# The code below modifies the individual names using the information encoded and is dataset specific.
temp<-rownames_to_column(as.data.frame(pca1$scores),"rownames")                               # converting the row names to the first column
temp$LAB_ID<-temp$rownames                                                                    # copying the first column to last
temp_v2 <- temp[,-1]                                                                          # creating a new dataframe without the first column
rownames(temp_v2)<-temp[,1]                                                                   # assigning the first column of one to the other
temp_v3 <- data.frame(do.call('rbind', strsplit(as.character(temp_v2$LAB_ID),'_',fixed=TRUE))) # create a new dataframe with the splitted rows
temp_v4 <- cbind(temp_v2[,1-4],temp_v3[,1-4])                                                 # joining both dataframes
temp_v5<-rename(temp_v4, c("X7"="Species","X4"="Site"))                                       # renaming


# Beautiful PC
myPCA<-data.frame(temp_v5$Species, temp_v5$Site, pca1$scores[, 1:4], temp_v5$LAB_ID, row.names = NULL) %>%
  rename(c("temp_v5.Species"="Species","temp_v5.Site"="Site"))

write_csv(myPCA, "./PCscores.csv")
#myPCA<-read.table("./PCscores.csv",sep=",", header=T)

# calculate PVE
pve <- (pca1$eig/sum(pca1$eig))*100

# make PCA scatterplot
a <- ggplot(myPCA, aes(PC1, PC2, colour = Site, pch = Species)) + geom_point(size=3)
b <- ggplot(myPCA, aes(PC2, PC3, colour = Site, pch = Species)) + geom_point(size=3)
# add PVE labels
a <- a + xlab(paste0("PC1 (", signif(pve[1], 2), "%)")) + ylab(paste0("PC2 (", signif(pve[2], 2), "%)"))
b <- b + xlab(paste0("PC2 (", signif(pve[2], 2), "%)")) + ylab(paste0("PC3 (", signif(pve[3], 2), "%)"))
# extract legend
c <- g_legend(a + theme(legend.position = "right"))
# make plot
grid.arrange(arrangeGrob(a + theme(legend.position = "none"),
                         b + theme(legend.position = "none"), ncol = 2),
             c, ncol = 2, widths = c(6, 1))
```
