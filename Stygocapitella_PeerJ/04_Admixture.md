We've completed a PCA and MDS plot, now we do an Admixture analysis.

First, we need to convert the vcf.
```
VCF=Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.indswith90missingnessRemoved.randomSNP.vcf

# Chromossomes cannot be numbers.. so we add "c" in front of all of them
sed "s/^/^c/g" $VCF > tmp.vcf

# We remove the "c" from the header, i.e. lines starting with the cardinal (#).
nano tmp.vcf 

# We convert tmp.vcf using plink
plink --vcf tmp.vcf --recode12 --out admixture.input --allow-extra-chr --double-id --make-bed

# and now remove the c (This sounds circular, but welcome to genomics)
sed -i "s/^c//g" admixture.input.bim

```

Now, we run admixture.

```
input=admixture.input.bed
for i in $(seq 2 6); do
K=${i}
echo "This is my K for today = $i"
cd K$i
mkdir run1 run2 run3 run4 run5
cd run1
admixture --cv $INPUT $K -j4 | tee log
cd ..
cd run2
admixture --cv $INPUT $K -j4 | tee log
cd ..
cd run3
admixture --cv $INPUT $K -j4 | tee log
cd ..
cd run4
admixture --cv $INPUT $K -j4 | tee log
cd ..
cd run5
admixture --cv $INPUT $K -j4 | tee log
cd ..
cd ..
done
```

And finally we need to plot, but we should format some files first:

```
# First, we create a file consisting of:
# <indivdual><tab><species><population>
# In this dataset it looks like this..
awk '{print $1}' admixture.input.fam | awk -F "_" '{print $(NF-3)"_"$(NF-2) "\t"  $(NF-1) "\t" $1 }' > inds.pop.species.tsv

# Then, we get the CV error.
cat ./K*/run*/log | grep "CV error" > cv.error.tsv

# and format a final file with:
# <K><tab><error_value>
# <2><tab><error_forK2>
# <3><tab><error_forK3>
# ...
# We format by using nano:
nano cverror.final
```

Now, for the great finalle, we plot. For the plots we need the cverror.final file we have created, as well as the files terminating in ".Q" produced by the admixture analysis.

We move to R:
```
## examining Admixture data ##
rm(list = ls())

library(tidyverse)
library(gridExtra)

setwd("C:/Users/Cerca/Desktop/ongoing_stuff_ambientedetrabalho/project3/17_admixture/03_plot")

# we need the ind.species.pop.tsv file which is a "awk-modified" file from plink.fam
# the .Q files (they contain admixture proportions)
# the cv_error file

# first read in and examine cross-validation error
cv_error<- tbl_df(read.table("./cverror.final", sep="\t", header=T))

# read in individual data
# ind<-as.character(read.table("./admixture.input.fam")[,1])
myData <- tbl_df(read.table("./inds.pop.species.tsv", col.names = c("ind","species","site")))
#The ind.pop.species was created with a custom awk script based on the .fam file

# plot
a <- ggplot(cv_error, aes(K, error))+ geom_point() + geom_line()
a + ylab("Cross-valudation error") + theme_light()

## K2
# set up the dataset
K2 <- tbl_df(read.table("./admixture.input.2.Q",
                        col.names = c("one", "two")))
# add individuals
K2 <- mutate(K2, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K2 <- gather(K2, key = cluster, value = q, -ind, -site, -species)


## K3
# set up the dataset
K3 <- tbl_df(read.table("./admixture.input.3.Q",
                        col.names = c("one", "two", "three")))

# add individuals
K3 <- mutate(K3, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K3 <- gather(K3, key = cluster, value = q, -ind, -site, -species)


## K4
# set up the dataset
K4 <- tbl_df(read.table("./admixture.input.4.Q",
                        col.names = c("one", "two", "three", "four")))

# add individuals
K4 <- mutate(K4, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K4 <- gather(K4, key = cluster, value = q, -ind, -site, -species)

## K5
# set up the dataset
K5 <- tbl_df(read.table("./admixture.input.5.Q",
                        col.names = c("one", "two", "three", "four", "five")))
# add individuals
K5 <- mutate(K5, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K5 <- gather(K5, key = cluster, value = q, -ind, -site, -species)


##K6
# set up the dataset
K6 <- tbl_df(read.table("./admixture.input.6.Q",
                        col.names = c("one", "two", "three", "four", "five","six")))
# add individuals
K6 <- mutate(K6, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K6 <- gather(K6, key = cluster, value = q, -ind, -site, -species)


# set colour palettes
k2_cols <- c("dodgerblue", "firebrick1")
#k3_cols <- c("dodgerblue", "forestgreen",  "firebrick1")
k4_cols <- c("dodgerblue", "forestgreen",  "firebrick1", "gold")
k5_cols <- c("dodgerblue", "forestgreen",  "firebrick1", "gold", "darkorchid")
k6_cols <- c("dodgerblue", "forestgreen",  "firebrick1", "gold", "darkorchid", "grey66")

k3_cols <- c("#9ACD31",
             "#1D90FF",
             "#FF8001")

#by sites!
# make plot
a <- ggplot(K2, aes(as.character(ind), q, fill = cluster)) + geom_bar(stat = "identity")
a <- a + scale_fill_manual(values = k2_cols)
a <- a + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
a <- a + facet_wrap(~site, scales = "free_x")
a <- a + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
a

a

# make plot
b <- ggplot(K3, aes(ind, q, fill = cluster)) + geom_bar(stat = "identity")
b <- b + scale_fill_manual(values = k3_cols)
b <- b + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
b <- b + facet_wrap(~site, scales = "free_x")
b <- b + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
b

grid.arrange(a + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()),
             b + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()))


######################################by species
# make plot
c <- ggplot(K2, aes(as.character(ind), q, fill = cluster)) + geom_bar(stat = "identity")
c <- c + scale_fill_manual(values = k2_cols)
c <- c + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
c <- c + facet_wrap(~species, scales = "free_x")
c <- c + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
c

# make plot
d <- ggplot(K3, aes(ind, q, fill = cluster)) + geom_bar(stat = "identity")
d <- d + scale_fill_manual(values = k3_cols)
d <- d + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
d <- d + facet_wrap(~species, scales = "free_x")
d <- d + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
                               
d

grid.arrange(b + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()),
             d + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()))
```
