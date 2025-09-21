First we concatenate all vcfs (we have a vcf per chr):

```
vcf_list=$(ls *gz | grep -v "SUPER_[XY]"  | grep -v "Scaffold" | grep -v "MT")

bcftools concat --threads 8 -n -O z -o reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75.vcf.gz ${vcf_list}
bcftools index reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75.vcf.gz
```

Then we explore linkage disequilibrium:

```
#We set the variable
vcf=reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75.vcf.gz

# Running plink
plink --vcf $vcf --recode --allow-extra-chr --const-fid 0 --allow-no-sex --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --out decay

# gzip, as it is required by Mark's script
echo "gz'ippin'"
gzip decay.ld

# Now running Mark's script
echo "Time for Mark's script"
# ld_decay.py - is in this folder
/cluster/projects/nn9244k/jose_cerca/029_reindeer/02_popgen/01_PCA/02_LDfilter/ld_decay.py -i decay.ld.gz -o ld_reindeer
```

We explore the files obtained using R:

```
library(tidyverse)

setwd("~/Downloads")

myData <- read.table(file="./ld_reindeer.ld_decay_bins", header = T, sep="\t")

# bin here

# mean R2 per bin
bin_data <- myData %>% group_by(distance) %>% summarise(R2 = mean(avg_R2))

colnames(myData) <- c("R2","dist")

ggplot(bin_data, aes(distance, R2)) +
  geom_point() +
  xlim(0, 25000) +
  geom_hline(yintercept=0.1075)


## Decision = 0.1075

```

And use this value to clean for ld:

```
#Prunning
ml PLINK/1.9b_6.13-x86_64
plink --vcf ../01_concatenate_vcf/reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75.vcf.gz --double-id  --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex --indep-pairwise 25 5 0.1075 --out reindeer_processing_tmp

plink --vcf ../01_concatenate_vcf/reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75.vcf.gz --double-id  --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex --extract reindeer_processing_tmp.prune.in --make-bed --recode vcf --out reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75_LD0.1075

#
ml BCFtools/1.19-GCC-13.2.0
bgzip *vcf
bcftools index reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75_LD0.1075.vcf.gz
```

For the PCA we use plink:

```
plink --vcf ../03_Ldprunning/reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75_LD0.1075.vcf.gz --double-id  --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex --pca --out pca_reindeer
```

And plot results using R:

```
setwd("~/Desktop/tmp/pca_reindeer/")

df <- read.table("./pca_reindeer_v2.eigenvec")
new_column_names <- c("population", "individual", paste0("PC", 1:20))
colnames(df) <- new_column_names

library(RcppCNPy)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggfortify)
library(randomcoloR)
library(plotly)


#### morphs #####
PC1_2 <- ggplot(data=df, aes(x=PC1, y=PC2, color=population)) +
  geom_point() +
  xlab(paste("PC1 (3.14204%)")) +
  ylab(paste("PC2 (2.34353%)")) +
  coord_equal() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme_bw()

PC1_3 <- ggplot(data=df, aes(x=PC1, y=PC3, color=population)) +
  geom_point() +
  xlab(paste("PC1 (3.14204%)")) +
  ylab(paste("PC3 (1.8628%)")) +
  coord_equal() +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

require(gridExtra)
grid.arrange(PC1_2, PC1_3, ncol=2)
#### morphs #####
PC1_2
PC1_3
```

To do the phylogeny we first convert the files using Edgardo M. Ortiz's vcf2phyl.py (https://github.com/edgardomortiz/vcf2phylip) ...
```
ml Python/3.12.3-GCCcore-13.3.0
./vcf2phyl.py -i ../01_data/reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75_LD0.1075.vcf.gz -n
```

... and run the phylogeny:
```
nexus=../02_convert/reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75_LD0.1075.min4.nexus

ml IQ-TREE/2.2.2.7-gompi-2023a
iqtree2 -st DNA -m GTR+ASC --threads-max 16 -B 1000 -s $nexus
```

For admixture analysis we:

```
## 01 - convert the file using plink
plink --vcf ../../03_Ldprunning/reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75_LD0.1075.vcf.gz --recode12 --out reindeer_admixture --allow-extra-chr --double-id --make-bed

# Remember that admixture wants integers as chrs - so we change "SUPER_" to chr
sed -i "s/SUPER_/chr/g" reindeer_admixture.bim
sed -i "s/chr//g" reindeer_admixture.bim

## 02 - We run admixture:

ml ADMIXTURE/1.3.0

## Variables
# k and iteration are exported
input=/cluster/projects/nn9244k/jose_cerca/029_reindeer/02_popgen/01_PCA_phyl_admix/07_admixture/01_plink/reindeer_admixture.ped

##
mkdir -p ./K${k}
cd ./K${k}
mkdir -p ./run_${iteration}
cd ./run_${iteration}
admixture --seed=$RANDOM --cv $input $k -j8 | tee log


# To calculate fit and values:

grep  "^CV" */*log
grep -A 1 "^Loglikelihood" */*log

# K1 Highest log likelihood is run_1
run_1/log:Loglikelihood: -106713288.803101
run_1/log-CV error (K=1): 0.39666

# K2 Highest log likelihood is run_5
run_5/log:Loglikelihood: -99459579.385929
run_5/log:CV error (K=2): 0.39916

# K3 Highest log likelihood is run_4
run_4/log:Loglikelihood: -96218323.764672
run_4/log:CV error (K=3): 0.43688

# K4 Highest log likelihood is run_1
run_1/log:Loglikelihood: -94056033.798305
run_1/log:CV error (K=4): 0.47673

# K5 Highest log likelihood is run_2
run_2/log:Loglikelihood: -92264209.764925
run_2/log:CV error (K=5): 0.54281

# K6 Highest log likelihood is run_3
run_3/log:Loglikelihood: -90222305.122066
run_3/log:CV error (K=6): 0.57596

# K7 Highest log likelihood is run_2
run_2/log:Loglikelihood: -88659088.943749
run_2/log:CV error (K=7): 0.57570
```

Now we plot Admixture results:

```
## examining Admixture data ##
rm(list = ls())

library(tidyverse)
library(gridExtra)

setwd("~/Desktop/tmp/admix_reindeer//")

# first read in and examine cross-validation error
cv_error<- tbl_df(read.table("./cv_error.tsv", sep="\t", header=T))

# read in sampleividual data
#sample<-as.character(read.table("./admixture.input.fam")[,1])
myData <- tbl_df(read.table("./popmap.tsv", header=T))# %>%
  #rename(pop=pop, sample=sample_ID)
#The sample.pop.species was created with a custom awk script based on the .fam file

# plot
a <- ggplot(cv_error, aes(K, error))+ geom_point() + geom_line()
a + ylab("Cross-valudation error") + theme_light()

## K2
# set up the dataset
K2 <- tbl_df(read.table("./reindeer_admixture.2.Q",
                        col.names = c("one", "two")))
# add sampleividuals
K2 <- mutate(K2, sample = myData$sample, pop = myData$pop)
# gather to plot
K2 <- gather(K2, key = cluster, value = q, -sample, -pop)


## K3
# set up the dataset
K3 <- tbl_df(read.table("./reindeer_admixture.3.Q",
                        col.names = c("one", "two", "three")))
# add sampleividuals
K3 <- mutate(K3, sample = myData$sample, pop = myData$pop)
# gather to plot
K3 <- gather(K3, key = cluster, value = q, -sample, -pop)


## K4
# set up the dataset
K4 <- tbl_df(read.table("./reindeer_admixture.4.Q",
                        col.names = c("one", "two", "three", "four")))

# add sampleividuals
K4 <- mutate(K4, sample = myData$sample, pop = myData$pop)
# gather to plot
K4 <- gather(K4, key = cluster, value = q, -sample, -pop)


## K5
# set up the dataset
K5 <- tbl_df(read.table("./reindeer_admixture.5.Q",
                        col.names = c("one", "two", "three", "four", "five")))
# add sampleividuals
K5 <- mutate(K5, sample = myData$sample, pop = myData$pop)
# gather to plot
K5 <- gather(K5, key = cluster, value = q, -sample, -pop)


##K6
# set up the dataset
K6 <- tbl_df(read.table("./reindeer_admixture.6.Q",
                        col.names = c("one", "two", "three", "four", "five","six")))
# add sampleividuals
K6 <- mutate(K6, sample = myData$sample, pop = myData$pop)
# gather to plot
K6 <- gather(K6, key = cluster, value = q, -sample, -pop)

##K7
# set up the dataset
K7 <- tbl_df(read.table("./reindeer_admixture.7.Q",
                        col.names = c("one", "two", "three", "four", "five","six","seven")))
# add sampleividuals
K7 <- mutate(K7, sample = myData$sample, pop = myData$pop)
# gather to plot
K7 <- gather(K7, key = cluster, value = q, -sample, -pop)


# set colour palettes
k2_cols <- c("dodgerblue", "firebrick1")
k3_cols <- c("dodgerblue", "gold",  "firebrick1")
k4_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen")
k5_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen", "darkorchid")
k6_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen", "darkorchid", "grey66")
k7_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen", "darkorchid", "grey66", "#8B4513")
k8_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen", "darkorchid", "grey66", "#8B4513", "#FFB6C1")
k9_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen", "darkorchid", "grey66", "#8B4513", "#FFB6C1", "#AFEEEE")
k10_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen", "darkorchid", "grey66", "#8B4513", "#FFB6C1", "black")


#by pops!
# make plot
a <- ggplot(K2, aes(as.character(sample), q, fill = cluster)) + geom_bar(stat = "identity")
a <- a + scale_fill_manual(values = k2_cols)
a <- a + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
a <- a + theme_light() + theme(legend.position = "none",
              panel.border = element_blank(),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 10),
              axis.text.x = element_text(size = 6, angle = 90),
              axis.ticks.x = element_blank(),
              strip.text.x = element_text(colour = "black", face = "bold"),
              strip.background = element_rect(colour = "black", fill = "white"))
a


b <- ggplot(K3, aes(as.character(sample), q, fill = cluster)) + geom_bar(stat = "identity")
b <- b + scale_fill_manual(values = k3_cols)
b <- b + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
b <- b + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
b


c <- ggplot(K7, aes(as.character(sample), q, fill = cluster)) + geom_bar(stat = "identity")
c <- c + scale_fill_manual(values = k7_cols)
c <- c + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
c <- c + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
c

grid.arrange(a,b,c,nrow=3)



###
##### with facet_wrap
a <- ggplot(K2, aes(as.character(sample), q, fill = cluster)) + geom_bar(stat = "identity")
a <- a + scale_fill_manual(values = k2_cols)
a <- a + facet_wrap(~pop, scales = "free_x", nrow=1)
a <- a + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
a <- a + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
a


b <- ggplot(K3, aes(as.character(sample), q, fill = cluster)) + geom_bar(stat = "identity")
b <- b + scale_fill_manual(values = k3_cols)
b <- b + facet_wrap(~pop, scales = "free_x", nrow=1)
b <- b + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
b <- b + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
b


c <- ggplot(K7, aes(as.character(sample), q, fill = cluster)) + geom_bar(stat = "identity")
c <- c + scale_fill_manual(values = k7_cols)
c <- c + facet_wrap(~pop, scales = "free_x", nrow=1)
c <- c + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
c <- c + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))

grid.arrange(a,b,c,nrow=3)
```
