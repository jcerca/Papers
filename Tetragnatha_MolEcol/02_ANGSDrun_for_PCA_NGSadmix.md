### 01 - let's do this! Let's get ANGSD runing

```
angsd -GL 2 -bam ../bam.list -out spider -minInd 38 -minQ 30 -minMapQ 30 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05 -doGlf 2 -doMajorMinor 1 -nThreads 40 -doGeno -4 -doPlink 2 -doPost 2 -doCounts 1
```

### 02 - Now we remove the variants (Genotype likelihoods) from repeat regions
```
# We use the whiteList1 crested in the (01 part of this code)
conda activate python2
python2 ReadBEDFilterBEAGLEv01.py bed_input=whiteList1.bed beagle_input=spider.beagle beagle_output=spider.repeatsRemoved.beagle

#
zcat spider.repeatsRemoved.beagle.gz | awk 'FNR > 1 {print $1}' > loci_after.repeats.tsv
grep -w -F --file loci_after.repeats.tsv ../00_ANGSD/spider.tped > spider.repeats.tped

# -w tells grep to match whole words only (i.e. so ABC123 won't also match ABC1234).
# -F search for fixed strings (plain text) rather than regular expressions
# -f genelist.txt read search patterns from the file

# and get the tfam file
cp ../00*/*tfam ./spider.repeats.tfam
```

### 03 - We pruned the data using PLINK, 0.11 linkage (see my paper on Amphilophus on how I calculate LD using plink)
```
#Prunning
plink --tfile ../01_MikesScript_RepeatRegionsRemoved/spider.repeats --double-id  --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex --indep-pairwise 25 5 0.11 --out spider
plink --tfile ../01_MikesScript_RepeatRegionsRemoved/spider.repeats --double-id  --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex --extract spider.prune.in --make-bed --recode vcf --out spider.repeats_cleaned.ldPrunned0.11
```

### 04 - We run NGS admix to get admixture plots
```
# Doing 1 by one:
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
K=5
mkdir -p K_$K
NGSadmix -likes spider.repeatsCleaned.LDprunned0.11.OUTGROUPremoved.bgl.gz -K $K -P 2 -o ./K_$K/K_${K}_run_1 -seed $RANDOM &

# and plot it
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(randomcoloR)
library(tidyverse)
library(dplyr)

# mfrow: fills in by row
# mfcol: fills in by column
# par(mfcol = c(2, 2))

n <- 16
col <- distinctColorPalette(n)

myData <- read_table("/Users/josepc/Library/CloudStorage/OneDrive-NTNU/tmp/structure_tet/ind_pop.tsv", col_names=F) %>%
  rename("ind"="X1", "species"="X2")

## K2
# set up the dataset
K2 <-tibble(read.table("/Users/josepc/Library/CloudStorage/OneDrive-NTNU/tmp/structure_tet/K_2.qopt",
                        col.names = c("one", "two")))
# add individuals
K2 <- mutate(K2, ind = myData$ind, species = myData$species)
# gather to plot
K2 <- gather(K2, key = cluster, value = q, -ind, -species)

# set up the dataset
K3<- tibble(read.table("/Users/josepc/Library/CloudStorage/OneDrive-NTNU/tmp/structure_tet/K_3.qopt",
                       col.names = c("one", "two", "three")))
# add individuals
K3 <- mutate(K3, ind = myData$ind, species = myData$species)
# gather to plot
K3 <- gather(K3, key = cluster, value = q, -ind, -species)

# set up the dataset
K15<- tibble(read.table("/Users/josepc/Library/CloudStorage/OneDrive-NTNU/tmp/structure_tet/K_15.qopt",
                        col.names = c("one", "two", "three",
                                      "four","five", "six",
                                      "seven", "eight", "nine",
                                      "ten", "eleven", "twelve",
                                      "thirteen", "forteen", "fifteen")))
# add individuals
K15 <- mutate(K15, ind = myData$ind, species = myData$species)
# gather to plot
K15 <- gather(K15, key = cluster, value = q, -ind, -species)


# k2_cols <- c("dodgerblue", "firebrick1")


#adding a lake column..
a <- ggplot(K2, aes(as.character(ind), q, fill = cluster)) + geom_bar(stat = "identity")
a <- a + scale_fill_manual(values = col)
a <- a + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0)
)
a <- a + facet_wrap(~species, scales = "free_x", nrow = 1)
a <- a + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white")
)
a




b <- ggplot(K3, aes(as.character(ind), q, fill = cluster)) + geom_bar(stat = "identity")
b <- b + scale_fill_manual(values = col)
b <- b + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0)
)
b <- b + facet_wrap(~species, scales = "free_x", nrow = 1)
b <- b + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white")
)
b

c <- ggplot(K15, aes(as.character(ind), q, fill = cluster)) + geom_bar(stat = "identity")
c <- c + scale_fill_manual(values = col)
c <- c + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0)
)
c <- c + facet_wrap(~species, scales = "free_x", nrow = 1)
c <- c + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white")
)
c

grid.arrange(a,b,c, nrow=3)

```

### 05 - and the PCA
```
pcangsd --beagle spider.repeatsCleaned.LDprunned0.115.bgl.gz -e 6 --threads 30 -o pca_6
```
