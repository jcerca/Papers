Okidoki - so I started with a VCF file generated from GATK, and we will be looking at lane, coverage effects using PCAs and sample information.

We start by getting the heterozygosity, missing data per individual, and depth using vcftools
```
vcftools --vcf $VCF --het --out heterozygosity.amphilophus
vcftools --vcf $VCF --missing-indv --out missing.amphilophus
vcftools --vcf $VCF --depth --out depth.amphilophus
```
This created a .het, imiss, and idepth files.

Before we move to R, I formated an excel file that looks like this:
Individual;Library;Sex;Species;Location;Morphotype.
This will allow us to put some of these variables on the PC and explore them together with depth (i.e. does depth vary with population?), library (i.e. are there library artefacts on the species level?)

```
head amph_sample_data_160218.csv
ind;library;sex;species;location;citrinellus_morph
13NIej18_U_C_LV;3;U;Citrinellus;LV;Green
13NIej15_U_C_LV;1;U;Citrinellus;LV;Green
13NIej27_U_C_LV;5;U;Citrinellus;LV;Green
13NIej30_U_C_LV;1;U;Citrinellus;LV;Green
13NIej34_U_C_LV;2;U;Citrinellus;LV;Green
13NIej29_F_C_LV;2;F;Citrinellus;LV;Green
13NIej20_U_C_LV;1;U;Citrinellus;LV;Green
13NIej16_U_C_LV;5;U;Citrinellus;LV;Green
13NIej22_U_C_LV;1;U;Citrinellus;LV;Green
```

Also, below we read the "PCscores.csv" - a CSV file with PCscores. Check 03_PCA to see how it's done. In any case, the file looks like this:
Species, site, PC1, PC2, PC3, PC4
```
head PCscores.csv
Species,Site,PC1,PC2,PC3,PC4
C,Om,-0.21842003651071404,-2.5389094851659704,-2.1716174142783196,-0.8750975875669117
C,Om,10.466172387978698,2.4278327240982436,1.6859045090934124,-0.2507413954093064
C,Om,-1.1052059971730377,-2.613541232545509,-0.20623229806478713,-1.7704471183535233!
```

Now, let's do R plotting:
```
rm(list = ls())

library(tidyverse)
setwd("C:/Users/josece_adm/Desktop/AMPHILOPHUS/DATA_exploration/")

het <- read_delim("./heterozygosity.amphilophus.het", delim = "\t", col_names = c("ind", "ho", "he", "n_sites", "f"), skip = 1)
miss <- read_delim("./missing.amphilophus.imiss", delim = "\t", col_names = c("ind", "n_data", "n_fil", "n_miss", "f_miss"), skip = 1)
depth <- read_delim("./depth.amphilophus.idepth", delim = "\t", col_names = c("ind", "n_sites", "mean_depth"), skip = 1)
extra <- read_csv2("./amph_sample_data_160218.csv")
pc <- read_csv("./PCscores.csv")

# First, we combine different datasets (heterozygosity and missing)
myData <- left_join(het, miss, by = "ind")
# Now we add the depth information
myData$depth <- depth$mean_depth
# and we add the additional information
myData <- left_join(myData, extra, by = "ind")

# For this particular case we will be altering individual names
myData$ind <- unlist(map(strsplit(myData$ind, "_"), 1))

# And we join datasets again!
myData <- left_join(myData, select(pc, PC1, PC2, ind), by = "ind")

# We get library as factor:
myData <- mutate(myData, library = as.factor(library))

# ... and plot
ggplot(myData, aes(f, f_miss, colour = species)) + geom_point()
ggplot(myData, aes(f, f_miss, colour = library)) + geom_point()

ggplot(myData, aes(depth, f_miss, colour = species)) + geom_point()
ggplot(myData, aes(depth, f_miss, colour = library)) + geom_point()

ggplot(myData, aes(depth, f, colour = species)) + geom_point()

ggplot(myData, aes(depth, PC1, colour = species)) + geom_point()
ggplot(myData, aes(depth, PC2, colour = species)) + geom_point()

ggplot(myData, aes(f, PC1, colour = species)) + geom_point()
ggplot(myData, aes(f, PC2, colour = species)) + geom_point()

ggplot(myData, aes(f_miss, PC1, colour = species)) + geom_point()
ggplot(myData, aes(f_miss, PC2, colour = species)) + geom_point()

ggplot(myData, aes(PC1, PC2, colour = library)) + geom_point()
ggplot(myData, aes(PC1, PC2, colour = f_miss)) + geom_point()
ggplot(myData, aes(PC1, PC2, colour = f)) + geom_point()

# investigate outlier data
filter(myData, PC2 < -10)
```

Alles gut!
