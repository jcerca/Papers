Ideograms - Jose Cerca 18/05/2021
```
# Following  https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html
setwd("~/Desktop/tmp/distmat/")
#install.packages("RIdeogram")
require(RIdeogram)
library(tidyverse)
library(scales)
show_col(viridis_pal()(3))

gene_density <- GFFex(input = "./annotation.subgenomeA.gff3", karyotype ="./chr_ids.subgenomeA.txt", feature = "gene", window = 2500000)
karyotype<-read.table("./chr_ids.subgenomeA.txt", sep="\t", header=T, stringsAsFactors = F)

# ideogram(karyotype = karyotype)
# convertSVG("chromosome.svg", device = "png")
# ideogram(karyotype = karyotype, overlaid = gene_density)
# convertSVG("chromosome.svg", device = "png")


TE_density <- GFFex(input = "./repeats.subgenomeA.gff3", karyotype ="./chr_ids.subgenomeA.txt", feature = "similarity", window = 2500000)
TE_density$Color<-"fc8d62"

LTR_density <- GFFex(input = "./LTRs.subgenomeA.gff3", karyotype ="./chr_ids.subgenomeA.txt", feature = "repeat_region", window = 2500000)
LTR_density$Color<-"8da0cb"

### Merging TE_density and LTR_density - basically we need to change the columns to Value_1 and Color_1 and then add the new info as Value_2, Color_2
# names(TE_density) <-c("Chr", "Start", "End", "Value_1", "Color_1")
# TE_density$Value_2<-LTR_density$Value
# TE_density$Color_2<-LTR_density$Color
# head(TE_density)

ideogram(karyotype = karyotype, overlaid = gene_density, label = TE_density, label_type = "polygon")
convertSVG("chromosome.svg", device = "png")


### SubgenomeB

gene_density <- GFFex(input = "./annotation.subgenomeB.gff3", karyotype ="./chr_ids.subgenomeB.txt", feature = "gene", window = 2500000)
karyotype<-read.table("./chr_ids.subgenomeB.txt", sep="\t", header=T, stringsAsFactors = F)

TE_density <- GFFex(input = "./repeats.subgenomeB.gff3", karyotype ="./chr_ids.subgenomeB.txt", feature = "similarity", window = 2500000)
TE_density$Color<-"fc8d62"

ideogram(karyotype = karyotype, overlaid = gene_density, label = TE_density, label_type = "polygon")
convertSVG("chromosome.svg", device = "png")
```
