### 01 - Data
```
# I got the same vcf as for D suite.
ln -s ../../00_vcf/02_renaming/spider.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.vcf .

# We also need the popmap used in D suite
popmap.tsv

# I ran D suite investigate as in the analysis 05, and established the following comparisons:

Maroons comparison - T_restricta_T_kamakou_T_perreirai_localFstats_maroon_50_25.txt
Small Brown comparison 1 - T_tantalus_polychromata_brevignatha_waikamoi_T_restricta_T_mohihi_localFstats_test_smallB_comparison1_50_25.txt
Small Brown comparison 2 - T_tantalus_polychromata_brevignatha_waikamoi_T_obscura_kukuiki_kikokiko_anuenue_T_mohihi_localFstats_test_smallB_comparison2_50_25.txt
Small Brown comparison 3 - T_tantalus_polychromata_brevignatha_waikamoi_T_restricta_T_obscura_kukuiki_kikokiko_anuenue_localFstats_test_smallB_comparison3_50_25.txt

# I'll only focus on Maroon from now on for this github page.
# Windows of 100 informative snps, moving forward 50 snps at a time.
Dsuite Dinvestigate -w 50,25 ../00_vcf/spider.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.vcf popmap.tsv ../01_popmap/test_maroons.set -n maroon
```

### 03 - R plotting and getting top 1% windows
```
awk '{print $6, "maroon"}' ../02_investigate/T_restricta_T_kamakou_T_perreirai_localFstats_maroon_50_25.txt > f_dM_maroon.tsv

# top 1% D_suite
# How many windows we have?
sort -k 6 T_restricta_T_kamakou_T_perreirai_localFstats_maroon_50_25.txt | wc -l
# 9306 total

# Extracting the top 93 (94 because of the header)
sort -k 6 T_restricta_T_kamakou_T_perreirai_localFstats_maroon_50_25.txt | tail -n 94 > ../04_top_1percent_f_dM_windows/f_dM_maroon_topWindows.tsv
```

### 05 - DXY
```
# 01 - converting to geno
/data/bigexpansion/jcerca/local_bin/genomics_general/VCF_processing/parseVCF.py -i ../../00_vcf/02_renaming/spider.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.vcf -o ./spider.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno

# 02 - maroons
window=10000
python /data/bigexpansion/jcerca/local_bin/genomics_general/popgenWindows.py -w $window \
-g ../../01_parsing/spider.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno -o maroon_maroon_green${window} -T 10 -f phased \
-p maroon_1 T_kamakou_008,T_kamakou_009,T_kamakou_010,T_kamakou_011,T_kamakou_012,T_kamakou_013,T_kamakou_014,T_kamakou_015,T_kamakou_068 \
-p maroon_2 T_perreirai_061,T_perreirai_024,T_perreirai_025 \
-p green T_kauaiensis_016,T_kauaiensis_017,T_kauaiensis_018,T_kauaiensis_019,T_kauaiensis_020

```
### 06 - Getting the windows. Whole genome comparison
```
# Getting the values for all DXY comparisons (whole genome)
head -n 1 ../../../05_dxy/02_maroons/maroon_maroon_green10000
scaffold,start,end,mid,sites,pi_maroon_1,pi_maroon_2,pi_green,dxy_maroon_1_maroon_2,dxy_maroon_1_green,dxy_maroon_2_green,Fst_maroon_1_maroon_2,Fst_maroon_1_green,Fst_maroon_2_green

# Here we are interested in column number 9 - maroon vs maroon.
awk -F ',' '{print $9}' ../../../05_dxy/02_maroons/maroon_maroon_green10000 | grep -v "nan" > DXY_maroons_all_windows.tsv

awk '{print $0 "\t01_maroon_vs_maroon_GenomeWide_DXY"}' DXY_maroons_all_windows.tsv > tmp
mv tmp DXY_maroons_all_windows.tsv
```

### 07 - Getting the DXY for the top 1% fDM regions
```
# I will obtain the average DXY on the top 1% regions with f_dM

head -n 1 ../../../05_dxy/02_maroons/maroon_maroon_green10000
scaffold,start,end,mid,sites,pi_maroon_1,pi_maroon_2,pi_green,dxy_maroon_1_maroon_2,dxy_maroon_1_green,dxy_maroon_2_green,Fst_maroon_1_maroon_2,Fst_maroon_1_green,Fst_maroon_2_green
# Again, as in the case above, we are comparing maroon vs maroon for the genome genome and top 1%. Here, we want column nr 9.

# For each region I ran a code similar to:
awk -F "," '$1=="ScpjasG_3236" && $2>1015476 && $3< 1079884' ../../../05_dxy/02_maroons/maroon_maroon_green10000  | awk -F ',' '{ total += $9 } END { print total/NR }'  >> DXY_maroons_top1percent_windows.tsv

## code explained:
# awk -F "," '$1=="ScpjasG_3236" && $2>1015476 && $3< 1079884' ../05_dxy/02_maroons/maroon_maroon_green10000 ### First, column 1 has to be the scaffold name, second column has to be bigger than the starting point of the region and the latter column smaller than the end of the region. I got the region numbers of D suite investigate output.
# awk -F ',' '{ total += $9 } END { print total/NR }' ### Above we established we wanted column 9. Here we get it
#  >> DXY_maroons_top1percent_windows.tsv ### Saving it to a new file.

# Then we make it R-friendly ...
awk '{print $0 "\t02_Top_1%_Fdm_Maroon_Maroon_DXY"}' DXY_maroons_top1percent_windows.tsv > tmp
mv tmp DXY_maroons_top1percent_windows.tsv


# I did this for all comparisons established and concatenated a new file
cat DXY_maroons_top1percent_windows.tsv DXY_maroons_all_windows.tsv

```

### 08 - plotting
```
setwd("/Users/josepc/Library/CloudStorage/OneDrive-NTNU/tmp/dsuiteTet/maroons")
library(tidyverse)
library(see)
library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)


dxy<-read_table(file="./fullData.tsv", col_names=FALSE)
colnames(dxy)<-c("dxy", "comparison")
dxy$dxy<-as.numeric(dxy$dxy)


plt <- ggbetweenstats(
  data = dxy,
  x = comparison,
  y = dxy
)

plt

# Does it look good? Then plot the 'reduced version' and save it.


dxy<-read_table(file="./reduced.tsv", col_names=FALSE)
colnames(dxy)<-c("dxy", "comparison")
dxy$dxy<-as.numeric(dxy$dxy)


plot <- ggbetweenstats(
  data = dxy,
  x = comparison,
  y = dxy
)

plot

```
### 09 - Wilcoxon test

```
# 01

# Wilcoxon test to compare means
library(tidyverse)
library(rstatix)
library(ggpubr)
library(DescTools)

setwd("~/OneDrive - Universitetet i Oslo/tmp/wilcox/")

df<-read_table("01_fullData.tsv", col_names = c("value", "level"))

## First comparison :
# 01_maroon_vs_maroon_GenomeWide_DXY
# 02_Top_1_Fdm_Maroon_Maroon_DXY

variable1<-subset(df, df$level==c("01_maroon_vs_maroon_GenomeWide_DXY"))
variable2<-subset(df, df$level==c("02_Top_1_Fdm_Maroon_Maroon_DXY"))

var1<-as.data.frame(variable1$value) %>%
  `colnames<-`("01_maroon_vs_maroon_GenomeWide_DXY")

var2<-as.data.frame(variable2$value) %>%
  `colnames<-`("02_Top_1_Fdm_Maroon_Maroon_DXY")

nrow(var1)
nrow(var2)

varX<-var1[sample(1:91), ]
varY<-var2[1:91,]

final_df <- data.frame(var1 = rep(NA, 91), var2 = rep(NA, 91))

final_df$maroon_vs_maroon_GenomeWide_DXY<-varX
final_df$Top_1_Fdm_Maroon_Maroon_DXY<-varY

# The test assumes that differences between paired samples should be distributed symmetrically around the median.
# So we compute the differences between pairs and create histograms:

long <- final_df %>%
  gather(key = "group", value = "weight", maroon_vs_maroon_GenomeWide_DXY, Top_1_Fdm_Maroon_Maroon_DXY)
head(long, 3)

bxp <- ggpaired(long, x = "group", y = "weight",
                ylab = "Value", xlab = "Level")
bxp

final_df <- final_df %>% mutate(differences = maroon_vs_maroon_GenomeWide_DXY - Top_1_Fdm_Maroon_Maroon_DXY)
gghistogram(final_df, x = "differences", y = "..density..",
            fill = "steelblue",bins = 5, add_density = TRUE)

stat.test <- long  %>%
  wilcox_test(weight ~ group, paired = TRUE) %>%
  add_significance()

stat.test
```
