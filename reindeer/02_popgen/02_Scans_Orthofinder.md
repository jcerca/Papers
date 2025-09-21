To get FST we started by converting the files to geno:

```
#Unfortunately I forgot to save the script. But here is the example from Simon Martin's github:
python VCF_processing/parseVCF.py -i input.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 max=50 -o output.geno.gz
```

Then we calculate FST:

```
# input= input folder
# output= geno file

python ~/local_bin/genomics_general/popgenWindows.py -w 100000 -m 25000 \
-g $input/$geno -o $output -T 10 -f phased \
-p Forollhogna F01,F02,F03,F04,F05,F07,F08,F09,F10,F11,F12,F13 \
-p wild_Knutshoe_Rondane_Snoehetta ST01,ST02,ST03,ST04,ST05,ST06,ST07,R01,R02,R03,R04,R05,R06,R07,KT01,KT02,KT03,KT04,KT05,KT06,KT07 \
-p Riast_Hyllingen RH01,RH02,RH03,RH04,RH05,RH06,RH07,RH08,RH09,RH10,RH11,RH12 \
```

To run orthofinder we did:
```
ml OrthoFinder/2.5.5-foss-2023a
orthofinder -f . -t 40
```

To phase the data, we used SHAPEIT5:
```
SHAPEIT5_phase_common --input ${vcf} --region ${window} --output ${window}_phased.bcf --thread 12

# Now we sort files:
#Sorting files
for super_file in $(ls SUPER_*_phased.bcf); do
    # Extract the SUPER_X prefix (e.g., SUPER_1) from each filename
    super_prefix=$(echo $super_file | cut -d':' -f1)

    # Gather all files with the same SUPER_X prefix and sort them numerically based on the number after ':'
    ls ${super_prefix}:*_phased.bcf | sort -t':' -k2,2n > "${super_prefix}.tsv"
done

SHAPEIT5_ligate --input ${SUPER}.tsv --output ${SUPER}_snps_qual20_mindp15_maxdp75_phased.bcf --index --thread 12
bcftools view -O z -o ${SUPER}_snps_qual20_mindp15_maxdp75_phased.vcf.gz ${SUPER}_snps_qual20_mindp15_maxdp75_phased.bcf

```

To run iHS and xp-EEH we used a Rscript:
```
rm(list = ls())

# install.packages(c("getopt", "tidyverse", "rehh", "gridExtra"))
library(tidyverse)
library(rehh)
library(getopt)
library(gridExtra)


# Specify command line options
spec=matrix(c(
  'window_start', 's', 1, "double", 'specify the start of the window',
  'window_end', 'e', 1, "double", 'specify the end of the window',
  'chr', 'c', 1, "character", 'specify chromosome',
  'gene', 'g', 1, "character", 'specify gene name'
),byrow=T, ncol=5)

# set command line options
opt = getopt(spec)

#gene<-"FUNCG00000001178"
#chr<-"SUPER_1"
#window_start<-as.numeric("80767308")
#window_end<-as.numeric("80768753")


# set variables for call
window_start <- opt$window_start
window_end <- opt$window_end
chr <- opt$chr
gene <-opt$gene


extended_start<-as.numeric(window_start)-10000000
extended_end<-as.numeric(window_end)+10000000


wild_ihs<-read.csv("wild_ihs.csv")
feral_ihs<-read.csv("feral_ihs.csv")
semidomestic_ihs<-read.csv("semidomestic_ihs.csv")

wild_semidomestic <- read.csv("xpehh_wild_semidomestic.csv")
semidomestic_feral <- read.csv("xpehh_semidomestic_feral.csv")



# plotting iHS results - whole chr
a <- ggplot(wild_ihs$ihs, aes(wild_ihs$ihs.POSITION, wild_ihs$ihs.IHS)) + geom_point(color="#87CEFA", alpha = 0.5) + theme_classic() +
 ggtitle("Wild") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5) + theme_classic()

b <- ggplot(feral_ihs$ihs, aes(feral_ihs$ihs.POSITION, feral_ihs$ihs.IHS)) + geom_point(color="#87CEFA", alpha = 0.5) + theme_classic() +
  ggtitle("Feral") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5) + theme_classic()

c <- ggplot(semidomestic_ihs$ihs, aes(semidomestic_ihs$ihs.POSITION, semidomestic_ihs$ihs.IHS)) + geom_point(color="#87CEFA", alpha = 0.5) + theme_classic() +
 ggtitle("Semidomestic") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5) + theme_classic()


require(gridExtra)


jpeg(file = paste0(gene,  "_iHS_wholechr_", chr, "_", window_start, "_", window_end, ".jpg"))
grid.arrange(a,b,c, ncol=1)
dev.off()

# plotting iHS p-value results - whole chr
d <- ggplot(wild_ihs$ihs, aes(wild_ihs$ihs.POSITION, wild_ihs$ihs.LOGPVALUE)) + geom_point(color="#87CEFA", alpha = 0.5) + theme_classic() +
  ggtitle("Wild") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5) + theme_classic()

e <- ggplot(feral_ihs$ihs, aes(feral_ihs$ihs.POSITION, feral_ihs$ihs.LOGPVALUE)) + geom_point(color="#87CEFA", alpha = 0.5) + theme_classic() +
  ggtitle("Feral") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5) + theme_classic()

f <- ggplot(semidomestic_ihs$ihs, aes(semidomestic_ihs$ihs.POSITION, semidomestic_ihs$ihs.LOGPVALUE)) + geom_point(color="#87CEFA", alpha = 0.5) + theme_classic() +
  ggtitle("West") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5) + theme_classic()

jpeg(file = paste0(gene, "_iHS_whole_chr_pvalue_", chr, "_", window_start, "_", window_end, ".jpg"))
grid.arrange(d,e,f, ncol=1)
dev.off()


# plotting iHS results - zooming
g <- ggplot(wild_ihs$ihs, aes(wild_ihs$ihs.POSITION, wild_ihs$ihs.IHS)) + geom_point(color="#87CEFA", alpha = 0.5) + theme_classic() +
  ggtitle("Wild") +
  xlim(extended_start, extended_end) +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5) + theme_classic()

h <- ggplot(feral_ihs$ihs, aes(feral_ihs$ihs.POSITION, feral_ihs$ihs.IHS)) + geom_point(color="#87CEFA", alpha = 0.5) + theme_classic() +
  ggtitle("Feral") +
  xlim(extended_start, extended_end) +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5) + theme_classic()

i <- ggplot(semidomestic_ihs$ihs, aes(semidomestic_ihs$ihs.POSITION, semidomestic_ihs$ihs.IHS)) + geom_point(color="#87CEFA", alpha = 0.5) + theme_classic() +
  ggtitle("West") +
  xlim(extended_start, extended_end) +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5) + theme_classic()

jpeg(file = paste0(gene, "_iHS_10mb_window_", chr, "_", window_start, "_", window_end, ".jpg"))
grid.arrange(g,h,i, ncol=1)
dev.off()


### xpEHH now

# xpEEH whole chr
j<-ggplot(wild_semidomestic, aes(POSITION, XPEHH_wild_semidomestic)) +
  theme_bw() +
  ggtitle("Wild vs Semidomestic") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5)

l<-ggplot(semidomestic_feral, aes(POSITION, XPEHH_semidomestic_feral)) +
  theme_bw() +
  ggtitle("Semidomestic vs Feral") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5)


jpeg(file = paste0(gene, "_xpEEH_", chr, "_", window_start, "_", window_end, ".jpg"))
grid.arrange(j, l, ncol=1)
dev.off()

# xpEEH whole chr p-values
m<-ggplot(wild_semidomestic, aes(POSITION, LOGPVALUE)) +
  theme_bw() +
  ggtitle("Wild vs Semidomestic") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5)

n<-ggplot(semidomestic_feral, aes(POSITION, LOGPVALUE)) +
  theme_bw() +
  ggtitle("Semidomestic vs Feral") +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5)


jpeg(file = paste0(gene, "_xpEEH_pvalue_", chr, "_", window_start, "_", window_end, ".jpg"))
grid.arrange(m, n, ncol=1)
dev.off()

# plotting xpEEH results - zooming
o<-ggplot(wild_semidomestic, aes(POSITION, XPEHH_wild_semidomestic)) +
  theme_bw() +
  ggtitle("Wild vs Semidomestic") +
  xlim (extended_start, extended_end) +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5)

p<-ggplot(semidomestic_feral, aes(POSITION, XPEHH_semidomestic_feral)) +
  theme_bw() +
  ggtitle("Semidomestic vs Feral") +
  xlim (extended_start, extended_end) +
  geom_vline(xintercept = c(window_start,window_end), color = "black", size=0.5) +
  geom_point(color="#87CEFA", alpha = 0.5)


jpeg(file = paste0(gene, "_xpEEH_", chr, "_", window_start, "_", window_end, ".jpg"))
grid.arrange(o, p, ncol=1)
dev.off()
```
