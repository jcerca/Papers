We now mark duplicates:

```
# ref= reference genome
# bam_file= alignment files
# TMP=`pwd`/tmp # Temporary folder

# first we sort using Picard
gatk SortSam \
-I ${bam_file} \
-O ${sample}_sort.bam \
-SO coordinate --TMP_DIR $TMP

# then we mark duplicates
echo "### Running Picard MarkDuplicates on $sample ###"
gatk MarkDuplicates \
-I ${sample}_sort.bam \
-O ${sample}_sort_dedup.bam \
--METRICS_FILE ${sample}_dedup_metrics.txt --TMP_DIR $TMP

# Now we index bams using PICARD
echo "### Running Picard sort on $sample ###"
gatk BuildBamIndex \
-I ${sample}_sort_dedup.bam --TMP_DIR $TMP
```

For variant calling we started by dividing the genome in 10mb blocks:

```
cut -f 1-2 mRanTar2.1.hap1.fa.fai > genome_size.txt

# Now we break the genome in 10 mb windows
ml BEDTools/2.31.0-GCC-12.3.0
bedtools makewindows -g genome_size.txt -w 10000000 | sed "s/\t/:/; s/\t/-/" > 10mb.tsv
```

Now we call variants:
```
## set variables ##
# reference_genome= reference genome
# bams= list with location of bamfiles
# region = region of 10mb along the genome

# So we will be having each chr runing on its own folder
chr=$(echo $region | sed "s/:.*//")
mkdir -p /cluster/work/users/josece/var_call/${chr}
cd /cluster/work/users/josece/var_call/${chr}

echo "Doing ${region}"

bcftools mpileup --threads 4 -d 8000 --ignore-RG -r ${region} -a AD,DP,SP -Ou -f ${reference_genome} -b ${bams} | \
bcftools call --threads 4 -f GQ,GP -mO z -o ${region}.vcf.gz

touch ${region}_complete
```

Now we will concatenate the different 10mb into chromosomes:

```
# chr = we define the chromosome
cd /cluster/work/users/josece/var_call/${chr}

echo "Doing ${chr}"

# In a folder with all the 10mb-vcfs for one chromsome
vcf_list=$(ls *gz | sort -t"-" -k2 -n)
bcftools concat --threads 4 -n -O z -o ${chr}_concat.vcf.gz ${vcf_list}
bcftools index ${chr}_concat.vcf.gz
```

Now we rename the vcfs:

```
bcftools query -l ${chr}_concat.vcf.gz | xargs -n 1 basename | awk -F '_' '{print $1}' > samples_${chr}
bcftools reheader -s samples_${chr} -o ${chr}.vcf.gz ${chr}_concat.vcf.gz
```

We now explore filters by:
```
# 1 - subsetting the vcf (it is more pratical like this):

bcftools view ${chr}.vcf.gz > ${chr}.vcf
ml vcflib/1.0.9-foss-2022a-R-4.2.1
vcfrandomsample -r 0.001 ${chr}.vcf > ${chr}_subset.vcf
touch ${chr}_subset_complete

# 2 - We index the subsets
bgzip ${chr}_subset.vcf
# index vcf
bcftools index ${chr}_subset.vcf.gz

# 3 - Concatenate the subsets
vcf_list=$(ls *gz | sort -t"-" -k2 -n)

bcftools concat --threads 4 -n -O z -o subset.vcf.gz ${vcf_list}
bcftools index subset.vcf.gz

# 4 - Explore filters:
echo "Calculate allele frequency"
vcftools --gzvcf ${vcf} --freq2 --out ${out} --max-alleles 2

echo "Calculate mean depth per individual"
vcftools --gzvcf ${vcf} --depth --out ${out}

echo "Calculate mean depth per site"
vcftools --gzvcf ${vcf} --site-mean-depth --out ${out}

echo "Calculate site quality"
vcftools --gzvcf ${vcf} --site-quality --out ${out}

echo "Calculate proportion of missing data per individual"
vcftools --gzvcf ${vcf} --missing-indv --out ${out}

echo "Calculate proportion of missing data per site"
vcftools --gzvcf ${vcf} --missing-site --out ${out}

echo "Calculate heterozygosity and inbreeding coefficient per individual"
vcftools --gzvcf ${vcf} --het --out ${out}

```

To plot and explore these filters, we use R:

```
# load tidyverse package
library(tidyverse)

setwd("/Users/josecer/Library/CloudStorage/OneDrive-UniversitetetiOslo/Science/Projects/2022_UiO_CEES/03_Projects/01_Reindeer_Atle/06_Bioinformatics/01_snp_quality_filters")

#Variant quality
var_qual <- read_delim("./reindeer.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)

a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0,25)

# will try minimum threshold of 30

#Depth
var_depth <- read_delim("./reindeer.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0,60)

summary(var_depth$mean_depth)

# Missing data
var_miss <- read_delim("./reindeer.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)

# Minor allele frequency
var_freq <- read_delim("./reindeer.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

var_freq <- var_freq %>%
  mutate(
    # Split the "a1" column at the tab "\t" and extract the second value or assign "0" if not present
    a2 = ifelse(grepl("\t", a1), sub(".*\t", "", a1), "0"),
    # Remove the tab character from "a1" (keeping only the first part)
    a1 = sub("\t.*", "", a1)
  )

var_freq$maf <- as.numeric(var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z)))

a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim (0.01,0.2)
summary(var_freq$maf)


# indv depth
ind_depth <- read_delim("./reindeer.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0,70)


# indv iss
ind_miss  <- read_delim("./reindeer.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 0.025)


# heterozygosity
ind_het <- read_delim("./reindeer.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
```

We normalize the vcfs:

```
# Normalizing files
bcftools view -V indels -e 'ALT="*" | N_ALT>1' ${chr}.vcf.gz | bcftools norm -D -O z -o ${chr}_norm.vcf.gz
```

To get an 'all-sites vcf' we do:

```
echo "Working on $chr"
cd /cluster/work/users/josece/var_call/norm

# Moving these to a new folkder, and indexing them
mkdir -p /cluster/work/users/josece/var_call/all_sites_vcfs

qual=20
min_dp=15
max_dp=75

vcftools --gzvcf ${chr}_norm.vcf.gz --remove-indels --remove-filtered-all \
 --max-alleles 2 \
 --minQ $qual \
 --min-meanDP ${min_dp} --max-meanDP ${max_dp} \
 --minDP ${min_dp} --maxDP ${max_dp} \
 --recode --recode-INFO-all --stdout | \
 bcftools view -e 'N_ALT>1' -O z -o  ${chr}_filtered_allsites_qual20_mindp15_maxdp75.vcf.gz

mv ${chr}_filtered_allsites_qual20_mindp15_maxdp75.vcf.gz /cluster/work/users/josece/var_call/all_sites_vcfs
cd /cluster/work/users/josece/var_call/all_sites_vcfs

bcftools index ${chr}_filtered_allsites_qual20_mindp15_maxdp75.vcf.gz
```

To get a snps-only vcf, we do:

```
qual=20
min_dp=15
max_dp=75

vcftools --gzvcf ${chr}_norm.vcf.gz --remove-indels --remove-filtered-all \
 --min-alleles 2 --max-alleles 2 \
 --minQ $qual \
 --min-meanDP ${min_dp} --max-meanDP ${max_dp} \
 --minDP ${min_dp} --maxDP ${max_dp} \
 --recode --recode-INFO-all --stdout | \
 bcftools view -e 'N_ALT>1' -O z -o ${chr}_snps_qual20_mindp15_maxdp75.vcf.gz

mv ${chr}_snps_qual20_mindp15_maxdp75.vcf.gz /cluster/work/users/josece/var_call/snp_vcfs
cd /cluster/work/users/josece/var_call/snp_vcfs

bcftools index ${chr}_snps_qual20_mindp15_maxdp75.vcf.gz
```
