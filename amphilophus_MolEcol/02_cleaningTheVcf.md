Let's get a high-quality dataset!

First, let's clean for coverage & minimum allele frequency. This is a RADseq dataset with wide variation in coverage, so we will be somehow permissive: We want coverages between 15-200 and a minimum allele frequency of 0.01 (notice, we have hundreds of individuals, so even if this value sounds very reduced, it should target mostly errors).
We also include a max-missing 0.25 filter to filter for loci with tons of missing data.
```
vcftools -vcf amphilophus_snp_filter_hardpass.vcf.gz \
        --maf 0.01 \
        --maxDP 200 \
        --max-meanDP 200 \
        --minDP 15 \
        --min-meanDP 15 \
        --max-missing 0.25 \
        --recode \
        --stdout > amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.vcf.gz
```

Now, we want to generate a Linkage Disequilibrium-prunned dataset for admixture and PCA analysis - two analyses where te variants should be independent. We do this using plink. First plink run is to obtain the LD-information (.ld file)

```
plink --vcf amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.vcf \
--recode --allow-extra-chr --const-fid 0 --allow-no-sex --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --out amphilophus_ldDecay
```

Now we have some LD-stats to plot. Let's go to R.

```
rm(list = ls())

setwd("C:/Users/josece_adm/Desktop/AMPHILOPHUS/ld_decay/newdataSet_20190427/")

library(readr)
library(dplyr)
library(getopt)

# specify command line options
spec <- matrix(c(
  'infile', 'i', 1, 'character', 'specify path to input files',
  'outfile', 'o', 1, 'character', 'specify path to output file'), ncol = 5, byrow = T)

# set command line options
opt = getopt(spec)

# show help if asked for
if (!is.null(opt$help)) {
  cat(paste(getopt(spec, usage=T),"\n"));
  q();
}

# set variables for call
infile <- opt$infile
outfile <- opt$outfile

# dummy variables for testing
infile <- ("./amphilophus_ldDecay.ld")
outfile <- ("./out_LDdecay.csv")

# This .ld file was created by the plink command above :)

# read in data
myData <- read_table(infile, skip = 1, col_names = F)[, c(-3, -6)]
# rename
colnames(myData) <- c("chrA", "bpA", "chrB", "bpB", "R2")

# create a diff
myData <- mutate(myData, dist = abs(bpA-bpB))

# create breaks, labels etc
breaks <- seq(0, signif(max(myData$dist), digits = 1), 1000)
labels <- seq(0, signif(max(myData$dist), digits = 1)-1000, 1000)
# bin here
myData$dist_bin <- cut(myData$dist, breaks = breaks, labels = labels)

# mean R2 per bin
bin_data <- myData %>% group_by(dist_bin) %>% summarise(mean_R2 = mean(R2))

write.csv(bin_data, outfile, quote = F, row.names = F)

library(ggplot2)
ggplot(bin_data, aes(dist_bin, mean_R2)) + geom_point() + geom_smooth(formula = newfile$file.dist(bp) ~ newfile$newr2, se=F)
ggplot(bin_data, aes(dist_bin, mean_R2, group = 1)) +
  geom_point() +
  geom_line()
```

Decision time! We opted to prune for a r2 > 0.1 (r2 is a measure of LD). Let's use plink to extract it.

```
plink --vcf amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex --indep-pairwise 50 10 0.1 --out amphilophus_REF_pruning_step1

#extracting the actual data..
plink --vcf 1_MAF_depth_coverage_missing/amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.vcf \
--double-id --allow-extra-chr --set-missing-var-ids @:# --extract amphilophus_REF_pruning_step1.prune.in --make-bed --out amphi_ref_ld_prune --recode vcf
```

OK. We're almost there. Let's extract individuals wth tons of missing data (more than 50% data missing)
```
# Getting missing data per individual:
vcftools --vcf amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.LDprunned.vcf --missing-indv

# Getting individuals with more than 50% missing data
cat out.imiss | awk ' > 0.5'
cat out.imiss | awk '$5 > 0.5 ' | awk '{print $1}' > toremove.tsv
vcftools --vcf amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.LDprunned.vcf \
--remove toremove.tsv --recode --recode-INFO-all --stdout \
> amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.LDprunned.vcf.IndsWith0.5MissingnessRemoved.vcf```


Important notice: These steps were repeated for different datasets. For instance, we did not include a minimum allele cut-off for the FastSimCoal analysis, and we did not include a LD filter for FST. See our paper for a pipeline :)
