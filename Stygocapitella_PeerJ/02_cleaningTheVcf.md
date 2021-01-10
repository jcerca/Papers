On 01_ we cleaned our .vcf for -r 0.5 and -p 8. This is not enough since it still violates the assumptions of some analysis.

For example, for PCA we need the data to be independent, and we should add coverage filters to remove repeated regions and regions with low coverage (we have low confidence on these variants).


We start by adding a minimum allele frequency, max depth and minimum depth cut-off:
```
vcftools --gzvcf Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.vcf.gz \
--recode --stdout --maf 0.05 --max-meanDP 100 --min-meanDP 10 > \
Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.vcf

# I know! The file name is huge. But it allows us knowing what filters we have:
# -r 0.5 -p 4 (from before
# minimum allele frequency of 0.05, max depth of 100 and minimum depth of 10.
```

Now, we analyse missing data.

```
# First, we use vcftools to get a grasp of % of missing data.
vcftools --missing-indv --vcf Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.vcf

# We can sort out the individuals with more than 0.90 missing data like using the output from the command above. Like this:
awk '$5>0.90' out.imiss | awk '{print $1}' > inds.with.missignessabove90.tsv

# Now, we remove these individuals from the vcf.. and give it even a worse name!
vcftools --vcf ../1_mafmaxmeanminmean/blue* \
--exclude inds.with.missingnessabove90.tsv \
--recode --stdout > \
Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.indswith90missingnessRemoved.vcf
```

OK. Final cleaning we need one random SNP per radseq locus. There are easier ways to do this (e.g. using the --write-single-snp flag in the module "populations"), but it wasn't working in the version of stacks I was using.

```
# First, we separate the header of the VCF (lines starting with #), and keep only a variant per chromosome. We save this on a file called "randomized.noheader.vcf"

VCF=Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.indswith90missingnessRemoved.vcf
grep -v "^#" $VCF | sort -R  | awk ' {if (count[$1] < 1) {count[$1]++; print $0}}' > randomized.noheader.vcf

# Now, we get the header isolated:
grep "^#" $VCF > header.tsv

# and get them together, using an even more appalling name:

cat header.tsv randomized.noheader.vcf > \
Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.indswith90missingnessRemoved.randomSNP.vcf
```
