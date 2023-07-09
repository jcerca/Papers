### 01 - We start by re-running ANGSD
```
angsd -GL 2 -bam cleaned.bam.list -out ./01_angsdRun/spider_ngsDistRun \
-minInd 38 -minQ 30 -minMapQ 30 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05 -doCounts 1 -doGlf 2 -doMajorMinor 1 \
-nThreads 20 -doGeno 8 -doPost 1 -geno_minDepth 7 -geno_maxDepth 30
```

### 02 - We processed it by removing repeat regions, as done before.
```
# For this analysis we need the geno file, but the snps are separated by tabs and not by underscores, so we change that.
zcat ../../01_angsdRun/spider_ngsDistRun.geno.gz | sed "s/\t/_/" > tmp.geno

grep -F -f snpsNotonRepeats.tsv tmp.geno > spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno

# We want to target the second "_". A easy way is to change the two "_", and then change the first:
sed "s/_/\t/g" spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno | sed "s/\t/_/" > spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno.tmp

# We revert the file
mv spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno.tmp spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno
gzip spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno

## Getting a sample list
cp ../../sample.list .

```
### 03 - and run NGSdist
```
03_ngsDIST
# With ngsDist we can compute pairwise genetic distances without relying on individual genotype calls.
# To get number of sites
zcat spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno.gz | wc -l
#16101029

/data/bigexpansion/jcerca/local_bin/ngsDist/ngsDist -verbose 1 -geno ../02_genoProcessing/spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno.gz \
-probs -n_ind 76 -n_sites 16101029 -labels ../02_genoProcessing/sample.list \
-o spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.dist -n_threads 40
```
### 04 - We run fastme

```
04_fastme
## This tree is without bootstrap replicates
fastme -D 1 -i ../03_ngsDIST/spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.dist -o fastme.tree -m b -n b
## Now for replicates:
# We will generate 20 replicates by randomly sampling with replacemente blocks of 20 SNPs
# (since we called SNPs earlier, otherwise will indicate the genomic length).
05_bootstrappedReplicates_ngsDIST
cd ../02*/

ngsDist -verbose 1 -geno spider.FastaRun.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.geno.gz -probs -n_ind 76 -n_sites 16101029 -labels sample.list -o spider.FastaRun.MinCoverage7.MaxCoverage30.HWE_LTRcleaned.repeatsRemoved.renamed.BOOTSTRAP.dist -n_threads 40 -n_boot_rep 100 -boot_block_size 20
fastme -D 101 -i spider.FastaRun.MinCoverage7.MaxCoverage30.HWE_LTRcleaned.repeatsRemoved.renamed.BOOTSTRAP.dist -o fastree_bootstrap.tree -m b -n b

```
### 05 - and place support on the trees with raxml.
```
head -n 1 fastree_bootstrap.tree > ../06_placing_the_support_onTrees/mainTree.tree
tail -n +2 fastree_bootstrap.tree > ../06_placing_the_support_onTrees/allotherTrees.tree

conda activate phyl
raxmlHPC -f b -t mainTree.tree -z allotherTrees.tree -m GTRCAT -n spider_ngsTree_bootstrap
```
