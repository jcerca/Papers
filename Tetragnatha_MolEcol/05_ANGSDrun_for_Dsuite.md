### 01 - angsd RUN

```
angsd -GL 2 -bam ../bam.list -out spiderFastaRun -minInd 38 -minQ 30 -minMapQ 30 -SNP_pval 1e-6 -doMaf 1 -minMaf 0.05 -doGlf 2 -doMajorMinor 1 -nThreads 20 -doGeno 4 -doPlink 2 -doPost 2 -doCounts 1 -geno_minDepth 7 -geno_maxDepth 30

```

### 02 - converting to vcf (before I cleaned it by removing repeats, as shown in previous code
```
####
plink --tfile spider.fastaRun.repeatFiltered --double-id --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex --make-bed --recode vcf  --out spider.fastaRun.repeatFiltered

cat ../../../07_NGSdist/cleaned.bam.list | awk -F "/" '{print $NF}' | sed 's/_dedup_mapQual_sorted.bam//' > sample.list

bcftools reheader -s sample.list spider.fastaRun.HWEfiltered.repeatFiltered.vcf > renamed.vcf
mv renamed.vcf spider.fastaRun.HWEfiltered.repeatFiltered.vcf
```

### 03 - D suite!
Following https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data

```
Dsuite Dtrios ../02_finalVCF/spider.fastaRun.repeatFiltered.vcf popmap.tsv -o trios

# Same as above, made a tree:
(Outgroup,(T_pilosa,((T_mohihi,T_kauaiensis),((T_obscura_kukuiki_kikokiko_anuenue,T_quasimodo),(T_tantalus_polychromata_brevignatha_waikamoi,(T_perreirai,(T_restricta,T_kamakou))))));


### Exploring on R ... all D- and F-stats are below 0.2
setwd("~/Desktop/tmp/tetrag/d_stats")
setwd("~/Desktop/tmp/tetrag/d_stats")

D_BBAA <- read.table("./trios_BBAA.txt",as.is=T,header=T)
plot(D_BBAA$Dstatistic, ylab="D",xlab="trio number")
D_BBAA[which(D_BBAA$Dstatistic > 0.1),]
plot(D_BBAA$p.value, ylab="p value",xlab="trio number",ylim=c(0,0.05))
plot(p.adjust(D_BBAA$p.value,method="BH"), ylab="p value",xlab="trio number",ylim=c(0,0.05))
plot(D_BBAA$f4.ratio, ylab="f4-ratio",xlab="trio number", ylim=c(0,1))
######

### Now getting some nice svg plots.
cut -f 2 popmap.tsv | sort | uniq > plot_order.txt
# plot_d.rb came from
https://github.com/millanek/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/plot_d.rb
# plot_f4ratio.rb came from
https://github.com/millanek/tutorials/blob/master/analysis_of_introgression_with_snp_data/src/plot_f4ratio.rb

## Note, I deleted the *Dmin.txt file without wanting to do so...

#Plotting f4 and Dstat
# 0.16 is the max D value - I saw this in the R code above.
ruby plot_d.rb trios_BBAA.txt plot_order.txt 0.16 ABBA_BABA_BBAA.svg
# 0.11 is the max f4 value
ruby plot_f4ratio.rb trios_BBAA.txt plot_order.txt 0.11 f4.svg

# Now doing the F branch
../../../../local_bin/Dsuite/Build/Dsuite Fbranch tree.nwk trios_tree.txt > Fbranch.txt
../../../../local_bin/Dsuite/utils/dtools.py Fbranch.txt tree.nwk

```
03 - D suite investigate

```
#  we need a popmap-like file.
popmap.tsv

# We also need a set tests file. They have it formated as:
# mbuna deep  Diplotaxodon

# On the tutorial they state:
# "The “TestTrios” file specifies that we want to investigate the admixture signal between the Diplotaxodon genus and the deep benthic group". So, we have to add:
# T_restricta\tT_perreirai\tT_kamakou
# did this as
nano test.sets

# Windows of 50 informative snps, moving forward 25 snps at a time.
../../../../local_bin/Dsuite/Build/Dsuite Dinvestigate -w 50,25 spider.fastaRun.repeatFiltered.vcf.gz popmap.tsv test.sets


### 252 seconds
T_restricta T_perreirai T_kamakou
D=-0.234851
f_d=-0.147053 -26721.7/181714
f_dM=-0.0761522 -26721.7/350899


# Windows of 50 informative snps, moving forward 5 snps at a time.
../../../../local_bin/Dsuite/Build/Dsuite Dinvestigate -w 50,5 spider.fastaRun.repeatFiltered.vcf.gz popmap.tsv test.sets
T_restricta T_perreirai T_kamakou
D=-0.234851
f_d=-0.147053 -26721.7/181714
f_dM=-0.0761522 -26721.7/350899

```
