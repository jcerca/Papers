Alright, for an FST analysis we need a list of individuals:

```
VCF=Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.indswith90missingnessRemoved.randomSNP.vcf

bcftools query -l $VCF | grep "blue" > S_subterranea.tsv
bcftools query -l $VCF | grep "green" > S_josemariobrancoi.tsv
bcftools query -l $VCF | grep "purple" > S_westheidei.tsv
```

And for the analysis:
```
vcftools --vcf --fst-window-size 50000 --fst-window-step 10000 --weir-fst-pop S_subterranea.tsv --weir-fst-pop S_josemariobrancoi.tsv --out Ss_Sj.tsv
vcftools --vcf --fst-window-size 50000 --fst-window-step 10000 --weir-fst-pop S_subterranea.tsv --weir-fst-pop S_westheidei.tsv --out Ss_Sw.tsv
vcftools --vcf --fst-window-size 50000 --fst-window-step 10000 --weir-fst-pop S_josemariobrancoi.tsv --weir-fst-pop S_westheidei.tsv --out Sj_Sw.tsv
```
Simple, eh? Happy this doesn't look anything like the fastsimcoal run!
