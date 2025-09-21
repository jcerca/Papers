First we concatenate all vcfs (we have a vcf per chr):

```
vcf_list=$(ls *gz | grep -v "SUPER_[XY]"  | grep -v "Scaffold" | grep -v "MT")

bcftools concat --threads 8 -n -O z -o reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75.vcf.gz ${vcf_list}
bcftools index reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75.vcf.gz
```

Then we prune for linkage disequilibrium:

```
#We set the variable
vcf=reindeer_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75.vcf.gz

# Running plink
plink --vcf $vcf --recode --allow-extra-chr --const-fid 0 --allow-no-sex --r2 --ld-window 100 --ld-window-kb 100 --ld-window-r2 0 --out decay

# gzip, as it is required by Mark's script
echo "gz'ippin'"
gzip decay.ld

# Now running Mark's script
echo "Time for Mark's script"
/cluster/projects/nn9244k/jose_cerca/029_reindeer/02_popgen/01_PCA/02_LDfilter/ld_decay.py -i decay.ld.gz -o ld_reindeer
```
