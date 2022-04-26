01_genomes
```
# all fasta alignments (bams)
```

02_psmcConvert
```
REF=../../01_S_atractyloidesGenomeDATA/00_refgenome/scalesia_atractyloides.orderedByChrSize.fasta
for i in bam; do bcftools mpileup -a DP -Ou -Q 30 -q 30 -f $REF $i |  bcftools call -c | vcfutils.pl vcf2fq -d 5 > ../02_psmcConvert/${i%bam}fq

for i in *bam; do echo $i; bcftools mpileup --threads 2 -a DP -Ou -Q 30 -q 30 -f $REF $i |  bcftools call -c --threads 2 | vcfutils.pl vcf2fq -d 5 > ../02_psmcConvert/${i%bam}fq; done
# This one is better to separate and do 10-20 loops with 2 individuals each

for i in *fq; do echo $i; ../../../../local_bin/psmc/utils/fq2psmcfa $i > ${i%fq}psmc; done
```

03_psmcRuns

```
for i in ../02*/*psmc; do j=$(echo ${i#../02_psmcConvert/}); echo $j;
../../../../local_bin/psmc/utils/splitfa $i > ./${j}fa; mkdir ${j};
seq 100 | xargs -P 10 -i ../../../../local_bin/psmc/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o ./${j}/diploid_round-{}.psmc ${j}fa;
cd $j; cat *psmc > combined.$j.psmc; cd ..; done
```

04_psmcPlots
```
ln -s ../03_psmcRuns/*/combined* .
for i in combined*; do ../../../../local_bin/psmc/utils/psmc_plot.pl -g 3 -Y 50 -X 20000000 -u 6e-09 ${i%.psmc.psmc} $i; done
```
