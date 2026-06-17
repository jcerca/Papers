We split the genome into windows

```
# Getting chromossome sizes
cut -f 1-2 /cluster/projects/nn9244k/jose_cerca/029_reindeer/01_dataprocessing/04_alignment/ref_genome/mRanTar2.1.hap1.fa.fai > genome_size.txt

# Breaking appart the genome in 10 mb windows
ml BEDTools/2.31.0-GCC-12.3.0
bedtools makewindows -g genome_size.txt -w 10000000 | sed "s/\t/:/; s/\t/-/" > 10mb.tsv
```

Variant calling
```
reference_genome=
bams=

# So each chr runs on its own folder - for chr
# chr=$(echo $region | sed "s/:.*//")
# mkdir -p /cluster/work/users/josece/var_call/${chr}
# cd /cluster/work/users/josece/var_call/${chr}

# So all scaffolds go on the same folder - for scaffold
chr=Scaffold
mkdir -p /cluster/work/users/josece/var_call/${chr}
cd /cluster/work/users/josece/var_call/${chr}


echo "Doing ${region}"

bcftools mpileup --threads 4 -d 8000 --ignore-RG -r ${region} -a AD,DP,SP -Ou -f ${reference_genome} -b ${bams} | \
bcftools call --threads 4 -f GQ,GP -mO z -o ${region}.vcf.gz

touch ${region}_complete```

```
