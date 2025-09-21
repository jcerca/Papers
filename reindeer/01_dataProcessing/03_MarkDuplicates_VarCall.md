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

And, finally, we normalize them:

```
# Normalizing files
bcftools view -V indels -e 'ALT="*" | N_ALT>1' ${chr}.vcf.gz | bcftools norm -D -O z -o ${chr}_norm.vcf.gz
```
