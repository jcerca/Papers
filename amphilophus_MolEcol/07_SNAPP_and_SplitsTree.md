SNAPP is one of my favourite programs. It can be clunky though.

First, we need to reduce the loci - I'll get 4,000 random loci. Here's a very quick and dirty way of doing it:
```
grep -v "^#" amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.IndsWith0.5MissingnessRemoved.vcf | cut -f1-3 | shuf | head -n 4000 | cut -f 3 > random.snps
vcftools --vcf amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.IndsWith0.5MissingnessRemoved.vcf  --snps random.snps --recode --out trimmed.vcf
```

Then, I selected 10 individuals per population.  I kept the first 10 inds of each population after shuffling (made that part manually).

```
bcftools query -l trimmed.vcf | shuf > shuffled_inds.tsv
```

Now, we convert the vcf to nexus, using the script vcf2nex.pl available through the SNAPP website.
```
perl ../vcf2nex.pl trimmed.vcf > trimmed.nex
```

We also need to add 'missing=.' on the nex file, so missing data is processed. Do that with nano.

Load this file on SNAPP and:
1) format populations manually
2) uncheck "include non-polymorphic sites"
3) calculate mutation rates.

Set running:
```
module load beast2
module load beagle

snappfile=/usit/abel/u1/josece/003_amphilophus_radseq/2_20190319_newdata/11_SNAPP/2_VCF_with_only4000snps_indsRemoved/trimmed.indsRemoved.xml
filename=$(echo $snappfile | awk -F "/" '{print $10}')

cp $snappfile  $SCRATCH

## Do some work:
cd $SCRATCH
chkfile "."

beast -seed 223456 -threads $OMP_NUM_THREADS -beagle $filename
```


OK. Splitstree is easy. We convert the vcf
We first convert the vcf (no need just to have 4,000 variants, but don't do too many otherwise it may take ages).
I used this script: https://github.com/edgardomortiz/vcf2phylip

```
python vcf2phylip.py -i ../../1_cleaningTheVCF/4_NOTE_thereIsAnExtraPopulation_HA_which_I _removed/*vcf --nexus
```
and load the nexus file on the splitstree program. Done!
