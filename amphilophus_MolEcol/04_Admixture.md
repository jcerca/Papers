Stop. [Hammer] Admixture time.


First, we will convert the data using plink v1.9. Remember we wantthe LD-filtered dataset.

```
VCF=amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.LDprunned.vcf.IndsWith0.5MissingnessRemoved.vcf
plink --vcf $VCF --recode12 --out admixture.input --allow-extra-chr --double-id --make-bed
```

Now; admixture wants chromossomes to be integers and not CCOE01 (like we have on this work) ... so we do some black magic, which will convert the chromosome Ids to numbers. This is dataset specific.

```
cat admixture.input.bim | sed -s "s/ENA|CCOE........|CCOE//" | sed -s "s/\.1//" | less > newinput.bim


head -n 5 newinput.bim # see the first column
01000001        ENA|CCOE01000001|CCOE01000001.1:83630   0       83630   A       C
01000001        ENA|CCOE01000001|CCOE01000001.1:128953  0       128953  C       G
01000001        ENA|CCOE01000001|CCOE01000001.1:130030  0       130030  T       G
01000001        ENA|CCOE01000001|CCOE01000001.1:154460  0       154460  C       T
01000001        ENA|CCOE01000001|CCOE01000001.1:154578  0       154578  G       T
mv newinput.bim admixture.input.bim 
```

Now, we create a new file for plotting:

```
awk -F "_" '{print $1"\t" $3 "\t" $4}' admixture.input.fam > ind.species.pop.tsv

head -n 6 ind.species.pop.tsv
13NIei01        C       Om
13NIei16        C       Om
13NIei17        C       Om
13NIei18        C       Om
13NIei19        C       Om
13NIei20        C       Om
```

Now the actual ADMIXTURE run:
````
INPUT=/usit/abel/u1/josece/003_amphilophus_radseq/2_20190319_newdata/3_admixture/1_fileconvertion/admixture.input.bed
MIX=/usit/abel/u1/josece/000_locally_installed_stuff/miniconda2/bin/admixture

for i in 2 3 4 5 6; do K=${i}; echo "This is my K for today = $i";
cd K$i ;
mkdir run1 run2 run3 run4 run5 ; 
cd run1 ;
$MIX --cv $INPUT $K -j4 | tee log ;
cd .. ; 

cd run2 ;
$MIX --cv $INPUT $K -j4 | tee log ;
cd .. ;

cd run3; 
$MIX --cv $INPUT $K -j4 | tee log; 
cd .. ;

cd run4 ;
$MIX --cv $INPUT $K -j4 | tee log ;
cd .. ;

cd run5 ;
$MIX --cv $INPUT $K -j4 | tee log ;
cd .. ;

cd .. ;
done
````
