We can either use the VCF or the phylogenetic matrix.

matrix-based:
```
# We convert 
# https://github.com/tkchafin/scripts/blob/master/fasta2nexus.pl
fasta2nexus.pl FcC_supermatrix.fas > FcC_supermatrix.nex
```

vcf-based:

```
VCF=Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.indswith90missingnessRemoved.randomSNP.vcf
# https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py

python vcf2phylip.py -i $VCF --nexus
```

Now load it on splitstree, which is a point and click program.
