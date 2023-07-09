### 01 - I got some fasta files with seeds for COI, 16s, and other mitochondrial genes.

### 02 - ran Novoplasty for the different seeds and saw whether I could get a mt genome.
```
NOVOPlasty4.2.pl -c config.txt
```

### 03 - concatenated them
```
# Processed the data first by making sure they all started in the same position.
mafft --auto --thread 30 mtgenomes.fasta > mtGenomes.fai
```

### 04 - removed ends and ran the tree!
```
iqtree -s all_genomes.EndsTrimmed.fai -nt 20 -B 1000
```
