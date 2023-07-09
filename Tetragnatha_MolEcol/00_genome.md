### 00 - I am using the version with non-ambiguousbasepairs!
```
mv tetragnatha_kauaiensis_11Sep2016_pjasG.NonAmbiguousBP.fasta.gz T_kauaiensis_ref.fasta
```

### 01 - I indexed the genome
```
bwa index T_kauaiensis_ref.fasta
```

### 02 - I created a reference dictionary
```
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=T_kauaiensis_ref.fasta O=T_kauaiensis_ref.dict
```

### 03 - create the .fai file
```
samtools faidx T_kauaiensis_ref.fasta
```
