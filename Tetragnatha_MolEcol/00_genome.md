# I am using the version with non-ambiguousbasepairs!
tetragnatha_kauaiensis_11Sep2016_pjasG.NonAmbiguousBP.fasta.gz

mv tetragnatha_kauaiensis_11Sep2016_pjasG.NonAmbiguousBP.fasta.gz T_kauaiensis_ref.fasta

# I indexed the genome
ml BWA/0.7.17-intel-2018b
bwa index T_kauaiensis_ref.fasta

# I created a reference dictionary
ml picard/2.18.27-Java-1.8
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=T_kauaiensis_ref.fasta O=T_kauaiensis_ref.dict

# create the .fai file
samtools faidx T_kauaiensis_ref.fasta
