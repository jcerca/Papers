### Isoseq Transcriptome assembly
#### Folder structure: 00_raw, 01_limaCleaning, 02_isoSeqRefine, 03_isoSeqClustering, 04_alignToRefGenome, 05_SQANTI3_qc, 06_SQANTI3_filter

00_raw
# primers.fasta contains the pacbio primers. I obtained them from https://www.pacb.com/wp-content/uploads/SMRT-Tools-Reference-Guide-v8.0.pdf
```
# We make a symbolic link to the genome
ln -s ../../../00_ScalesiaGenome/00_fromDovetail/scalesia_atractyloides.fasta ./scalesia_atractyloides.fasta

# Load the tools
conda activate isoseq

# Get indexes for the genome
samtools faidx scalesia_atractyloides.fasta
```
cd 01_limaCleaning/
samtools faidx scalesia_atractyloides.fasta```

```
# Now we'll be cleaning the raw (transcriptome) reads with lima, by trimming primers
lima --isoseq --dump-clips --peek-guess --num-threads 10 ../00_raw/R28_r64067_200927_001.ccs.bam ../00_raw/primers.fasta s_atractyloides.IsoSeq.cleanedByLima.bam
less *summary
# The first indicators indicate how many reads passed your filters.
```

02_isoSeqRefine/ 
```
# Now we will remove polyAAAAAA tails and artifictial concatemers
isoseq3 refine --require-polya ../01_limaCleaning/s_atractyloides.IsoSeq.cleanedByLima.primer_5p--primer_3p.subreadset.xml \
../00_raw/primers.fasta scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.bam

less *json
# Tells you num_reads_fl (input), and then "num_reads_flnc_ploya" (output; reads that have polyAtails and not concatemeters)
```

cd 03_isoSeqClustering
```
# Clustering takes the bam and clusters by isoform level. Not parallelizable since it reads all reads at once; in other words, memory hungry!

isoseq3 cluster ../02_isoSeqRefine/scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.bam \
scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.polished --verbose --use-qvs --num-threads 10

# After clustering we are identifying transcripts which are similar, not exploring which ones belong to the same gene (we havent added the genome, so we dont know where they fall)
head *csv
# each row is a different read. E-g-:
# transcript/0,m64014_190506_005857/69535114/ccs,FL  < reads belonging to the first cluster, named "transcript/0"
# transcript/0,m64014_190506_005857/93062567/ccs,FL  < reads belonging to the first cluster, named "transcript/0"
# transcript/1,m64014_190506_005857/166200685/ccs,FL < reads belonging to the second cluster, named "transcript/1"
```

04_alignToRefGenome/
```
pbmm2 align /data/bigexpansion/jcerca/013_ScalesiaGenome/00_ScalesiaGenome/00_fromDovetail/scalesia_atractyloides.fasta \
../03_isoSeqClustering/scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.polished.transcriptset.xml \
scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.polished.alignedToRef.bam --preset ISOSEQ --sort -j 24 --log-level INFO

isoseq3 collapse scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.polished.alignedToRef.bam \
scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.polished.alignedToRef.colapsed.gff

# This error on the second comand is totally fine:
#>|> 20201111 15:22:35.113 -|- WARN -|- Run -|- 0x7fab48c4d7c0|| -|- Transcripts do not contain quality values, will not output scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.polished.alignedToRef.colapsed.fastq!
```
05_SQANTI3_qc 

```
source activate SQANTI3.env
ln -s ../../../../local_bin/SQANTI3/sqanti3_qc.py .
export PYTHONPATH=$PYTHONPATH:/data/bigexpansion/jcerca/local_bin/SQANTI3/cDNA_Cupcake/
export PYTHONPATH=$PYTHONPATH:/data/bigexpansion/jcerca/local_bin/SQANTI3/cDNA_Cupcake/sequence/

# The programs needs symbolic links on the folder, otherwise it'll mess up.
ln -s ../04_alignToRefGenome/scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.polished.alignedToRef.colapsed.abundance.txt ./abundance.txt
ln -s ../00_raw/scalesia_atractyloides.fasta ./ref_genome.fa
ln -s ../04_alignToRefGenome/scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.polished.alignedToRef.colapsed.fasta ./isoforms.fasta
ln -s ../04_alignToRefGenome/scalesia_atractyloides.cleanedByLima.polyAAAtailsRemoved.polished.alignedToRef.colapsed.gff ./annotations.gff

python sqanti3_qc.py -o S_atractyloides -n4 -t 8 --fl_count --polyA_motif polyA.txt \
--fl_count abundance.txt \
./isoforms.fasta annotations.gff ref_genome.fa

# This gives a lot of warnings.. they're ignoring them and advising one to ignore them.
# isoforms.fasta                                 - long read defined transcripts
# referenceTranscriptome.gtf                     - reference transcriptome
# scalesia_atractyloides.fasta                   - reference genome
# -o S_atractyloides                             - output prefix and directory
# --fl_count abundance.txt                       - abundance.txt file obtained from AnalysisPT1!!!!
# --polyA_motif polyA.txt                        - based on the literature

# See the pdf
# see the classification.txt file
```

06_SQANTI3_filter 
```
ln -s ../../../../local_bin/SQANTI3/sqanti3_RulesFilter.py
/data/bigexpansion/jcerca/013_ScalesiaGenome/01_annotation/03_Isoseq_transcriptome/06_SQANTI3_filter$ python sqanti3_RulesFilter.py ../05_SQANTI3_qc/S_atractyloides_classification.txt ../05_SQANTI3_qc/S_atractyloides_corrected.fasta ../05_SQANTI3_qc/S_atractyloides_corrected.gtf -c 3 -a 0.6
# Tells you how many isoforms were removed
```

Done!
