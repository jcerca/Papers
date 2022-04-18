Folder structure:
00_genomes, 01_sub_genomes, 02_headersSanitized, 03_selectingTheLongestIsoforms, 04_mergedDatasets_Cafe_Kinfin, 05_Orthofinder_all_ScalesiaWholeGenome, 06_Orthofinder_all_ScalesiaSubgenomes, 07_Orthofinder_SugenomesARABI, 08_Orthofinder_subgenomes_HELIA_MIKAN_LSATI, 09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ, 10_Orthofinder_subgenomes_HELIA, 11_Orthofinder_wholegenome_HELIA_MIKAN_LSATI_CONYZ

00_genome
```
# Cleaning up headers and files by renaming them.
cd 00_genome
mv arabidopsis_thaliana.faa ARABI.arabidopsis_thaliana.faa
mv cynara_cardunculus.faa CYNAR.cynara_cardunculus.faa
mv helianthus_annuus.faa HELIA.helianthus_annuus.faa
mv lactuca_salina.faa LSALI.lactuca_salina.faa
mv lactuca_sativa.faa LSATI.lactuca_sativa.faa
mv mikania_micrantha.faa MIKAN.mikania_micrantha.faa
mv scalesia_atractyloides.faa SCALE.scalesia_atractyloides.faa

# We clean the headers.
for i in *faa; do /data/bigexpansion/jcerca/local_bin/kinfin/scripts/filter_fastas_before_clustering.py  -f $i > tmp; mv tmp ../02_headersSanitized/${i%.faa}_cleanedHeaders.faa; done
```
01_sub_genomes
```
# Same with the subgenomes
cd 01_sub_genomes
mv scalesia_subgenomeA.faa SCSGA.scalesia_subgenomeA.faa
mv scalesia_subgenomeB.faa SCSGB.scalesia_subgenomeB.faa

for i in *faa; do /data/bigexpansion/jcerca/local_bin/kinfin/scripts/filter_fastas_before_clustering.py  -f $i > tmp; mv tmp ../02_headersSanitized/${i%.faa}_cleanedHeaders.faa; done

```

We need/want the longest isoforms (NCBI files should already be formatted in such a way), but not all files.
CREDIT https://bioinformatics.stackexchange.com/questions/595/how-can-longest-isoforms-per-gene-be-extracted-from-a-fasta-file

03_selectingTheLongestIsoforms
```
cd 03_*
# We make a directory
for i in ARABI CYNAR HELIA LSALI LSATI MIKAN SCALE SCSGA SCSGB; do mkdir $i; done

# Index the assemblies
for i in ARABI CYNAR HELIA LSALI LSATI MIKAN SCALE SCSGA SCSGB; do cd $i; cp ../../02_headersSanitized/${i}*faa . ; samtools faidx ${i}*faa; cd ..;  done

# Then, we get the biggest isoforms
for i in ARABI CYNAR HELIA LSALI LSATI MIKAN SCALE SCSGA SCSGB; do cd $i; awk -F'[\t.]' '{print $1,$2,$3,$4}' *fai | sort -k4nr,4 | sort -uk1,2 | cut -f1-3 -d' '| tr ' ' '.' > selection.$i.ls; cd ..; done

# Now, we extract them
for i in ARABI CYNAR HELIA LSALI LSATI MIKAN SCALE SCSGA SCSGB  ; do
cd $i; counter=0; wc -l selection.$i.ls;
while read contig; do
counter=`expr $counter + 1`; echo $counter;
samtools faidx ${i}*_cleanedHeaders.faa $contig >> ${i}.sanitizedHeaders.longestisoforms.faa;
done < selection.$i.ls; cd ..;
done

```

Orthofinder run example

```
# What I love about OrthoFinder, is that it is easy.
orthofinder -f . -t 10
```
