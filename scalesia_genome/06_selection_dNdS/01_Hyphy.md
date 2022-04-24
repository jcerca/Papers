Signatures of selection on the whole genome.

We need CDS files (not aa sequences)

Folder structure:
00_genomeCDS
01_headers_cleaned
02_longestIsoforms
03_Orthofinder
04_OrthofinderParsingOfResults
05_singleCopy
06_multipleCopy_all
07_multipleCopy_notAll
08_topGOFormat
09_920selection_ARABI

00_genome_CDS
```
#here we will extract CDS from the genome (fa) and features file (gff3)
conda activate agat
# Extracting cds from subgenome A
perl ../../../local_bin/bitacora/Scripts/gff2fasta_v3.pl  ../../03_circosplots/02_Subgenomes_A_vs_B/00_data/scalesia_subgenomeA.fasta  ../../01_annotation/08_subgenomeGFF/annotation.subgenomeA.gff3  scalesia_subgenomeA
# Extracting cds from subgenome B
perl ../../../local_bin/bitacora/Scripts/gff2fasta_v3.pl  ../../03_circosplots/02_Subgenomes_A_vs_B/00_data/scalesia_subgenomeB.fasta  ../../01_annotation/08_subgenomeGFF/annotation.subgenomeB.gff3  scalesia_subgenomeB

ln -s ../../02_otherGenomes/01_Helianthus_annuus/*CDS* .
mv Ha412v1r1_CDS_v1.0.fasta HELIA.cds.fa

ln -s ../../02_otherGenomes/02_Arabidopsis/*cds* ./ARABI.cds.fa
perl ../../../local_bin/bitacora/Scripts/gff2fasta_v3.pl ../../02_otherGenomes/03_Conyza_canadensis/Conyza_canadensis_V1.final.fasta ../../02_otherGenomes/03_Conyza_canadensis/Conyza_canadensis_V1.final.fasta.gff CONYZ
ln -s ../../02_otherGenomes/04_Mikania_micrantha/*cds* ./MIKAN.cds.fa
ln -s ../../02_otherGenomes/05_Cynara_cardunculus/*cds* ./CYNAR.cds.fa
ln -s ../../02_otherGenomes/06_Lactuca_sativa/*cds*  ./LSATI.cds.fa
```

01_headers_cleaned
```
# Ran the kinfin script to remove very short sequences and sanitize headers. Did this for every genome.
# e.g.
# filter_fastas_before_clustering.py -f ARABI.cds.fa > ../02_longestIsoforms/ARABI.cds_cleanedHeaders.fa
```

02_longestIsoforms

```
# Arabidopsis and Conyz have isoforms, here we remove them
# ARABI
conda activate compgen
cp ../01_headers_cleaned/ARABI.cds_cleanedHeaders.fa . ; samtools faidx ARABI*fa
awk -F'[\t.]' '{print $1,$2,$3,$4}' ARABI*.fai | sort -k4nr,4 | sort -uk1,2 | cut -f1-3 -d' '| tr ' ' '.' > selection.ARABI

while read contig; do
counter=`expr $counter + 1`; echo $counter;
samtools faidx ARABI*fa $contig >> ARABI.cleanedHeaders_longestIsoforms.fa;
done < selection.$i.ls; cd ..;
done

#De-interleaving:
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ARABI.cleanedHeaders_longestIsoforms.fa > tmp
mv tmp ARABI.cleanedHeaders_longestIsoforms.fa

# CONYZ
conda activate compgen
cp ../01_headers_cleaned/CONYZ.cds_cleanedHeaders.fa . ; samtools faidx CONYZ*fa
awk -F'[\t.]' '{print $1,$2,$3,$4}' CONYZ*.fai | sort -k4nr,4 | sort -uk1,2 | cut -f1-3 -d' '| tr ' ' '.' > selection.CONYZ

while read contig; do
counter=`expr $counter + 1`; echo $counter;
samtools faidx CONYZ.cds_cleanedHeaders.fa $contig >> CONYZ.cleanedHeaders_longestIsoforms.fa;
done < selection.CONYZ; cd ..;
done

#De-interleaving
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CONYZ.cleanedHeaders_longestIsoforms.fa > tmp
mv tmp CONYZ.cleanedHeaders_longestIsoforms.fa
```

03_Orthofinder
```
# Orthofinder with the -d (for DNA) - copied from 013*/05*/09*/cds/
ls ../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/
orthofinder -d -f . -t 10
```


04_Orthofinder
```
#parsing of results
# An e-mail from JVizueta recommened me establishing three groups:
#"Regarding your question about selection, that sounds really strange to select one gene randomly for the analyses. What I have typically read and done in my project, is to classify the different orthogroups (obtained from Orthofinder) into different categories: single copy (1 copy in all species), multiple copy or gene families (more than 1 copy present in all species), or other groups which could be missing in certain species (for instance in 3 species, it could be 1,1 and no copy..). Then, you can conduct for every group selection analyses, but the most informative and used would be the single copy. Then, the fact of selecting orthogroups with more than 1 copy and discarding extra copies sounds really weird and I don’t think this is correct (just saying from my point of view)."

# Getting 01_singleCopy.tsv
# To get single orthologues is easy.
ls ../../05_Orthofinder/06_Orthofinder_all_ScalesiaSubgenomes/OrthoFinder/Results_Feb08/Single_Copy_Orthologue_Sequences/ > 01_singleCopy.tsv

# Getting a list for single copy: 01_singleCopy.tsv
ls ../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Single_Copy_Orthologue_Sequences/ > 01_singleCopy.tsv


# Getting multipleCopy present on all species
# To get orthologues present in all genomes, but removing only those single copy. Adding ".fa", so we can compare with 01
grep "CONYZ" ../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroups/Orthogroups.txt | grep "LSATI" | grep "HELIA" | grep "SCSGA" | grep "SCSGB" | grep "MIKAN" | awk -F ":" '{print $1".fa"}' > tmp


#checking if it worked..
comm -12 01_singleCopy.tsv tmp | wc -l
1354
wc -l 01_singleCopy.tsv
1354 01_singleCopy.tsv

# Now supressing the single-copy ones:
comm -13 01_singleCopy.tsv tmp > 02_multipleCopy.allsp.tsv

# and checking whether it worked.
wc -l tmp 02_multipleCopy.allsp.tsv
  7273 tmp
  5919 02_multipleCopy.allsp.tsv

# It did!

# Getting  All the others (not single copy, not multiple copy in all sp)
###
ls ../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/ > tmp
ls ../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/ | wc -l
53,846 # 53,846, to which we have to remove the ones which are single copy (01) and those that are multiple copy, and in all sp (02).

man comm
comm -23 tmp 01_singleCopy.tsv > tmptmp
wc -l tmptmp
52,492 tmptmp
comm -23 tmptmp 02_multipleCopy.allsp.tsv > tmptmptmp
wc -l tmptmptmp
46,573 tmptmptmp # Perfect! It is 53,846 (all) = 46,573 + 5,919 (02) + 1,354 (01)

mv tmptmptmp 03_multipleCopy.notInAllgenomes.tsv


# First, we get the genes belonging to orthogroups, so we can grep them out.
while read orthoG; do tmp=$(echo $orthoG | sed "s/.fa//"); grep "$tmp" ../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroups/Orthogroups.tsv >> tmp_OrthogroupsWithGenes ; done  < 03_multipleCopy.notInAllgenomes.tsv

# Decision, keep only groups:
#	1. the subgenomes are present
grep "SCSG" tmp_OrthogroupsWithGenes > tmp_OrthogroupsWithGenes_firstclean

#	2. either of the subgenomes AND another sp
grep "CONYZ\|HELIA\|MIKAN\|LSATI" tmp_OrthogroupsWithGenes_firstclean > tmp_OrthogroupsWithGenes_secondclean

#	3. removing orthogroups with 4 or less genes
head -n 5629 tmp_OrthogroupsWithGenes_secondclean > tmp_OrthogroupsWithGenes_thirdclean

tmp_OrthogroupsWithGenes_thirdclean
awk '{print $1".fa"}' tmp_OrthogroupsWithGenes_thirdclean > 03_multipleCopy.notInAllgenomes.selected.tsv

```

Now, the analysis is the same for:
05_singleCopy
06_multipleCopy_all
07_multipleCopy_notAll

I'll only show the ones for 07_multipleCopy_notAll
Inside 07_multipleCopy_notAll we have the following folders:
	01_orthogroups
	02_Zorro
	03_ZorroResults
	04_bestOGfromZorro_namesCleaned
	05_Trees
	06_stopCodonsRemoved
	07_HyPhy
	08_parsingOutData
	99_manualinstallHyPhy

07_multipleCopy_notAll/01_orthogroups
```
# This takes a while, so splitting them like this.
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[0].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[1].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[2].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[3].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[4].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[5].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[6].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[7].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[8].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
while read GENE; do echo $GENE; if [[ ${GENE} =~ ^OG.*[9].fa$ ]] ; then prank -d=../../../05_Orthofinder/09_Orthofinder_subgenomes_HELIA_MIKAN_LSATI_CONYZ/cds/OrthoFinder/Results_Feb19/Orthogroup_Sequences/${GENE} -o=${GENE%.fa}.aligned.fas -codon; fi; done < ../../04_OrthofinderParsingOfResults/03_multipleCopy.notInAllgenomes.selected.tsv
```

07_multipleCopy_notAll/02_Zorro
```
conda activate phyl
for i in *0.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
for i in *1.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
for i in *2.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
for i in *3.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
for i in *4.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
for i in *5.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
for i in *6.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
for i in *7.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
for i in *8.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
for i in *9.aligned.fas.best.fas; do echo $i; ../../../../local_bin/zorro/zorro.sh $i > ../02_Zorro/${i%.fas.best.fas}.zorro; done
```

07_multipleCopy_notAll/03_ZorroResults
```
for i in *zorro; do avg=$(awk '{ sum += $1 } END { if (NR > 0) print sum / NR }' $i); echo "${i} ${avg}"; done > ../03_ZorroResults/zorro.avg
while read gene; do sample=$(echo ${gene} | sed "s/.aligned.*//"); echo $sample; cp ../01_orthogroups/${sample}* .; sed -i "s/CONYZ.*/CONYZ/; s/SCSGA.*/SCSGA/; s/HELIA.*/HELIA/; s/SCSGB.*/SCSGB/; s/MIKAN.*/MIKAN/; s/LSATI.*/LSATI/" ${sample}*; done <  ../03_ZorroResults/zorro.bestAlign.avg

# 07_multipleCopy_notAll/04_bestOGfromZorro_namesCleaned
# We clean the headers and just copy these
cd 04_bestOGfromZorro_namesCleaned

#I used a VARIANT of this code (forgot to save it on time..)
while read gene; do sample=$(echo ${gene} | sed "s/.aligned.*//"); echo $sample; cp ../01_orthogroups/${sample}.aligned.fas.best.fas .; sed -i '/^>/s/-/_/g; s/\./_/g' ${sample}.aligned.fas.best.fas; mv ${sample}.aligned.fas.best.fas ${sample}.aln.best.cleaned.fa; done < ../03_ZorroResults/zorro.bestAlign.avg

#### IN THIS ONE WE WILL RUN THE TREES ONLY AFTER CLEANING FOR STOP CODONS BECAUSE HYPHY ALSO REMOVED SOME OF THE TIPS.
````

07_multipleCopy_notAll/05_stopCodonsRemoved
```
### We need to clean stop codons.
# To do so, I made an image of HyPhy, by doing:
cd ../99_*
git clone https://github.com/veg/hyphy.git
cp -r hyphy/res/TemplateBatchFiles/TemplateModels .

# This basically allows us to run hyphy like this:
hyphy CleanStopCodons.bf -h # to remove codons

# So the first we do is to go back to 04_* and remove codons. But the problem is that files will be created only for those that stop codons are removed.
cd ../04*
conda activate compgen
for FILE in `ls *.fa`; do echo $FILE; cp $FILE ../99_manualinstallHyPhy/; cd ../99_manualinstallHyPhy/; echo "Running CleanStopCodons on $FILE"; echo `(echo "1"; echo $FILE; echo "5"; echo ${FILE%.fas}.cleaned) | HYPHYMPI CleanStopCodons.bf`; mv *cleaned ../05_stopCodonsRemoved; rm $FILE; cd ../04_bestOGfromZorro_namesCleaned; done

# code explained
# for FILE in `ls *.fas`; do #### Open a loop on every fasta
# echo $FILE; cp $FILE ../99_manualinstallHyPhy/; cd ../99_manualinstallHyPhy/; ### copy the file to the folder where we can run the CleanStopCodons, and move there
# echo "Running CleanStopCodons on $FILE"; #### Just an Echo
# echo `(echo "1"; echo $FILE; echo "5"; echo ${FILE%.fas}.cleaned) | HYPHYMPI CleanStopCodons.bf`;  ## Actual hyphy command for this case
# mv *cleaned ../05_stopCodonsRemoved; rm $FILE; cd ../04_bestOGfromZorro_namesCleaned; done #### move to 05 and remove the copy we created

# Now we go back to 05_stop codon removed, and copy all the original files there.
cd ../05*
cp ../04_bestOGfromZorro_namesCleaned/*fa .


# We ran a "if" statement to see how many files to remove IF there is its counterpart
# checking how many we would remove.. - 598
for i in *fa; do if [ -f "${i}".cleaned ] ; then     echo "rm $i"; fi; done | wc -l
2125

# checking how many "cleaned" were generated.. - 598
ls *cleaned | wc -l
2125


# They match perfectly, so we can move on
for i in *fa; do
        if [ -f "${i}".cleaned ] ; then
    rm $i
fi;
done

#Let's also homogeneize the names
for i in *; do echo ${i%.aln.*}; mv $i ${i%.aln.*}.InputForSelection.fas;  done


```
07_multipleCopy_notAll/06_Tree
```
cd ../05*

for i in *fas; do
        cp $i ../06_Trees;
        cd ../06_Trees;
        echo "Tree for $i"
        iqtree -nt 1 -s $i -bb 1000;
        cd ../05*;
done
```

07_multipleCopy_notAll/08_parsingOutData
```
cd ../07*
grep "^\* SCSG" *selection > ../08_parsingOutData/ScalesiaGenesUnderSelection.tsv

cd ../08*
awk -F "." '{print $1}' ScalesiaGenesUnderSelection.tsv  | sort | uniq -c  | wc -l 386 OGs under selection
349 orthogroups with evidence for selection

grep SCSGB ScalesiaGenesUnderSelection.tsv | sed "s/.*SCSGB/SCSGB/; s/,.*//" > underSelectionB.tsv
grep SCSGA ScalesiaGenesUnderSelection.tsv | sed "s/.*SCSGA/SCSGA/; s/,.*//" > underSelectionA.tsv
```

Now, we do the topGo.
08_topGOFormat
Folder structure within 08_topGOFormat

08_topGOFormat
```
cat ../05_singleCopy/08_parsingOutData/underSelection.sub[AB].tsv ../06_multipleCopy_all/08_parsingOutData/underSelection* ../07_multipleCopy_notAll/08_parsingOutData/underSelection* > allgenes_selection.tsv

wc -l allgenes_selection.tsv
920

grep SCSGA allgenes_selection.tsv  | sed "s/SCSGA[_.]//; s/\-/_/g" > subA.selection.tsv
grep SCSGB allgenes_selection.tsv  | sed "s/SCSGB[_.]//; s/\-/_/g" > subB.selection.tsv


cp ../../15_swissProt/07_TopGO_format/scalesia.tsv .

sed "s/SCALE.//; s/_HRSCAF_[0-9]\+//; s/\-/_/g" scalesia.tsv > scalesia.TopGo.tsv ; rm scalesia.tsv
sed -i "s/\./_/" scalesia.TopGo.tsv


sed -i "s/\./_/" subA.selection.tsv
sed -i "s/\./_/" subB.selection.tsv


while read gene; do grep $gene scalesia.TopGo.tsv; done < subA.selection.tsv > subA.selection.TopGo.tsv
while read gene; do grep $gene scalesia.TopGo.tsv; done < subB.selection.tsv > subB.selection.TopGo.tsv


wc -l sub*
   391 subA.selection.TopGo.tsv
   478 subA.selection.tsv
   356 subB.selection.TopGo.tsv
   442 subB.selection.tsv

Why different numbers? Simply becasue some genes have no GOs.

cat subA.selection.TopGo.tsv subB.selection.TopGo.tsv > scalesia.selection.TopGo.tsv

mv *TopGo.tsv ../02*
for i in *selection*; do echo $i; awk '{print $1}' $i > ${i%.tsv}.genelist.tsv; done
```

Now, we get the Scalesia vs Arabidopsis orthologs.
09_920selection_ARABI/

```
# getting the genes under selection, and cleaning names:
cat ../08_topGOFormat/01_processingData/allgenes_selection.tsv | sed "s/-/_/g; s/\./_/g" > 920genes_selection_namesHomogeneized.tsv

# Getting orthoglogues, and cleaning names:
cat ../../05_Orthofinder/07_Orthofinder_SugenomesARABI/OrthoFinder/Results_Feb11/Orthogroups/Orthogroups.tsv | sed "s/-/_/g; s/\./_/g" | sed "s/=[0-9]\+_/_/g; s/_HRSCAF//g" > arabiOrthogroups.tsv

grep -F -w -f 920genes_selection_namesHomogeneized.tsv arabiOrthogroups.tsv

#Making the orthogroup list
grep -F -w -f 920genes_selection_namesHomogeneized.tsv arabiOrthogroups.tsv | grep ARABI > OGs_arabidopsis_scalesiaCorrespondent.tsv

# Making a final list of Arabidopsis genes
grep -F -w -f 920genes_selection_namesHomogeneized.tsv arabiOrthogroups.tsv | grep ARABI | tr "\t" "\n" | tr "," "\n" | grep ARABI | sed "s/ARABI_//" | sed  "s/_/\./" | sed "s/ //" > arabidopsis_scalesiaCorrespondent.tsv

# Retrieving information
grep -w -F -f arabidopsis_scalesiaCorrespondent.tsv ../../02_otherGenomes/02_Arabidopsis/Araport11_genes.201606.pep.fasta  > gene_ids.tsv
```
