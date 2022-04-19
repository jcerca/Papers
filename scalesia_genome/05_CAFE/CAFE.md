Folder structure:
00_data_headers_sanitized_longestIsoforms, 01_allByallBLAST, 02_clustering, 03_Tree, 04_cafe_filteredData, 05_genes

Let's go!
00_data_headers_sanitized_longestIsoforms
```
conda activate cafe
cd 00_data_headers_sanitized_longestIsoforms/
# First, we concatenate all the aminoacid files on a file called:
CONYZ_HELIA_LSATI_MIKAN_SCALE.faa
```

01_allByallBLAST
```
# We make a blast database
makeblastdb -in ../00_data_headers_sanitized_longestIsoforms/CONYZ_HELIA_LSATI_MIKAN_SCALE.faa -dbtype prot -out blastdb

# And do a all-by-all Blast
blastp -num_threads 30 \
-db blastdb \
-query ../00_data_headers_sanitized_longestIsoforms/CONYZ_HELIA_LSATI_MIKAN_SCALE.faa \
-outfmt 7 \
-seg yes > blast_output.txt

```
02_clustering
```
# Now, we cluster
cd 02_clustering

# First, we get an .abc file.
grep -v "#" ../01_allByallBLAST/blast_output.txt | cut -f 1,2,11 > blast_output.abc

#Then, create a network..
mcxload -abc blast_output.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o blast_output.mci -write-tab blast_output.tab

#Then create a dictionary..
mcl blast_output.mci -I 3
mcxdump -icl out.blast_output.mci.I30 -tabr blast_output.tab -o dump.blast_output.mci.I30

#parse the files
python ../../00_pythonScripts/cafetutorial_mcl2rawcafe.py -i dump.blast_output.mci.I30 -o unfiltered_cafe_input.txt -sp "CONYZ HELIA LSATI MIKAN SCALE"

# remove gene families with large variance (>100 gene copies) from the dataset.
python ../../00_pythonScripts/cafetutorial_clade_and_size_filter.py -i unfiltered_cafe_input.txt -o filtered_cafe_input.txt -s

	# BASICALLY:
	# filtered_cafe_input.txt # was filtered for gene-families with astonishing variation
	# unfiltered_cafe_input.txt # contains_all
	# large_filtered_cafe_input.txt - families with >100 variation of genes
	# filtered_families - families with less than <100 variation
```
03_Tree
```
#OK - goal is the make the tree, we have a tree from the OrthoFinder run

ln -sf ../../../05_Orthofinder/11_Orthofinder_wholegenome_HELIA_MIKAN_LSATI_CONYZ/OrthoFinder/Results_Feb26/Species_Tree/SpeciesTree_rooted.txt .

# Now we use the Dendroscope graphical user interface.
# Load up your tree on dendroscope (copy-pasting the newik format works)
# Select the node which is supposed to be rooted. Click edit>reroot
# On each bootstrap support (displayed on the nodes), click on it and "edit node label", and delete bootstrap support
# Export it as newik (I named it as tree.bootstrapremoved.treefile)

# Here is how mine looks:
(CONYZ:0.0877205,(((SCALE:0.0703455,HELIA:0.11523):0.0730229,MIKAN:0.137656):0.0461635,LSATI:0.149095):0.0877205);

# Calculated on time tree LSATI and HELIA - 39 million year divergence.
python ../../00_pythonScripts/cafetutorial_rep_r8s.py -i tree.tree -o r8s_ctl_file.txt -s 430148 -c '23' -p 'HELIA,MIKAN'

        #-o for output
        #-s for the number of sites, I estimated using the number of sites in Scalesia: 145471
        #-p is the species pairs we are constraining for the time of divergence - I selected CYRTO and PRIMA and "11.1", as the number of millions of years, as reported in: http://timetree.org/

        # Now, to run r8s ... (ps - installing it was painful, and after trying on 4 different clusters, I found a docker image that worked.

# It also worked doing:
#       1. docker pull shkao/r8s:1.81
#       2. opening docker
#       3. There is no SCP so we have to do echo ""

r8s -b -f r8s_ctl_file.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > r8s_ultrametric.txt

#scp back!


#BEFORE WE PROCEEED! r8s adds a NODE-name based on the node that you "calibrated" the tree at. It added "LIAKAN" on my run. I created a new file:
cp r8s_ultrametric.txt r8s_ultrametric.modified.txt
nano r8s_ultrametric.modified.txt # removed "LIAKAN"
```

04_cafe_filteredData
```
mkdir reports

# run cafe
cafe

#inside cafe type:
#NOTE! THE TREE CAN'T HAVE BOOTSTRAP SUPPORT. Load up the tree in dendroscope to see the error.
load -i ../01_allByallBLAST/filtered_cafe_input.txt -t 50 -l reports/log_run1.txt

# Now we need the tree:
# cat ../../04_Orthologs_and_tree/05_making_theTree_Ultrametric/r8s_ultrametric.modified.txt

tree ((((SCALE:13.853014,HELIA:13.853014):9.146986,MIKAN:23.000000):5.395652,LSATI:28.395652):10.253031,CONYZ:38.648683);


# We need to specify that the LAMBDA PARAMETER will be similar for the whole tree..
# We do this by replacing every node and tip by a 1. It has to have the same topology as the tree above.
lambda -s -t ((((1,1)1,1)1,1)1,1);

#result: lambda         0.02582077447885

#And get the file files..
report reports/report_run1

# Now

# This code is now different!
python  ../../00_pythonScripts/cafetutorial_report_analysis.py -i reports/report_run1.cafe -o summary_run1
#It is easy to get the results from the ./05_CAFE/reports folder (see tutorial pages 10-11

#(This does no longer work..)
# python ../../00_pythonScripts/cafetutorial_draw_tree.py -i reports/report_run1.cafe -t '(((((SCALE:13.5703,HELIA:13.5703):9.42972,MIKAN:23):5.72475,(CYNAR:23.7437,LSATI:23.7437):4.98103):4.27841,CON
YZ:33.0032):17.8782,ARABI:50.8813)' -d '(((((SCALE<0>,HELIA<2>)<1>,MIKAN<4>)<3>,(CYNAR<6>,LSATI<8>)<7>)<5>,CONYZ<10>)<9>,ARABI<12>)<11>' -o summary_run1_Tree_rapid.png -y Rapid
# python ../../00_pythonScripts/cafetutorial_draw_tree.py -i reports/report_run1.cafe -t '(((((SCALE:13.5703,HELIA:13.5703):9.42972,MIKAN:23):5.72475,(CYNAR:23.7437,LSATI:23.7437):4.98103):4.27841,CON
YZ:33.0032):17.8782,ARABI:50.8813)' -d '(((((SCALE<0>,HELIA<2>)<1>,MIKAN<4>)<3>,(CYNAR<6>,LSATI<8>)<7>)<5>,CONYZ<10>)<9>,ARABI<12>)<11>' -o summary_run1_Tree_contractions.png -y Contractions
# python ../../00_pythonScripts/cafetutorial_draw_tree.py -i reports/report_run1.cafe -t '(((((SCALE:13.5703,HELIA:13.5703):9.42972,MIKAN:23):5.72475,(CYNAR:23.7437,LSATI:23.7437):4.98103):4.27841,CON
YZ:33.0032):17.8782,ARABI:50.8813)' -d '(((((SCALE<0>,HELIA<2>)<1>,MIKAN<4>)<3>,(CYNAR<6>,LSATI<8>)<7>)<5>,CONYZ<10>)<9>,ARABI<12>)<11>' -o summary_run1_Tree_expansions.png -y Expansions

```

OK. Let's process the data. Inside ./05_genes we have several folders:
01_expansions_contractions, 02_family, 03_HeliaGenes, 04_GOs, 05_final_files, 06_FilesForTopGo, 07_Revigo


01_expansions_contractions

```
cd 01_expansions_contractions

cat ../../../04_cafe_filteredData/summary_run1_fams.txt | grep "^SCALE" | tr "]" "\n" | tr "\t" "\n" | grep "+" | sed "s/,//; s/\[+/\t/; s/*//" | sed '1 i\geneFamily\tExpansion' > scalesia.expansions
cat ../../../04_cafe_filteredData/summary_run1_fams.txt | grep "^SCALE" | tr "]" "\n" | tr "\t" "\n" | grep "-" | sed "s/,//; s/\[-/\t/; s/*//" | sed '1 i\geneFamily\tContraction' > scalesia.contractions
```

02_family
```
cd 02_family
for i in 184 213 226 266 283 310 317 371 390 403 410 416 448 463 482 490 527 550 609 621 635 674 695 720 722 730 765 785 852 1014 1022 1209 1308 1317 1406 1428 1716; do echo "This is gene family number $i"; head -n $i ../../../01*/dump.blast_output.mci.I30 | tail -1; done > scalesia.contractions
for i in 217 288 386 396 426 436 466 484 575 617 624 659 884 1055 1113 1315 1362 1692 1954 2018 2070 2128 2626 3639 3742 4046; do echo "This is gene family number $i"; head -n $i ../../../01*/dump.blast_output.mci.I30 | tail -1; done > scalesia.expansions
```

03_HeliaGenes
```
cd 03_HeliaGenes
## Processing Helianthus because if a family is missing from scalesia, how can we really know the GOs?
cat ../02_family/scalesia.contractions | tr "\t" "\n" | grep "HELIA" > ./Helia.ScalesiaContractions.tsv

grep "$j" ../../../../../02_otherGenomes/01_Helianthus_annuus/Ha412v1r1_genes_v1.0.gff3 | grep "mRNA"; done < Helia.ScalesiaContractions.tsv > Helia.ScalesiaContractions.fullGeneInfo.tsv
grep "Ontology_term" Helia.ScalesiaContractions.fullGeneInfo.tsv | sed "s/.*Ontology_term\=//" | tr "," "\n" > HELIA.scalesiaContractions.GO.tsv
```


04_GOs
```
cd ../02_family
cat scalesia.contractions | tr "\t" "\n" | grep "SCALE" > ../04_GOs/gene.contractions.scalesia
cat scalesia.expansions | tr "\t" "\n" | grep "SCALE" > ../04_GOs/gene.expansions.scalesia

sed -i "s/\=/_/" *


while read gene; do grep $gene ../../../../../15_swissProt/07_TopGO_format/scalesia.tsv; done < gene.contractions.scalesia > contracted.GOs
while read gene; do grep $gene ../../../../../15_swissProt/07_TopGO_format/scalesia.tsv; done < gene.expansions.scalesia > expanded.GOs

awk -F "\t" '{print $2}' expanded.GOs | tr "," "\n" | sed "s/ //g" | grep "GO" > scalesia.expandedGOs.tsv
awk -F "\t" '{print $2}' contracted.GOs | tr "," "\n" | sed "s/ //g" > scalesia.contractedGOs.tsv
```

05_final_files

```
cd 05_final_files
ln -s ../04_GOs/*tsv .
ln -s ../03_HeliaGenes/HELIA.scalesiaContractions.GO.tsv .

## I get for contractions Helianthus and Scalesia GOs for one reason: gene families cotnracted in Scalesia have no genes ... right? So they won't show up.
```

06_FilesForTopGo
```
ln -s ../04_GOs/gene.* .
ln -s ../../../../../15_swissProt/07_TopGO_format/scalesia.tsv .
```

07_Revigo
```
I pushed the expansion BP file through revigo, keeping only SIGNIFICANT GOs
BP_Elim_Scalesia_CAFEexpansions.txt
```
