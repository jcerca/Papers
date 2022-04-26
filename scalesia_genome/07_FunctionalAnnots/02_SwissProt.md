00_db - Downloading the database

```
# trembl db - the one which is uncurated - is just too big.
# sprot db - the one which is curated - we use this one
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_plants.dat.gz

# Retrieving the information we want:
cat uniprot_sprot_plants.dat | grep "^ID\|^DR[ ]\+GO;\|^OS" > headers
# Basically we want lines starting with ID (ID), DR (for the GOs), and OS (species)
cat uniprot_sprot_plants.dat | sed -n  "/^SQ/,/\/\//p" | grep -v "^//$" > sequences # This code parses anything starting with SQ and finishing in // (where the sequences are)

grep -c "^ID" uniprot_sprot_plants.dat # all entries
21320124
grep -c "^SQ" sequences
21320124

# Now we modify the headers so we get the go terms together with the
cat headers | sed "s/ID[ ]\+/>/" | sed '/^>/s/[ ]\+/_/g; /^>/s/;/_/g;  /^>/s/__/_/g;  /^>/s/_AA\./_AAlength/g;' | sed '/^OS/s/ (.*//g' | sed '/^OS/s/OS[ ]\+//' | sed '/^DR/s/.*GO; //' | sed '/^GO:/s/;.*//' | sed '/^GO/s/GO/ GO/' | awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' | sed '/^>/{N; s/\n/ /;}' | sed '/>/s/\.//g' > newheaders

### Code explained
cat headers | # open the file
sed "s/ID[ ]\+/>/" | # On lines starting with ID, remove the first spaces and replace it with >
sed '/^>/s/[ ]\+/_/g; /^>/s/;/_/g;  /^>/s/__/_/g;  /^>/s/_AA\./_AAlength/g;' | # On lines starting with > do remove further spaces ; and __ and replace them with a single undersore. Also, change AA to aalength
sed '/^OS/s/ (.*//g' | # On the line starting with OS (species line) remove anything after the first ( - because it is the common name
sed '/^OS/s/OS[ ]\+//' | # On the line starting with OS, remove OS and spaces
sed '/^DR/s/.*GO; //' | # On the line starting with DR, remove the first GO;
sed '/^GO:/s/;.*//' | # On the line now starting with GO: remove all the GO information (other than the GO ID)
sed '/^GO/s/GO/ GO/' | # On the line starting with GO, remove the space before GO
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' | # Deinterleave any line not starting with >
sed '/^>/{N; s/\n/ /;}' | # If the line starts with >, remove new line
sed '/>/s/\.//g' > newheaders # cleaning dots and saving to a file

grep -c ">" newheaders
43549

## Now we clean up the sequence file (still has spaces)

cat sequences | sed '/^SQ  /!s/ //g'  | sed '/^SQ  /s/SQ/>SQ/' | sed "/>/s/;.*//g" | awk '{if(NR==1) {print $0} else {if($0 ~ />SQ/) {print "\n"$0} else {printf $0}}}' > newsequences

cat sequences | # open up the file
sed '/^SQ  /!s/ //g' | # For lines not starting with SQ space, replace spaces
sed '/^SQ  /s/SQ/>SQ/' | #Lines starting with SQ space, add a >
sed "/>/s/;.*//g" | # clean up the header on lines starting with >
awk '{if(NR==1) {print $0} else {if($0 ~ />SQ/) {print "\n"$0} else {printf $0}}}' > newsequences ## indent

grep -c ">" newsequences
43549

grep -v "^>" newsequences > newnewsequences

# Code explained
paste -d "\n" newheaders newnewsequences > db.sprot.fa
grep -c "^>" db.sprot.fa
43549
grep -cv "^>" db.sprot.fa
43549


### Decision - I'll be only using the SPROT db. The other one is too big ... and uncurated.

```
01_scalesiaBLAST
```
makeblastdb -in db.sprot.fa -out sprot_DB -parse_seqids -dbtype prot

01_scalesiaBLAST
blastp -out scalesia_plant_Sprot -db ../00_db/sprot_DB -query SCALE.scalesia_atractyloides_cleanedHeaders.faa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -num_threads 40
awk '$11<1e-10' scalesia_plant_Sprot  | awk '{print $1,$2}' > scalesia_plant_Sprot_evaluebelow1e-10.tsv
```

02_parsingOutGOs

```
cd ../01*

## Doing subgenomes separately, and loose chromosomes ... - now we have a list.
grep --no-group-separator "ScDrF4C_25_\|ScDrF4C_16_\|ScDrF4C_10_\|ScDrF4C_30_\|ScDrF4C_15_\|ScDrF4C_14_\|ScDrF4C_9_\|ScDrF4C_1633_\|ScDrF4C_1634_\|ScDrF4C_17_\|ScDrF4C_18_\|ScDrF4C_19_\|ScDrF4C_4_\|ScDrF4C_21_\|ScDrF4C_5_\|ScDrF4C_6_\|ScDrF4C_24_" scalesia_plant_Sprot_evaluebelow1e-10.tsv > ../02*/scalesia_plant_Sprot_evaluebelow1e-10.subgenomeA.tsv
grep --no-group-separator "ScDrF4C_12_\|ScDrF4C_1_\|ScDrF4C_116_\|ScDrF4C_11_\|ScDrF4C_13_\|ScDrF4C_23_\|ScDrF4C_1632_\|ScDrF4C_20_\|ScDrF4C_28_\|ScDrF4C_26_\|ScDrF4C_8_\|ScDrF4C_27_\|ScDrF4C_2_\|ScDrF4C_22_\|ScDrF4C_29_\|ScDrF4C_7_\|ScDrF4C_3_" scalesia_plant_Sprot_evaluebelow1e-10.tsv > ../02*/scalesia_plant_Sprot_evaluebelow1e-10.subgenomeB.tsv
grep -v --no-group-separator "ScDrF4C_12_\|ScDrF4C_1_\|ScDrF4C_116_\|ScDrF4C_11_\|ScDrF4C_13_\|ScDrF4C_23_\|ScDrF4C_1632_\|ScDrF4C_20_\|ScDrF4C_28_\|ScDrF4C_26_\|ScDrF4C_8_\|ScDrF4C_27_\|ScDrF4C_2_\|ScDrF4C_22_\|ScDrF4C_29_\|ScDrF4C_7_\|ScDrF4C_3_\|ScDrF4C_25_\|ScDrF4C_16_\|ScDrF4C_10_\|ScDrF4C_30_\|ScDrF4C_15_\|ScDrF4C_14_\|ScDrF4C_9_\|ScDrF4C_1633_\|ScDrF4C_1634_\|ScDrF4C_17_\|ScDrF4C_18_\|ScDrF4C_19_\|ScDrF4C_4_\|ScDrF4C_21_\|ScDrF4C_5_\|ScDrF4C_6_\|ScDrF4C_24_" scalesia_plant_Sprot_evaluebelow1e-10.tsv > ../02*/scalesia_plant_Sprot_evaluebelow1e-10.looseScaffolds.tsv

while read gene; do i=$(echo $gene | awk '{print $1}'); j=$(echo $gene | awk '{print $2}'); k=$(grep $j ../00_db/db.sprot.fa); echo -e "$i \t $k"; done < scalesia_plant_Sprot_evaluebelow1e-10.subgenomeA.tsv > scalesiaWithGOs.subA.tsv
while read gene; do i=$(echo $gene | awk '{print $1}'); j=$(echo $gene | awk '{print $2}'); k=$(grep $j ../00_db/db.sprot.fa); echo -e "$i \t $k"; done < scalesia_plant_Sprot_evaluebelow1e-10.subgenomeB.tsv > scalesiaWithGOs.subB.tsv
while read gene; do i=$(echo $gene | awk '{print $1}'); j=$(echo $gene | awk '{print $2}'); k=$(grep $j ../00_db/db.sprot.fa); echo -e "$i \t $k"; done < scalesia_plant_Sprot_evaluebelow1e-10.looseScaffolds.tsv > scalesiaWithGOs.looseScaffolds.tsv

# These files still have some garbadge. Cleaning:
cat scalesiaWithGOs.looseScaffolds.tsv | sed -e 's/\t.*[a-z] //' | sed 's/[a-z]\+) //' | grep -v "thaliana\|commersonii\|mays\|pendula\|carota\|max\|ananassa" | sed "s/natansCCMP621))//" > scalesiaWithGOs.looseScaffolds.CLEANED.tsv > scalesiaWithGOs.looseScaffolds.CLEANED.tsv
# I had to remove thaliana and others which were remnants from the previous code
cat scalesiaWithGOs.subA.tsv | sed -e 's/\t.*[a-z] //' | sed 's/[a-z]\+) //' | sed "s/Arabidopsis thaliana //" | grep "GO" | sed "s/natansCCMP621)) //; s/Chrysosplenium americanum //; s/Solanum tuberosum //; s/cardiophyllumL//" | grep -v "^>" > scalesiaWithGOs.subA.CLEANED.tsv
cat scalesiaWithGOs.subB.tsv | sed -e 's/\t.*[a-z] //' | sed 's/[a-z]\+) //' | grep "GO" | sed "s/Arabidopsis thaliana //; s/natansCCMP621)) //; s/cardiophyllumL//; s/Chrysosplenium americanum //; s/Solanum tuberosum //" | grep -v "^>" > scalesiaWithGOs.subB.CLEANED.tsv

# I debugged the cove above by removing the first column (scalesia gene, and adding:
# awk '{print $2,$3,$4,$5}' | grep "[a-z]"
```
03_interproGOs_cleaned
```
cat ../../14_interpro/SCALE/SCALE.sanitizedHeaders.longestisoforms.faa.tsv  | awk '$9<1e-10' | grep "GO:" | awk -F "\t" '{print $1, $14}' | tr "|" " " > interproGOs.tsv
grep --no-group-separator "ScDrF4C_25_\|ScDrF4C_16_\|ScDrF4C_10_\|ScDrF4C_30_\|ScDrF4C_15_\|ScDrF4C_14_\|ScDrF4C_9_\|ScDrF4C_1633_\|ScDrF4C_1634_\|ScDrF4C_17_\|ScDrF4C_18_\|ScDrF4C_19_\|ScDrF4C_4_\|ScDrF4C_21_\|ScDrF4C_5_\|ScDrF4C_6_\|ScDrF4C_24_" interproGOs.tsv > scalesia_intrproGOs.subgenomeA.tsv
grep --no-group-separator "ScDrF4C_12_\|ScDrF4C_1_\|ScDrF4C_116_\|ScDrF4C_11_\|ScDrF4C_13_\|ScDrF4C_23_\|ScDrF4C_1632_\|ScDrF4C_20_\|ScDrF4C_28_\|ScDrF4C_26_\|ScDrF4C_8_\|ScDrF4C_27_\|ScDrF4C_2_\|ScDrF4C_22_\|ScDrF4C_29_\|ScDrF4C_7_\|ScDrF4C_3_" interproGOs.tsv > scalesia_intrproGOs.subgenomeB.tsv
grep -v --no-group-separator "ScDrF4C_12_\|ScDrF4C_1_\|ScDrF4C_116_\|ScDrF4C_11_\|ScDrF4C_13_\|ScDrF4C_23_\|ScDrF4C_1632_\|ScDrF4C_20_\|ScDrF4C_28_\|ScDrF4C_26_\|ScDrF4C_8_\|ScDrF4C_27_\|ScDrF4C_2_\|ScDrF4C_22_\|ScDrF4C_29_\|ScDrF4C_7_\|ScDrF4C_3_\|ScDrF4C_25_\|ScDrF4C_16_\|ScDrF4C_10_\|ScDrF4C_30_\|ScDrF4C_15_\|ScDrF4C_14_\|ScDrF4C_9_\|ScDrF4C_1633_\|ScDrF4C_1634_\|ScDrF4C_17_\|ScDrF4C_18_\|ScDrF4C_19_\|ScDrF4C_4_\|ScDrF4C_21_\|ScDrF4C_5_\|ScDrF4C_6_\|ScDrF4C_24_" interproGOs.tsv > scalesia_intrproGOs.looseScaffolds.tsv
```

04_genesSorted
```
cat 02_parsingOutGOs/scalesiaWithGOs.looseScaffolds.CLEANED.tsv 03_interproGOs_cleaned/scalesia_intrproGOs.looseScaffolds.tsv | awk '{print $1}' | sort | uniq > 04_genesSorted/scalesia.looseScaffolds.genelist
cat 02_parsingOutGOs/scalesiaWithGOs.subA.CLEANED.tsv 03_interproGOs_cleaned/scalesia_intrproGOs.subgenomeA.tsv | awk '{print $1}' | sort | uniq > 04_genesSorted/scalesia.subgenomeA.genelist
cat 02_parsingOutGOs/scalesiaWithGOs.subB.CLEANED.tsv 03_interproGOs_cleaned/scalesia_intrproGOs.subgenomeB.tsv | awk '{print $1}' | sort | uniq > 04_genesSorted/scalesia.subgenomeB.genelist
```


05_geneByGene
```
mkdir loose subA subB
cd loose
while read gene; do echo $gene; grep $gene ../../02_parsingOutGOs/scalesiaWithGOs.looseScaffolds.CLEANED.tsv >> $gene; grep $gene ../../03_interproGOs_cleaned/scalesia_intrproGOs.looseScaffolds.tsv >> $gene; mv $gene $(echo $gene | sed "s/=/_/"); done < ../../04_genesSorted/scalesia.looseScaffolds.genelist

cd subA
while read gene; do echo $gene; grep $gene ../../02_parsingOutGOs/scalesiaWithGOs.subA.CLEANED.tsv >> $gene; grep $gene ../../03_interproGOs_cleaned/scalesia_intrproGOs.subgenomeA.tsv >> $gene; mv $gene $(echo $gene | sed "s/=/_/"); done < ../../04_genesSorted/scalesia.subgenomeA.genelist
cd subB
while read gene; do echo $gene; grep $gene ../../02_parsingOutGOs/scalesiaWithGOs.subB.CLEANED.tsv >> $gene; grep $gene ../../03_interproGOs_cleaned/scalesia_intrproGOs.subgenomeB.tsv >> $gene; mv $gene $(echo $gene | sed "s/=/_/"); done < ../../04_genesSorted/scalesia.subgenomeB.genelist
```

06_GOsPerGENE
```
cd ../05*/loose
for i in *; do echo $i; cat $i | sed "s/.*mRNA\-1 //" | tr " " "\n" | sort -u | grep "\S" > ../../06_GOsPerGENE/loose/$i; done
for i in *; do echo $i; cat $i | sed "s/.*mRNA\-1 //" | tr " " "\n" | sort -u | grep "\S" > ../../06_GOsPerGENE/subA/$i; done
for i in *; do echo $i; cat $i | sed "s/.*mRNA\-1 //" | tr " " "\n" | sort -u | grep "\S" > ../../06_GOsPerGENE/subB/$i; done
```

07_TopGO_format
```
cd ../06*/loose
for i in *; do echo $i; j=$(cat $i | tr "\n" "," | sed 's/,/, /g; s/,$//'); echo $i $j >> ../../07_TopGO_format/loose.tsv; done
cd ../../07*;
sed -i 's/,$// ; s/\-mRNA\-1 /\-mRNA\-1\t/' loose.tsv

cd ../06*/subA
for i in *; do echo $i; j=$(cat $i | tr "\n" "," | sed 's/,/, /g; s/,$//'); echo $i $j >> ../../07_TopGO_format/subA.tsv; done
cd ../../07*;
sed -i 's/,$// ; s/\-mRNA\-1 /\-mRNA\-1\t/' subA.tsv


for i in *; do echo $i; j=$(cat $i | tr "\n" "," | sed 's/,/, /g; s/,$//'); echo $i $j >> ../../07_TopGO_format/subB.tsv; done
cd ../../07*
sed -i 's/,$// ; s/\-mRNA\-1 /\-mRNA\-1\t/' subB.tsv


cat subA.tsv subB.tsv loose.tsv > scalesia.tsv
```
