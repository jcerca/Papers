01_LTRdigest
```
# To install LTR digest we need to go to genome tools and wget the latest version
# when making make, we need to do
make cairo=no

# Following http://avrilomics.blogspot.com/2015/09/ltrharvest.html
~/genometools-1.6.1/bin/gt suffixerator -db scalesia_atractyloides.chrOnly.fasta \
-indexname scalesia_atractyloides.fasta -tis -suf -lcp -des -ssp -sds -dna

# I'd like to keep the sequence name.  ... but we have to add another flag to run LTR digest. Code to keep sequence names (not running the -tabout no flag)
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt ltrharvest -index scalesia_atractyloides.fasta -out LTR.fa -gff3 LTR.gff3  -seqids yes

# We need to clean the LTR files just obtained with LTR digest, and we will add some HMM profiles for this.
# From the blog post: "specify a list of pHMM files in HMMER2 format for domain search,
# e.g. pHMMs from Pfam (eg. HMM1.hmm, HMM2.hmm, HMM3.hmm). (...)
# the LTRdigest paper by Steinbiss et al (2009) gives a list of LTR retrotransposon Pfam domains in tables B1 and B2. For example, PF03732 = Retrotrans_gag; we can download the HMM from Pfam using:

# I checked the Table B1 from the supplementary from the Steinbiss et al and think it may be worth downloading:
# PF07253 - gypsy
# PF03732 - Retrotransposon gag protein
# PF00077 - RVP Retroviral aspartyl protease
# PF08284 - RVP_2 Retroviral aspartyl protease
# PF00078 - RVT_1 Reverse transcriptase
# PF07727 - RVT_2 Reverse transcriptase
# PF06817 - RVT_thumb Reverse transcriptase thumb domain
# PF06815 - RVT_connect Reverse transcriptase connection domain
# PF00075 - Rnase H Ribonuclease H domain
# PF00552 - Integrase Integrase DNA binding domain
# PF02022 - Integrase_Zn Integrase Zinc binding domain
# PF00665 - rve Integrase core domain
# PF00098 - zf-CCHC Zinc knuckle domain
# PF00385 - Chromo ‘chromo’ (CHRomatin Organisation MOdifier) domain
# PF01393 - Chromo_shadow Chromo shadow domain
# PF00692 - dUTPase dUTPase domain
# PF01021 - TYA TYA transposon protein
# PF03078 - ATHILA ATHILA ORF-1 family
# PF04094 - DUF390 Protein of unknown function (DUF390)
# PF08330 - DUF1723 Protein of unknown function (DUF1723)
# PF04195 - Transposase_28 Putative gypsy type transposon
# PF05380 - Peptidase_A17 Pao retrotransposon peptidase
# PF01140 - Gag_MA Matrix protein (MA), p15
# PF02337 - Gag_p10 Retroviral GAG p10 protein
# PF01141 - Gag_p12 Gag polyprotein, inner coat protein p12
# PF00607 - Gag_p24 gag gene protein p24 (core nucleocapsid protein)
# PF02093 - Gag_p30 Gag P30 core shell protein
# PF00692 - dUTPase dUTPase domain
# PF00077 - RVP Retroviral aspartyl protease
# PF00078 - RVT_1 Reverse transcriptase
# PF06817 - RVT_thumb Reverse transcriptase thumb domain
# PF00552 - Integrase Integrase DNA binding domain
# PF02022 - Integrase_Zn Integrase Zinc binding domain
# PF00665 - rve Integrase core domain
# PF00075 - Rnase H Ribonuclease H domain
# PF00429 - TLV_coat ENV polyprotein (coat polyprotein)
# PF08791 - Viral_env Viral envelope protein
# PF09590 - Env-gp36 Env-gp36 lentivirus glycoprotein
# PF03408 - Foamy_virus_ENV Foamy virus envelope protein
# PF00516 - GP120 Envelope glycoprotein GP120
# PF03056 - GP36 Env gp36 protein (HERV/MMTV type)
# PF00517 - GP41 Envelope Polyprotein GP41
# PF00098 - zf-CCHC Zinc knuckle domain

# On PFAM, I searched for LTR, http://pfam.xfam.org/search/keyword?query=LTR&submit=Submit ; I also searched for Gypsy and  Copia
# PF14223	Retrotran_gag_2
# PF14244	Retrotran_gag_3
# PF12384	Peptidase_A2B	Ty3 transposon peptidase
# PF14244	Retrotran_gag_3	gag-polypeptide of LTR copia-type
# PF14223	Retrotran_gag_2	gag-polypeptide of LTR copia-type

mkdir hmmProfiles; cd hmmProfiles
wget http://pfam.xfam.org/family/PF03732/hmm

for i in PF07253 PF03732 PF00077 PF08284 PF00078 PF07727 PF06817 PF06815 PF00075 PF00552 PF02022 PF00665 PF00098 PF00385 PF01393 PF00692 PF01021 PF03078 PF04094 PF08330 PF04195 PF05380 PF01140 PF02337 PF01141 PF00607 PF02093 PF00692 PF00077 PF00078 PF06817 PF00552 PF02022 PF00665 PF00075 PF00429 PF08791 PF09590 PF03408 PF00516 PF03056 PF00517 PF00098 PF14223 PF14244 PF12384 PF14244 PF14223 ;
do
	wget http://pfam.xfam.org/family/$i/hmm
done

### We need to convert these.

conda activate compgen
for i in *; do hmmconvert -2 $i > ${i}.cleaned; done
# Cleaning names
for i in *cleaned; do j=$(grep ACC $i | awk '{print $2}'); echo $j; mv $i $j.hmm; done

# We will also download the HMMs from GyDB:
# https://gydb.org/index.php/Collection_HMM
#
wget https://gydb.org/extensions/Collection/collection/db/GyDB_collection.zip
unzip GyDB_collection.zip
cp GyDB_collection/profiles/*hmm .
rm -fr GyDB_collection
```

02_LTRdigest
```
## OK now cleaning the files obtained in the run above with LTR digest. First, we sort the gff3:
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt gff3 -sort LTR.gff3 > LTR.sorted.gff3

### Now we use LTR digest to curate the LTR models
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt ltrdigest -hmms hmmProfiles/*hmm -aaout -outfileprefix scalesia_ltrdigest -seqfile scalesia_atractyloides.chrOnly.fasta -matchdescstart <  LTR.sorted.gff3 > Scalesia_ltrdigest_output.gff

# We will remove all LTR retrotransposon candidates that don't have any domain hit at all.
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt select -rule_files ../filter_protein_match.lua -- < Scalesia_ltrdigest_output.gff > Scalesia_ltrdigest_output.filtered.gff

wc -l *gff
   5076450 Scalesia_ltrdigest_output.filtered.gff
   5310601 Scalesia_ltrdigest_output.gff

cd ../01*/
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt extractfeat -type LTR_retrotransposon -matchdescstart -seqid -encseq ../01_run/scalesia_atractyloides.chrOnly.fasta ../02_LTRdigest/Scalesia_ltrdigest_output.filtered.gff > ../02_LTRdigest/LTR_completesequence.filtered.fa
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt extractfeat -type long_terminal_repeat -matchdescstart -seqid -encseq ../01_run/scalesia_atractyloides.chrOnly.fasta ../02_LTRdigest/Scalesia_ltrdigest_output.filtered.gff > ../02_LTRdigest/LTR_LTRsequenceOnly.filtered.fa

# complete transposon
LTR_completesequence.filtered.fa
# LTR sequence once
LTR_LTRsequenceOnly.filtered.fa
```

03_LTRharvestHELIA
```
# Following http://avrilomics.blogspot.com/2015/09/ltrharvest.html
ln -s ../../02_otherGenomes/01_Helianthus_annuus/Ha412HOv2.0-20181130.fasta .
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt suffixerator -db Ha412HOv2.0-20181130.fasta -indexname Ha412HOv2.0-20181130.fasta -tis -suf -lcp -des -ssp -sds -dna

/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt ltrharvest -index Ha412HOv2.0-20181130.fasta -out LTR.HELIA.fa -gff3 LTR.HELIA.gff3  -seqids yes



# LTR DIGEST for Helianthus, the outgroup fhr the analysis.
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt ltrdigest -hmms ../01_run/hmmProfiles/*hmm -aaout -outfileprefix HELIA_ltrdigest -seqfile Ha412HOv2.0-20181130.fasta  -matchdescstart <  LTR.HELIA.sorted.gff3 > HELIA_ltrdigest_output.gff
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt gff3 -sort LTR.HELIA.gff3 > LTR.HELIA.sorted.gff3
cd ../01*/
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt ltrdigest -hmms ../01_run/hmmProfiles/*hmm -aaout -outfileprefix HELIA_ltrdigest -seqfile Ha412HOv2.0-20181130.fasta  -matchdescstart <  LTR.HELIA.sorted.gff3 > HELIA_ltrdigest_output.gff

# We will remove all LTR retrotransposon candidates that don't have any domain hit at all.
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt select -rule_files ../filter_protein_match.lua -- < HELIA_ltrdigest_output.gff > HELIA_ltrdigest_output.filtered.gff

wc -l *gff
   4919698 HELIA_ltrdigest_output.filtered.gff
   5196841 HELIA_ltrdigest_output.gff


cd ../01*/
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt extractfeat -type LTR_retrotransposon -matchdescstart -seqid -encseq ../03_LTRharvestHELIA/Ha412HOv2.0-20181130.fasta ../04_LTRdigestHELIA/HELIA_ltrdigest_output.filtered.gff > ../04_LTRdigestHELIA/LTR_completesequenceHELIA.filtered.fa
/data/bigexpansion/jcerca/local_bin/genometools-1.6.1/bin/gt extractfeat -type long_terminal_repeat -matchdescstart -seqid -encseq ../03_LTRharvestHELIA/Ha412HOv2.0-20181130.fasta ../04_LTRdigestHELIA/HELIA_ltrdigest_output.filtered.gff > ../04_LTRdigestHELIA/LTR_LTRsequenceOnlyHELIA.filtered.fa

# complete transposon
LTR_completesequenceHELIA.filtered.fa
# LTR sequence once
LTR_LTRsequenceOnlyHELIA.filtered.fa
```

04_AssigningSubgenomes

```
# Separating both subgenomes
ln -sf ../02_LTRdigest/*fa .
grep -A 1 --no-group-separator "ScDrF4C_25]\|ScDrF4C_16]\|ScDrF4C_10]\|ScDrF4C_30]\|ScDrF4C_15]\|ScDrF4C_14]\|ScDrF4C_9]\|ScDrF4C_1633]\|ScDrF4C_1634]\|ScDrF4C_17]\|ScDrF4C_18]\|ScDrF4C_19]\|ScDrF4C_4]\|ScDrF4C_21]\|ScDrF4C_5]\|ScDrF4C_6]\|ScDrF4C_24]" LTR_LTRsequenceOnly.filtered.fa > LTRsequenceOnly_subgenomeA.fa
grep -A 1 --no-group-separator "ScDrF4C_25]\|ScDrF4C_16]\|ScDrF4C_10]\|ScDrF4C_30]\|ScDrF4C_15]\|ScDrF4C_14]\|ScDrF4C_9]\|ScDrF4C_1633]\|ScDrF4C_1634]\|ScDrF4C_17]\|ScDrF4C_18]\|ScDrF4C_19]\|ScDrF4C_4]\|ScDrF4C_21]\|ScDrF4C_5]\|ScDrF4C_6]\|ScDrF4C_24]" LTR_completesequence.filtered.fa > LTRcompleteSequence_subgenomeA.fa

grep -A 1 --no-group-separator "ScDrF4C_12]\|ScDrF4C_1]\|ScDrF4C_116]\|ScDrF4C_11]\|ScDrF4C_13]\|ScDrF4C_23]\|ScDrF4C_1632]\|ScDrF4C_20]\|ScDrF4C_28]\|ScDrF4C_26]\|ScDrF4C_8]\|ScDrF4C_27]\|ScDrF4C_2]\|ScDrF4C_22]\|ScDrF4C_29]\|ScDrF4C_7]\|ScDrF4C_3]" LTR_LTRsequenceOnly.filtered.fa > LTRsequenceOnly_subgenomeB.fa
grep -A 1 --no-group-separator "ScDrF4C_12]\|ScDrF4C_1]\|ScDrF4C_116]\|ScDrF4C_11]\|ScDrF4C_13]\|ScDrF4C_23]\|ScDrF4C_1632]\|ScDrF4C_20]\|ScDrF4C_28]\|ScDrF4C_26]\|ScDrF4C_8]\|ScDrF4C_27]\|ScDrF4C_2]\|ScDrF4C_22]\|ScDrF4C_29]\|ScDrF4C_7]\|ScDrF4C_3]" LTR_completesequence.filtered.fa > LTRcompleteSequence_subgenomeB.fa

ln -s ../02_LTRdigest/scalesia_ltrdigest_[35]ltr.fas .

grep -A 1 --no-group-separator "ScDrF4C_25\|ScDrF4C_16\|ScDrF4C_10\|ScDrF4C_30\|ScDrF4C_15\|ScDrF4C_14\|ScDrF4C_9\|ScDrF4C_1633\|ScDrF4C_1634\|ScDrF4C_17\|ScDrF4C_18\|ScDrF4C_19\|ScDrF4C_4\|ScDrF4C_21\|ScDrF4C_5\|ScDrF4C_6\|ScDrF4C_24" scalesia_ltrdigest_3ltr.fas > scalesia_3ltr_subgenomeA.fa
grep -A 1 --no-group-separator "ScDrF4C_25\|ScDrF4C_16\|ScDrF4C_10\|ScDrF4C_30\|ScDrF4C_15\|ScDrF4C_14\|ScDrF4C_9\|ScDrF4C_1633\|ScDrF4C_1634\|ScDrF4C_17\|ScDrF4C_18\|ScDrF4C_19\|ScDrF4C_4\|ScDrF4C_21\|ScDrF4C_5\|ScDrF4C_6\|ScDrF4C_24" scalesia_ltrdigest_5ltr.fas > scalesia_5ltr_subgenomeA.fa

grep -A 1 --no-group-separator "ScDrF4C_12\|ScDrF4C_1\|ScDrF4C_116\|ScDrF4C_11\|ScDrF4C_13\|ScDrF4C_23\|ScDrF4C_1632\|ScDrF4C_20\|ScDrF4C_28\|ScDrF4C_26\|ScDrF4C_8\|ScDrF4C_27\|ScDrF4C_2\|ScDrF4C_22\|ScDrF4C_29\|ScDrF4C_7\|ScDrF4C_3" scalesia_ltrdigest_3ltr.fas > scalesia_3ltr_subgenomeB.fa
grep -A 1 --no-group-separator "ScDrF4C_12\|ScDrF4C_1\|ScDrF4C_116\|ScDrF4C_11\|ScDrF4C_13\|ScDrF4C_23\|ScDrF4C_1632\|ScDrF4C_20\|ScDrF4C_28\|ScDrF4C_26\|ScDrF4C_8\|ScDrF4C_27\|ScDrF4C_2\|ScDrF4C_22\|ScDrF4C_29\|ScDrF4C_7\|ScDrF4C_3]" scalesia_ltrdigest_5ltr.fas > scalesia_5ltr_subgenomeB.fa

```
05_Orthofinder

```
### I was experimenting and did four separate runs:
01 - using the full, filtered transposon sequence
02 - using only the full, filtered transposon long-repeat-terminal (the LTR)
03 - using only the transposon long-repeat-terminal (the LTR), 5-line, before filtering
04 - using only the transposon long-repeat-terminal (the LTR), 3-line, before filtering

# Before, I went through every folder and cleaned names.
for i in *fa; do sed -i  "s/\s\[/_/g; s/\s/_/g;  s/\]//g" $i; done (otherwise we will have a problem with the names)

05_Orthofinder/01_fulltransposon
ln -s ../../05_AssigningSubgenomes/LTRcompleteSequence_subgenome* .
ln -s ../../04_LTRdigestHELIA/LTR_completesequenceHELIA.filtered.fa .

orthofinder -f . -d -t 30 &> runlog.tsv

05_Orthofinder/02_onlyLTRseq/

ln -s ../../05_AssigningSubgenomes/LTRsequenceOnly_subgenome* .
ln -s ../../04_LTRdigestHELIA/LTR_LTRsequenceOnlyHELIA.filtered.fa .

orthofinder -f . -d -t 30

05_Orthofinder/03_ltr5/
ln -s ../../05_AssigningSubgenomes/*5ltr*fa .
ln -s ../../04_LTRdigestHELIA/*5ltr.fas .

orthofinder -f . -d -a 20 -t 20 &> runlog.tsv

# LTR sequence once

05_Orthofinder/04_ltr3/
ln -s ../../05_AssigningSubgenomes/*3ltr*fa .
ln -s ../../04_LTRdigestHELIA/*3ltr.fas .

orthofinder -f . -d -a 20 -t 20 &> runlog.tsv

```
07_parsingOutResults
```
# I decided to go with the OrthoFinder run 02 - using only the full, filtered transposon long-repeat-terminal (the LTR).
# Sorting the output for more than 75 genes, and more than 0 in Helia/Subgenome A or B
awk '$5>75 && $2>0 && $3 > 0 && $4 > 0' Orthogroups.GeneCount.tsv | less
awk '$5>75' Orthogroups.GeneCount.tsv | less


orthogroup      LTR_LTRsequenceOnlyHELIA.filtered       LTRsequenceOnly_subgenomeA      LTRsequenceOnly_subgenomeB      Total
OG0000007       2211    19      5       2235	# 1. SubA dominant
OG0000014       1805    6       37      1848	# 2. SubB dominant
OG0000041       716     18      512     1246	# 2. SubB dominant
OG0000042       1       1061    184     1246	# 1. SubA dominant
OG0000043       0       485     713     1198	# 3. (arguably equal)
OG0000019       0       839     869     1708	# 4. Equal numbers
OG0000093       0       352     416     768		# 4. Equal numbers
OG0000106       652     14      7       673		# 2. SubB dominant
OG0000120       181     50      336     567		# 2. SubB dominant
OG0000127       468     0       65      533		# 2. SubB specific
OG0000132       0       283     243     526		# 4. Equal numbers
OG0000135       1       365     149     515		# 2. SubB dominant
OG0000142       4       0       485     489		# 2. SubB dominant
OG0000143       328     62      95      485		# 4. Equal numbers - in all genomes
OG0000153       227     76      138     441		# 2. SubB dominant - in all genomes
OG0000165       380     34      0       414		# 1. SubA dominant
OG0000167       0       399     2       401		# 1. SubA dominant
OG0000173       0       182     200     382		# 4. Equal numbers
OG0000181       300     56      0       356		# 1. SubA specific
OG0000187       110     230     1       341		# 1. SubA dominant
OG0000188       1       340     0       341		# 1. SubA specific
OG0000190       322     7       2       331		# 1. SubA dominant - in all genomes
OG0000194       0       148     177     325		# 4. Equal numbers
OG0000197       92      200     31      323		# 1. SubA dominant - in all genomes
OG0000202       6       120     189     315		# 3. (arguably equal) - in all genomes
OG0000211       61      208     14      283		# 1. SubA dominant - in all genomes
OG0000213       231     4       44      279		# 2. SubB dominant - in all genomes
OG0000217       7       263     3       273		# 1. SubA dominant - in all genomes
OG0000223       201     39      21      261		# 4. Equal numbers - in all genomes
OG0000232       91      99      53      243		# 1. SubA dominant - in all genomes
OG0000233       0       146     97      243		# 4. Equal numbers
OG0000237       178     25      29      232		# 4. Equal numbers - in all genomes
OG0000241       220     2       5       227		# 4. subA dominant - in all genomes
OG0000246       68      66      90      224		# 4. Equal numbers - in all genomes
OG0000247       10      209     4       223		# 1. SubA dominant - in all genomes
OG0000255       0       68      148     216		# 1. SubA dominant
OG0000256       83      93      39      215		# 1. SubA dominant - in all genomes
OG0000260       114     61      36      211		# 1. SubA dominant - in all genomes
OG0000261       0       111     99      210		# 4. Equal numbers
OG0000262       0       65      143     208		# 1. SubA dominant
OG0000270       3       198     0       201		# 1. SubA dominant
OG0000274       13      0       184     197		# 3. SubB dominant
OG0000276       88      93      15      196		# 1. SubA dominant - in all genomes
OG0000283       155     0       29      184		# 2. SubB dominant
OG0000298       0       86      82      168		# 4. Equal numbers
OG0000307       0       77      81      158		# 4. Equal numbers
OG0000308       0       158     0       158		# 1. subA specific
OG0000310       23      12      122     157		# 2. subB dominant
OG0000313       0       78      78      156		# 4. Equal numbers
OG0000322       1       113     36      150		# 2. subB dominant
OG0000324       72      70      6       148		# 1. subA dominant
OG0000332       115     6       21      142		# 2. subB dominant
OG0000335       5       132     2       139		# 1. subA dominant - in all genome
OG0000338       96      29      12      137		# 1. subA dominant - in all genome
OG0000340       0       68      69      137		# 4. Equal numbers
OG0000346       128     2       1       131		# 1. Equal numbers - in all genome
OG0000348       10      51      70      131		# 4. Equal numbers - in all genome
OG0000350       1       129     0       130		# 1. subA specific
OG0000353       111     12      4       127		# 1. subA dominant - in all genome
OG0000357       34      34      58      126		# 3. (arguably equal) - in all genome
OG0000366       82      4       35      121		# 2. subB dominant - in all genome
OG0000371       39      34      46      119		# 4. Equal numbers - in all genome
OG0000380       115     0       2       117		# 1. subA dominant - in all genome
OG0000381       59      54      4       117		# 1. subA dominant - in all genome
OG0000392       0       2       110     112		# 2. subB dominant - in all genome
OG0000396       0       42      69      111		# 4. Equal numbers
OG0000401       70      15      24      109		# 4. Equal numbers - in all genome
OG0000402       91      6       12      109		# 2. subB dominant - in all genome
OG0000409       0       105     2       107		# 1. subA dominant
OG0000410       31      18      57      106		# 2. subB dominant - in all genome
OG0000411       1       105     0       106		# 1. subA dominant
OG0000416       0       53      52      105		# 4. Equal numbers
OG0000421       10      0       92      102		# 2. subB specific
OG0000427       93      3       2       98		# 4. Equal numbers - in all genome
OG0000428       84      8       6       98		# 4. Equal numbers - in all genome
OG0000434       1       94      0       95		# 1. subA dominant
OG0000437       1       0       93      94		# 2. subB specific
OG0000443       0       90      2       92		# 1. subA dominant
OG0000457       0       1       86      87		# 2. subB dominant
OG0000460       53      33      0       86		# 1. subA specific
OG0000461       80      4       2       86		# 1. subA dominant
OG0000465       40      4       41      85		# 2. subB dominant - in all genome
OG0000470       0       60      24      84		# 1. subA dominant
OG0000472       2       81      0       83		# 1. subA specific
OG0000473       39      24      20      83		# 4. Equal numbers - in all genome
OG0000488       45      10      24      79		# 2. subB dominant - in all genome
OG0000489       37      7       35      79		# 2. subB dominant - in all genome
OG0000492       0       1       75      76		# 2. subB dominant


###
# DECISION: CHOICE OF ORTHOGROUPS. ORIGINALLY I HAD CHOSEN ONLY ORTHOGROUPS IN EQUAL NUMBERS ON BOTH SUBGENOMES
# HOWEVER, AFTER CLEANING FOR LENGTHS, MANY OF THESE GET BIASED ... SO I'LL GET ALSO THOSE THAT MAY, AT FIRST, BE UNBALANCED ...
# AND ONLY ACCESS IF THEY ARE UNBALANCED AFTER CLEANING EVERYTHING.

07_parsingOutResults/01_InBothSubgenomes_andAlsoPresentInHELIA/01_alignment

#added after understanding that cleaning may bias
OG0000120       181     50      336     567		# 2. SubB dominant
OG0000153       227     76      138     441		# 2. SubB dominant - in all genomes
OG0000197       92      200     31      323		# 1. SubA dominant - in all genomes
OG0000202       6       120     189     315		# 3. (arguably equal) - in all genomes
OG0000232       91      99      53      243		# 1. SubA dominant - in all genomes
OG0000256       83      93      39      215		# 1. SubA dominant - in all genomes
OG0000260       114     61      36      211		# 1. SubA dominant - in all genomes
OG0000276       88      93      15      196		# 1. SubA dominant - in all genomes
OG0000310       23      12      122     157		# 2. subB dominant
OG0000324       72      70      6       148		# 1. subA dominant
OG0000332       115     6       21      142		# 2. subB dominant
OG0000338       96      29      12      137		# 1. subA dominant - in all genome
OG0000357       34      34      58      126		# 3. (arguably equal) - in all genome
OG0000402       91      6       12      109		# 2. subB dominant - in all genome
OG0000410       31      18      57      106		# 2. subB dominant - in all genome
OG0000488       45      10      24      79		# 2. subB dominant - in all genome
OG0000489       37      7       35      79		# 2. subB dominant - in all genome

for i in OG0000120 OG0000153 OG0000197 OG0000202 OG0000232 OG0000256 OG0000260 OG0000276 OG0000310 OG0000324 OG0000332 OG0000338 OG0000357 OG0000402 OG0000410 OG0000488 OG0000489 ; do
	mafft  --thread 20 --auto ../../../06_Orthofinder/02_onlyLTRseq/OrthoFinder/Results_Apr30/Orthogroup_Sequences/$i.fa > $i.fai;
done

06_parsingOutResults/01_InBothSubgenomes_andAlsoPresentInHELIA/02_gblocks
# I pushed the alignments through the online server, and selecting less stringent selection parameters, i.e.:
#	allow smaller final blocks
#	allow gap positions within the final blocks
#	allow less strict flanking positions

for i in *fai; do awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i | sed "s/ //g" > ${i%.fai}.cleaned.fai; done

06_parsingOutResults/01_InBothSubgenomes_andAlsoPresentInHELIA/03_lengthTrimming
cd ../02*
sed "/>/! s/a//g; />/! s/t//g; />/! s/c//g; />/! s/g//g;" $i | awk '{print $0, length()}'
sed "/>/! s/a//g; />/! s/t//g; />/! s/c//g; />/! s/g//g;" $i | awk '{print $0, length()}' | less
sed "/>/! s/a//g; />/! s/t//g; />/! s/c//g; />/! s/g//g;" $i | awk '{print $0, length()}' | awk '{print $2}' | sort -g

# Length of seq (alignment length w/ aliview)

#			Lenth of sequence 		Nr of seq 			seq removed (<50% data)
# OG0000143 745 positions			485 sequences 		56 sequences (429 left)
# OG0000223 311 positions			261 sequences 		9  sequences (252 left)
# OG0000237 364 positions			232 sequences 		11 sequences (221 left)
# OG0000241 725 positions 			227 sequences 		26 sequences (201 left)
# OG0000246 586 positions 			224 sequences 		31 sequences (193 left)
# OG0000348 355 positions 			131 sequences 		3 sequences (128 left)
# OG0000371 691 positions 			119 sequences 		19 sequences(100 left)
# OG0000401 683 positions			109 sequences 		28 sequences (81 left)
# OG0000427 892 positions			98  sequences 		20 sequences (78 left)
# OG0000428 671 positions			98  sequences 		12 sequences (86 left)
# OG0000473 307 positions			83  sequences 		10 sequences (73 left)
# OG0000120 534 positions 			567 sequences 		145sequences (422 left)
# OG0000153 653 positions 			441 sequences 		159sequences (282 left)
# OG0000197 491 positions 			323 sequences 		48 sequences (275 left)
# OG0000202 269 positions 			315 sequences 		0  sequences (315 left)
# OG0000232 664 positions 			243 sequences 		33 sequences (210 left)
# OG0000256 579 positions 			215 sequences 		45 sequences (170 left)
# OG0000260 425 positions 			211 sequences 		11 sequences (200 left)
# OG0000276 586 positions 			196 sequences 		45 sequences (151 left)
# OG0000310 388 positions 			157 sequences 		23 sequences (134 left)
# OG0000324 656 positions 			148 sequences 		28 sequences (120 left)
# OG0000332 521 positions 			142 sequences 		20 sequences (122 left)
# OG0000338 229 positions 			137 sequences 		18 sequences (119 left)
# OG0000357 597 positions 			126 sequences 		15 sequences (109 left)
# OG0000402 378 positions 			109 sequences 		13 sequences (96 left)
# OG0000410 572 positions 			106 sequences 		16 sequences (90 left)
# OG0000488 556 positions 			79  sequences 		 4 sequences (75 left)
# OG0000489 560 positions 			79  sequences 		 0 sequences (79 left)


Important note: This step "took out" some of the balance between subgenomes. I.e. After removing the 31 sequences on OG0000246, no sequence was left from subgenome A

06_parsingOutResults/01_InBothSubgenomes_andAlsoPresentInHELIA/04_realign
### realigning
cd ../03*
for i in *fai; do
	mafft  --thread 10 --auto $i > ../04_realign/${i%.fai}.realign.fai;
done

06_parsingOutResults/01_InBothSubgenomes_andAlsoPresentInHELIA/05_tree

for i in *fai; do
	iqtree -s $i -T 10 -B 1000; done



#			Notes on Tree:
OG0000143; calculate JC distances on this one.
			# Scalesia monophyletic
			# the scalesia genes are on both subgenomes.
			# 2 HELIA have really long branches

# OG0000223; Biased subset of genes
			# Scalesia NOT monophyletic.
			# One of the subgenomes is over represented
			# 2 HELIA have really long branches

OG0000237; calculate JC distances on this one
			# Scalesia monophyletic
			# the scalesia genes are on both subgenomes.

# OG0000241;  Biased subset of genes
			# Scalesia NOT monophyletic.
			# One of the subgenomes is represented at all

# OG0000246;
			# Scalesia monophyletic
			# the scalesia genes are only on subgenome B.


# OG0000348;
			# Scalesia not monophyletic, but only by a single helianthus, which I think we can remove: (Ha 20694)

# OG0000371;
			# Scalesia monophyletic, but subgenomes separated - not really active on both subgenomes.

# OG0000401;
			# Scalesia not monophyletic
# OG0000427;
			# Only 3 scalesia. Also, Not monophyletic.

# OG0000428;
			# Scalesia not monophyletic
# OG0000473;
			# Scalesia monophyletic, but subgenomes segregated.

# OG0000120
			# Scalesia non monophyletic
			# both subgenomes separated - maybe a good example of branch lengths when they were separated.
			# If I show this as an example I should remove the 3 long branches from the Helianthus.

# OG0000153
			# Scalesia monophyletic
			# both subgenomes separated

# OG0000197
			# Scalesia non monophyletic
			# both subgenomes separated - maybe a good example of branch lengths when they were separated.

OG0000202; do JC distances
			# Scalesia monophyletic
			# But remove the long branch with the following IDs:
				# 74815; 74816; 165175; 116474; 169569
				# 165176; 21697; 21698; 116473;
				# 169570; 100148; 100147; 79956

# OG0000232
			# Scalesia non monophyletic
			# both subgenomes separated

# OG0000256
			# Scalesia non monophyletic
			# both subgenomes separated

# OG0000260
			# Scalesia monophyletic
			# both subgenomes separated

OG0000276; calculate JC distances
			# Be aware that this one is most present in one of the subgenomes.
			# Remove the long branch with 141925 and 141926

# OG0000310
			# Scalesia non monophyletic
			# both subgenomes mostly separated

# OG0000324
			# Scalesia non monophyletic
			# both subgenomes mostly separated

OG0000332; calculate JC distances
			#
# OG0000338
			# Scalesia non monophyletic
			# both subgenomes mostly separated

# OG0000357
			# Scalesia non monophyletic

# OG0000402
			# Scalesia non monophyletic

# OG0000410
			# Scalesia non monophyletic

# OG0000488
			# Scalesia non monophyletic
			# both subgenomes separated
OG0000489; calculate JC distances


####
06_parsingOutResults/01_InBothSubgenomes_andAlsoPresentInHELIA/06_cleaning
#			Notes on Tree:

OG0000143; calculate JC distances on this one.
			# 2 HELIA have really long branches
cp ../04*/OG0000143.gblocks.cleaned.realign.fai .
nano # removed 16943 and 16944

######
OG0000237; calculate JC distances on this one
nano # removed a long branch with 96518

######
OG0000202; do JC distances
nano # and removed the long branch with the following IDs:
				# 74815; 74816; 165175; 116474; 169569
				# 165176; 21697; 21698; 116473;
				# 169570; 100148; 100147; 79956

#####
OG0000276; calculate JC distances
			# Be aware that this one is most present in one of the subgenomes.
nano # remove the long branch with 141925 and 141926

#####
OG0000332; calculate JC distances
nano	# removed the long branch with 16332 / 163784 / 187157 / 187158

OG0000489; calculate JC distances

conda activate phyl
for i in *fai; do
	distmat -nucmethod 1 -sequence ${i} -outfile ${i%.final.fai}.distmat;
done

07_parsingOutResults/01_InBothSubgenomes_andAlsoPresentInHELIA/07_distMatONLYSCALESIA/
cd ../06*

# Making a distmatrix of only Scalesia values
for i in *fai; do
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i > ${i%.fai}.deinterleaved.fai
	mv ${i%.fai}.deinterleaved.fai ../07*/
	cd ../07*
	grep -A 1 "Sc" ${i%.fai}.deinterleaved.fai > ${i%.fai}.deinterleaved.onlyscalesia.fai
	distmat -nucmethod 1 -sequence ${i%.fai}.deinterleaved.onlyscalesia.fai -outfile ${i%.fai}.distmat;
	cd ../06*;
done

## This one, I copy pasted the results on excel and extracted average

#### Making a distmatrix for every scalesia value VS all helianthus ...
07_parsingOutResults/01_InBothSubgenomes_andAlsoPresentInHELIA/08_scalesiaVShelia/

# The code below creates one file for EACH scalesia gene vs Helianthus. so we can get an estimate of Scalesia VS helianthus. I did it this way so we don't have just a single matrix and I go value-by-value.
for i in *fai; do
	echo $i
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i > ${i%.fai}.deinterleaved.fai
	mv ${i%.fai}.deinterleaved.fai ../08*/
	cd ../08*
	mkdir -p ${i%.final.fai}
	cd ${i%.final.fai}
	mv ../${i%.fai}.deinterleaved.fai .
	grep "Sc" ${i%.fai}.deinterleaved.fai > tmp
	grep -A 1 "Ha" ${i%.fai}.deinterleaved.fai > allHelia.fai
	while read gene; do
	grep -A 1 $gene ${i%.fai}.deinterleaved.fai > ${gene#>}.scalesia;
	cat ${gene#>}.scalesia allHelia.fai > ${gene#>}.fai
	rm *.scalesia
	distmat -nucmethod 1 -sequence ${gene#>}.fai -outfile ${gene#>}.distmat
	done < tmp
	rm tmp allHelia.fai ${i%.fai}.deinterleaved.fai
	cd ../../06*;
done

# To get the Scalesia one results (since Scalesia sits on row 1, it's only one row)
for i in *; do cd $i; grep "Sc" *> ../${i}.tsv; cd ..; done
for i in *tsv; do sed -i "s/.*0.00\t//" $i; done
for i in *tsv; do sed -i "s/^long.*//; /^$/d" $i; done # each tsv file will have scalesia divergence for each OG.




							OG000489	OG000332	OG0000276	OG0000237	OG0000202	OG0000143
Average Scalesia branch 			: 	7,31		14,05		6,24		4,60		8,57		14,98
Max	Scalesia branch				:	16,8		27,04		51,99		10,24		24,01		63,73
Average Scalesia vs Helia branch		:	16,65		24,25		9,63		13,20		40,26		20,94
Max	Scalesia vs Healia branch		:	22,88		39,85		53,36		23,26		81,15		88,92


							OG000489	OG000332	OG0000276	OG0000237	OG0000202	OG0000143
How proportionally big are the averages ?2,28x		1,73x		1,54x		2,87x		4,70x		1,40x
How proportionally big are the max values ?1,36x	1,47x		1,02x		2,27x		3,38x		1,40x
(e.g. for OG0000143: 20,94 / 14,98 = 1,40x)


OG000489:
If the divergence in Helianthus is = 16,65 JC, and it diverged 15,57 MYA
Then, the divergence in Scalesia is 7,31 JC and thus x (x = 6,83 MYA)

OG000332
If the divergence in Helianthus is = 24,25  JC, and it diverged 15,57 MYA
Then, the divergence in Scalesia is 14,05  JC and thus x (x = 9,02 MYA)

OG0000276
If the divergence in Helianthus is = 9,63  JC, and it diverged 15,57 MYA
Then, the divergence in Scalesia is 6,24  JC and thus x (x = 10,09 MYA)

OG0000237
If the divergence in Helianthus is = 13,20  JC, and it diverged 15,57 MYA
Then, the divergence in Scalesia is 4,60  JC and thus x (x = 5,43 MYA)

OG0000202 ### This one I had to remove a long branch with ~10 Scalesia samples - so maybe we remove?
If the divergence in Helianthus is = 40,26 JC, and it diverged 15,57 MYA
Then, the divergence in Scalesia is 8,57  JC and thus x (x = 3,31 MYA)

OG0000143
If the divergence in Helianthus is = 20,94  JC, and it diverged 15,57 MYA
Then, the divergence in Scalesia is 14,98  JC and thus x (x = 11,13 MYA)


####
09_histograms
cd ../08*
for i in *tsv; do cat $i | tr "\t" "\n" | sed "s/^ //" | grep -v "long" | sed "/^$/d" | sed 's/ //g' > ../09_histograms/${i%.tsv}.Scalesia_vs_helia.tsv; done

cd ../07*
for i in *distmat; do tail -n +9 $i | tr "\t" "\n" | grep -v "long" | sed '/^$/d' | sed 's/ //g' > ../09_histograms/${i%.final.distmat}.Scalesia_vs_scalesia.tsv; done


##### Further pruning:

OG0000143 # Looks good on histogram

Scalesia vs Scalesia 14.5 // Scalesia vs Helianthus 18.25

If the divergence in Helianthus is = 18,25 JC, and it diverged 6,14 MYA
Then, the divergence in Scalesia is 14,5 JC and thus x (x = 4.87 MYA)


OG0000202 # Looks weird on histogram - not plotted

OG0000237 # Looks good on histogram

Scalesia vs Scalesia 4.35 // Scalesia vs Helianthus 12.8

If the divergence in Helianthus is = 12,8 JC, and it diverged 6,14 MYA
Then, the divergence in Scalesia is 4,35 JC and thus x (x = 2,09 MYA)


OG0000276 # Looks good on histogram

Scalesia vs Scalesia 4.75 // Scalesia vs Helianthus 7.75

If the divergence in Helianthus is = 7,75 JC, and it diverged 6,14 MYA
Then, the divergence in Scalesia is 4,75 JC and thus x (x = 3,76 MYA)


OG000332  # Looks OK on histogram

Scalesia vs Scalesia 13.5 // Scalesia vs Helianthus 25.5

If the divergence in Helianthus is = 25,5 JC, and it diverged 6,14 MYA
Then, the divergence in Scalesia is 13,5 JC and thus x (x = 3,25 MYA)


OG000489 # Looks OK on histogram

Scalesia vs Scalesia 6.15 // Scalesia vs Helianthus 18.25

If the divergence in Helianthus is = 18,25 JC, and it diverged 6,14 MYA
Then, the divergence in Scalesia is 6.15 JC and thus x (x = 2,07 MYA)
