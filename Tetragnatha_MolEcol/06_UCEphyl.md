### 01 - First we get a consensus fasta from the alignments.

```
# set ref genome, bam and the output
REF=T_kauaiensis_ref.fasta
BAM=${SAMPLE}_dedup_mapQual_sorted.bam
FASTA=${SAMPLE}.fa

# call variants with bcftools
bcftools mpileup -a DP -Ou -f $REF $BAM | bcftools call -mv -O z -o ${SAMPLE}_calls.vcf.gz

# index
bcftools index ${SAMPLE}_calls.vcf.gz

# normalise
bcftools norm -f $REF ${SAMPLE}_calls.vcf.gz -O b -o ${SAMPLE}_calls.norm.bcf

# filter - need to set specific filters
bcftools filter --IndelGap 4 ${SAMPLE}_calls.norm.bcf \
-i '%QUAL>20 & FMT/DP>1 & FMT/DP<30' \
-O z -o ${SAMPLE}_call.norm.flt-indels.vcf.gz

# index
bcftools index ${SAMPLE}_call.norm.flt-indels.vcf.gz

# next create a consensus
cat $REF | bcftools consensus -s $SAMPLE ${SAMPLE}_call.norm.flt-indels.vcf.gz > $FASTA
gzip $FASTA
```

### 02 - now we use Phyluce

```
# Following the following tutorial:
# https://phyluce.readthedocs.io/en/latest/tutorial-three.html

# First we download the UCE files for spiders
wget 'https://ndownloader.figshare.com/files/6042078'

# Second, we need to create twobit files. That's the input for the program.
for i in *fa; do
	echo "### Working on $i";
	faToTwoBit $i ${i%.fa}.2bit;
	twoBitInfo  ${i%.fa}.2bit ${i%.fa}.tab;
	mkdir ../../02_phyluce/01_twoBitFiles/${i%.fa};
	mv ${i%.fa}.2bit ../../02_phyluce/01_twoBitFiles/${i%.fa};
	mv ${i%.fa}.tab ../../02_phyluce/01_twoBitFiles/${i%.fa};
done


# Third, we align the probes to the genomes.

# To get a scaffold list:
ls ../../01_bcftoolsSNPcall/results/*fa -l | sed "s/.*\/T_/T_/; s/\.fa//" | tr "\n" " "


phyluce_probe_run_multiple_lastzs_sqlite --db tetragna.sqlite  \
--output 02_latz  \
--scaffoldlist T_anuenue_046 T_anuenue_047 T_anuenue_048 T_anuenue_049 T_anuenue_060 T_anuenue_066 T_brevignatha_007 T_brevignatha_050 T_brevignatha_062 T_brevignatha_067 T_filiciphila_074 T_filiciphila_075 T_kamakou_001 T_kamakou_008 T_kamakou_009 T_kamakou_010 T_kamakou_011 T_kamakou_012 T_kamakou_013 T_kamakou_014 T_kamakou_015 T_kamakou_068 T_kauaiensis_016 T_kauaiensis_017 T_kauaiensis_018 T_kauaiensis_019 T_kauaiensis_020 T_kikokiko_072 T_macracantha_073 T_macracantha_076 T_mohihi_021 T_mohihi_022 T_mohihi_023 T_mohihi_051 T_mohihi_052 T_perreirai_024 T_perreirai_025 T_perreirai_026 T_perreirai_027 T_perreirai_061 T_pilosa_028 T_pilosa_029 T_pilosa_030 T_pilosa_031 T_pilosa_032 T_polychromata_053 T_polychromata_054 T_polychromata_055 T_polychromata_056 T_polychromata_057 T_quasimodo_002 T_quasimodo_003 T_quasimodo_004 T_quasimodo_005 T_quasimodo_006 T_quasimodo_033 T_quasimodo_034 T_quasimodo_035 T_quasimodo_036 T_quasimodo_058 T_quasimodo_063 T_quasimodo_064 T_quasimodo_065 T_restricta_037 T_restricta_038 T_restricta_069 T_restricta_070 T_restricta_071 T_tantalus_039 T_tantalus_040 T_waikaimoi_059 T_waikamoi_041 T_waikamoi_042 T_waikamoi_043 T_waikamoi_044 T_waikamoi_045 \
--genome-base-path ./01_twoBitFiles  \
--probefile /data/bigexpansion/jcerca/010_Tetragnatha/05_PhylogenUCE/02_phyluce/00_UCEs/Arachnida-UCE-1.1K-v1.fasta \
--cores 80

# NOTICE: run without \
# NOTICE: it doesn't like "-" in the names..

# Now we have to make a genomes.conf file as specified in the tutorial
nano genomes.conf

cat genomes.conf
# [scaffolds]
# CVH07_1:/data/bigexpansion/jcerca/011_SparrowPhylogeny/03_UCEs/02_twobitConvert/CVH07_1/CVH07_1.2bit
# CVH09_1:/data/bigexpansion/jcerca/011_SparrowPhylogeny/03_UCEs/02_twobitConvert/CVH09_1/CVH09_1.2bit
# ...

phyluce_probe_slice_sequence_from_genomes \
    --lastz 02_latz \
    --conf genome.conf \
    --flank 100 \
    --name-pattern "Arachnida-UCE-1.1K-v1.fasta_v_{}.lastz.clean" \
    --output 03_fasta

## Run without the \

#Now we transition to the second part of the tutorial
#  https://phyluce.readthedocs.io/en/latest/tutorial-one.html#uceextraction
# But it's actually confusing, so now I am doing this on my own:

cd 04_namescCleaned/
#First, we clean names so we get uces
cp ../03*/* .
for i in *fasta; do sed "s/>.*uce/>uce/" $i | sed "s/|.*//" > ${i%.fasta}.corrected.fa; done
rm *fasta

### Still on 05, but sending output to 06_deinterleavingFastas
for i in *fa; do awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $i > ../05_deinterleavingFastas/${i%.corrected.fa}.deinterleaved.fa; done

cd 05_deinterleavingFastas
## Now we get a list of uces present in more than 50 individuals
grep ">" * | awk -F ":" '{print $2}' | sort | uniq -c | awk '$1>10' | awk '{print $2}' | sed "s/>//" > ../uce_list_moreThan10Inds.tsv


# code explained
# grep ">" * | ### grep ">" in all files on 06_deinterleavingFastas
# awk -F ":" '{print $2}' | ### print the second field, i.e. first field will be name, second the uce
# sort | uniq -c | ### Sort, unique them and count.
# awk '$1>50' | awk '{print $2}' | ### print only those uces present in 50 individuals or more, and print the uce list
# sed "s/>//" > ../uce_list_moreThan50Inds.tsv ### remove the ">" that may be annoying. Get this on a file.

###
cd 06_separatingTheData
# Now we get all the uce sets to separate files:
while read uce; do echo "This is my uce $uce"; grep -A 1 --no-group-separator "$uce$" 05_deinterleavingFastas/* > ./06_separatingTheData/$uce.fasta; done  < uce_list_moreThan10Inds.tsv

# code explained
# while read uce; do ### opening the while loop
#echo "This is my uce $uce"; ### Just to see how it progresses.
# grep -A 1 --no-group-separator "$uce$" 06_deinterleavingFastas/* > ./07_separatingTheData/$uce.fasta; ### get the $uce ... LAST $ IS VERY IMPORTANT BECAUSE grep uce-10 would get i.e. uce-100 uce-1013. The $ limits this
# done  < uce_list_moreThan50Inds.tsv #list of UCEs

###
cd 07_cleaningUptheheaders
#Now we need to clean up the headers.
#on 07_* run the following
for i in *fasta; do sed -e " />/! s/.*fa//" $i | sed "s/.*\//>/; s/\.deinter.*//" > ../07_cleaningUptheheaders/${i%.fasta}.headerCleaned.fasta; done

#code explained
# for i in *fasta; do### open loop
# sed -e " />/! s/.*fa//" $i | # on lines without ">", remove everything before fa (i.e. clean everything before the sequence)
# sed "s/.*\//>/; s/\.deinter.*//" ../08_cleaningUptheheaders/${i%.fasta}.headerCleaned.fasta  # On the remaining header clean everything so we only have the individual ID

```
### 03 - Alingment and tree
```
for i in *fasta; do echo "They see me workin $i"; mafft --thread 20 --auto $i > ../08_alignment/${i%.fasta}.fai; done


##
cd 08_alignment

## Remember Fasconcat needs files to be called .fas
for i in *fai; do mv $i ${i%.fai}.fas; done
perl ../../../../local_bin/fasconcat/FASconCAT-G_v1.04.pl -s -p -n
mv FcC_* ../09_concatenation/

## making the partitions file
cat FcC_info.xls| grep "uce-" | awk '{print $1" = "$2"-"$3";"}' | sed 's/^/\tcharset /' | sed '1 i\begin sets;' | sed '1 i\#nexus' | sed '$ a\end;' > partitions.nex
```


###
cd 10_tree

ln -s ../09_concatenation/FcC_supermatrix.fas .
ln -s ../09_concatenation/partitions.nex .

iqtree -s FcC_supermatrix.fas -T 20 -B 1000 -p partitions.nex
```

```
