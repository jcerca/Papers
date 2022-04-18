Folder structure:

00_UCEset  01_twobitConvert  02_latz  03_fasta  04_results

00_UCEset

```
## I got the set based on a recommendation from Jen Mandel. It's on my e-mail but I can download it here:
wget https://github.com/Smithsonian/Compositae-COS-workflow/raw/master/COS_probes_phyluce.fasta
```

Brilliant tutorial:
# https://phyluce.readthedocs.io/en/latest/tutorial-three.html

```
conda activate phyl
## The program  does not take files that start with numbers. This may lead to an error.
```

01_twobitConvert
```
# First, we start by creating twobit files. That's the input for the program.
cd 01_twobitConvert/
# We get the genome
ln -s ../../../../013_ScalesiaGenome/00_ScalesiaGenome/00_fromDovetail/scalesia_atractyloides.fasta .

i=scalesia_atractyloides.fasta
faToTwoBit $i ${i%.fa}.2bit;
twoBitInfo ${i%.fa}.2bit ${i%.fa}.tab;

# Notice the files must be each on a folder by itself.. hence the command above. Command explained:
# faToTwoBit $i ${i%.fasta}.2bit;                         ## convert the file from fasta to 2bit
# twoBitInfo  ${i%.fasta}.2bit  ${i%.fasta}.tab;          ## get a tab file with the size of the chromosomes
```
02_latz

```
phyluce_probe_run_multiple_lastzs_sqlite --db scalesia.sqlite  \
--output ./../02_lastz  \
--scaffoldlist scalesia_atractyloides \
--genome-base-path ./01_twobitConvert \
--probefile /data/bigexpansion/jcerca/014_ScalesiaPhylogenomics/01_30x_data/06_phyluce/00_UCEset/COS_probes_phyluce.fasta \
--cores 10
```

03_fasta
```
# To obtain paralogs, I modified two lines of the phyluce_probe_slice_sequence_from_genomes.py script.

# Specifically, on line 132, I changed 1 to 2.
# dupe_set = set([uce for uce, contigs in matches.iteritems() if len(contigs) > 2])

# Then, on line 306 I modifed == 1 to just be > 0
#assert len(matches.keys()) > 0, "There are multiple UCE matches"

# I saved the pythonscript with a new name
### RUNNING THIS ONE WITH THE MODIFIED FILE!!!!!!!!
./phyluce_probe_slice_sequence_from_genomesMODIFIED \
    --lastz 02_lastz \
    --conf genomes.conf \
    --flank 500 \
    --name-pattern "COS_probes_phyluce.fasta_v_{}.lastz.clean" \
    --output 04_scalesia_genomes_fasta &> 04_scalesia_genomes_fasta.txt
    --output 04_scalesia_genomes_fasta &> 04_scalesia_genomes_fasta.txt

```
Now, I formated a list with the 34 biggest contigs (chromosomes), called contig.ids, containing:
contig, two dots, and contig name

03_fasta

```
less contig.ids
contig:ScDrF4C_1
contig:ScDrF4C_10
# ... etc

# and run a while loop to separate UCEs:
while read GENE; do sed "s/>.*contig/>contig/" scalesia_atractyloides.fasta \
| sed "s/slice.*uce/uce/" | sed "s/|match.*//" | grep "$GENE;" | sed "s/.*uce/uce/" | sort > $GENE.uces; done < contig.ids

#Then cleaned the names..
for i in contig*; do mv $i ${i#contig:}; done

# and do this one by one :)
for i in ScDrF*; do echo $i; comm -12 $i ScDrF4C_1.uces | wc -l; done
for i in ScDrF*; do echo $i; comm -12 $i ScDrF4C_10.uces | wc -l; done

#And formatted the pairwise table in Cerca et al supplementary to get the correspondence.

```
Conclusions:
Pairs (notice, this still has the old chromosome names):
pair01: ScDrF4C_12-ScDrF4C_25
pair02: 01-16
pair03: 116-10
pair04: 11-30
pair05: 15-13
pair06: 14-23
pair07: 9-1632
pair08: 20-1633
pair09: 28-1634
pair10: 26-17
pair11: 8-18
pair12: 27-19
pair13: 4-2
pair14: 21-22
pair15: 29-5
pair16: 6-7
pair17: 3-24
