For this analysis we need the fasta generated with populations (see 01). I did this for every population (it's population genetics for a reason!) but will only demonstrate a single population.

Let's prepare the data first.
```
# data: populations.samples.fa

# First, we need  files separated by species (each species is a colour here:
for i in blue green purple; do
grep -A 1 --no-group-separator "$i" 00_raw/populations.samples.fa > all.$i.sequences.fa
done

for i in blue green purple; do mkdir $i/01_all_${i}_loci; mv $i/*fa; done


# 2. Now, we further split them on populations.
# S. subterranea (blue) populations:

color=blue
for i in musselburough ardtoe keitum morsum little ile_callot nairn glenancross hausstrand; do
grep -A 1 --no-group-separator "$i" ../01_all_${color}_loci/all.${color}.sequences.fa > ${i}.fa
done


# 3. Alright, now run the script to get a locus-by-locus

for i in musselburough ardtoe keitum morsum little ile_callot nairn glenancross hausstrand; do
cd $i
echo "Working on $i"
perl Split_STACKS_fasta_file.pl ${i}.fa
echo "Removing Alleles"
rm *Alle*; cd ..
done

# 4. Now we get the loci with least missing data.
# Number of specimens for each population:
# blue
#        ardtoe          7
#        glenancross     4
#        hausstrand      1               #Not done! Too few.
#        ile callot      3
#        keitum          6
#        little gru      1               #Not done! Too few.
#        morsum          3
#        musselbu        3
#        nairn           2               #Not done! Too few.


##Showing only for ardtoe. 7 individuals= 14 alleles max. 14 alleles = 0 missing data. 12 alleles = 2 missing data.
# blue ardtoe
grep -c ">" *fas | tr ":" "\t" | sort -n -k 2 | less

grep -c ">" *fas | tr ":" "\t" | sort -n -k 2 | awk '$2> 13' | wc -l
# 122 loci have 14 alleles! (no missing data). Let's get them.

grep -c ">" *fas | tr ":" "\t" | sort -n -k 2 | awk '$2> 13' | cut -f 1 > 122loci_with_NOleast_missingData.ardtoeBlue.tsv

# 5. now we extract them to a specific folder
while read i
do cp $i ../../04_lociwithleastmissingdata/ardtoe/
done < 122loci_with_NOleast_missingData.ardtoeBlue.tsv

```

After this you need to have them all in the same folder and run DNAsp whichh will generate Tajima's D, pi, and S.
