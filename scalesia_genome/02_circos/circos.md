Folder structure
00_data 01_alignment 02_circos

Example with MASKED subgenomes. The genome is too repeat-rich to obtain results.

00_data
```
conda activate circos
### First, we get only the chromossomes for the Scalesia subgenomes. I have the code in data:
cd 00_data

ln -s ../../../01_annotation/09_masked_subgenomes/*fasta .
ls
# masked_scalesia_subgenomeA.fasta masked_scalesia_subgenomeB.fasta
```

We use MUMmer to align.

01_alignment
```
# promer will run the alignment
promer -p Scalesia_maskedsubgenomes.promer_def ../00_data/scalesia_subgenomeA.fa \
../00_data/scalesia_subgenomeB.fa &> mummer.log

# and show-coords will make the format reaable.
show-coords -r Scalesia_maskedsubgenomes.promer_def.delta > Scalesia_maskedsubgenomes.promer_def.coords

# Cleaning the ";"" and "="" from the scalesia genome because circos doesn't like it.
sed -i "s/_[A-Z].*=[0-9]*\t/\t/" Scalesia_maskedsubgenomes.promer_def.coords
sed -i "s/_[A-Z].*=[0-9]*//" Scalesia_maskedsubgenomes.promer_def.coords
```

Now, the plotting.
02_circos

```
# First, we start by making karyotypes.
cd 02_circos/config_files

nano 01_MakingKaryotes.py
# Change indir to your input directory
# Change outdir to your out directory
# Doublecheck that the fasta files end up with "fasta", and not "fas", or "fa"

python 01_MakingKaryotes.py

Output folder: cd 02_circos/01_karyotype

# the file should look something like this (and not have any "=" or ";" characters ):
head *txt -n 3
chr - ScDrF4C_27 1 0  70807918 green
chr - ScDrF4C_6 2 0  73591239 green
chr - ScDrF4C_1633 3 0  80823814 green

# To clean it, I did:
sed -i "s/_[A-Z].*=[0-9]*//" *txt

# Second, we draw the links.
cd 02_circos/02_links

nano 02_links.py
# Change coords_fname to the file we have created with mummer
# Change links_fname (output)
python 02_links.py

Output folder: cd 02_circos/02_links

## Third, we need to filter the links now.
cd 02_circos/config_files

nano ./03_filterLinks.py
# Change coords_fname
# Change min_length (min length 10,000 = did not work well)
# Change links_fname (output)

python ./03_filterLinks.py

# output: cd 02_circos/03_Filteredlinks
#  24288 filteredLinks_ScalA_vs_ScalB.min_length1250bp.txt
#  12658 filteredLinks_ScalA_vs_ScalB.min_length1500bp.txt
#   3947 filteredLinks_ScalA_vs_ScalB.min_length2000bp.txt
#     49 filteredLinks_ScalA_vs_ScalB.min_length5000bp.txt

### Four, we add colours.
cd 02_circos/config_files

nano 04_colour.py
# Remember to change the variables inside inside the script.
# You can get the scaffold IDs by just looking on the 02_circos/01_karyotype/*txt files, pick one of the comparisons to do.
chr2color = {}
chr2color['ScDrF4C_6'] = 'lblue'
chr2color['ScDrF4C_1633'] = 'blue'
chr2color['ScDrF4C_1634'] = 'lred'
chr2color['ScDrF4C_9'] = 'red'
chr2color['ScDrF4C_4'] = 'lyellow'
chr2color['ScDrF4C_25'] = 'yellow'
chr2color['ScDrF4C_21'] = 'lpurple'
chr2color['ScDrF4C_19'] = 'purple'
chr2color['ScDrF4C_10'] = 'lgreen'
chr2color['ScDrF4C_17'] = 'green'
chr2color['ScDrF4C_14'] = 'lgrey'
chr2color['ScDrF4C_24'] = 'grey'
chr2color['ScDrF4C_30'] = 'lorange'
chr2color['ScDrF4C_5'] = 'orange'
chr2color['ScDrF4C_18'] = 'black'
chr2color['ScDrF4C_15'] = 'lgrey'
chr2color['ScDrF4C_16'] = 'grey'

# also change fasta_fname
# out_fname
# coords_fname
# links_fname
# and to change the size of the links you want. In this case, we're going with 5000 bp.
python 04_colour.py
output: cd 02_circos/04_coloring



cd 02_circos/04_coloring
sed -i "s/;[A-Z].*=[0-9]*//" *txt

sed "s/;[A-Z].*=[0-9]*//"
### Fifth, we make the output.
cd 02_circos/config_files

nano 05_circosConfigFile.config
#change..
# karyotype
# file (under link)
#
circos -conf 05_circosConfigFile.config

### Notice, for the final plotting I changed the colour coding files.
# this one: masked_scalesia_subgenomeA.color-coded.txt
# By adding: "chr number after chr - ScDrF4C_4 4 (this last 4 is what will be plotted)
chr - ScDrF4C_4 4 0  90640368 lyellow
# I also changed the order so subgenomes would come together :)
`
