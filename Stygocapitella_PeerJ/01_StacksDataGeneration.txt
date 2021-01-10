The first step of the pipeline involves cleaning the raw fastq files using Stacks, using process radtags which will demultiplex the illumina files, remove adapters and read the enzme overhang

```
PRT="/home/jose10/catchenlab/local/bin/process_radtags"

barcodes_lane01="/home/jose10/research/001_stygoRAD/info/barcodes_lane01.txt"
barcodes_lane02="/home/jose10/research/001_stygoRAD/info/barcodes_lane02.txt"

lane01_F="/home/jose10/research/001_stygoRAD/raw/01_lane/1-Jose-lib-18X-2XP-L1_S18_L005_R1_001.fastq.gz"
lane01_R="/home/jose10/research/001_stygoRAD/raw/01_lane/1-Jose-lib-18X-2XP-L1_S18_L005_R2_001.fastq.gz"

lane02_F="/home/jose10/research/001_stygoRAD/raw/02_lane/2-L2-18X-Amp-2XP_S19_L006_R1_001.fastq.gz"
lane02_R="/home/jose10/research/001_stygoRAD/raw/02_lane/2-L2-18X-Amp-2XP_S19_L006_R2_001.fastq.gz"

output_folder="/home/jose10/research/001_stygoRAD/cleaned"

$PRT -i gzfastq -1 $lane01_F -2 $lane01_R -b $barcodes_lane01 -o $output_folder --renz_1 pstI --renz_2 mseI -q -r -c
```

