The first thing we did was to identify adapters using Adapter Removal:

```
# F=Forward read
# R=Reverse read
AdapterRemoval --file1 $in/$F --file2 $in/$R --identify-adapters --threads 4
```

After identifying adapters, we used trimmomatic, feeding it with the specific adapter information for each individual:

```
# F=Forward read
# R=Reverse read
# sample = name of sample
# ADAPTER = adapter information
java -jar trimmomatic-0.39.jar PE -phred33 -threads 5 ${F} ${R} \
${sample}_F.cleaned.fastq.gz ${sample}_F.cleaned_unpaired.fastq.gz \
${sample}_R.cleaned.fastq.gz ${sample}_R.cleaned_unpaired.fastq.gz \
ILLUMINACLIP:${ADAPTER}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
```
