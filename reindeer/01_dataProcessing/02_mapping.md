We used bwa to map the data and obtain statistics:

```
# OUTDIR= base directory for output files
# REF= reference genome
# INPUT_FOLDER= input folder
# FORWARD= forward (paired) reads
# REVERSE= reverse (paired) reads
# FORWARD_UNPAIR= forward (unpaired) reads - notice, trimmomatic separates paired from unpaired reads, hence we have four input files.
# REVERSE_UNPAIR= reverse (unpaired) reads - notice, trimmomatic separates paired from unpaired reads, hence we have four input files.
# PREFIX = prefix for each individual
# ID = id for each individual
# PLATFORM=Illumina
# LIBRARY= library used
PU_DATA=$(zcat $FORWARD | head -1 | cut -d ":" -f 3,4)

# construct read group
READGROUP="@RG\tID:${ID}\tPL:${PLATFORM}\tLB:${LIBRARY}\tSM:${PREFIX}\tPU:${PU_DATA}"

# First, we map the paired-ended reads
bwa mem -M -t 10 \
-R $READGROUP \
$REF $FORWARD $REVERSE | samtools view -b | samtools sort -T ${PREFIX} > ${PREFIX}_pe.bam

# Now we map unpaired, forward reads
echo "Aligning $PREFIX, unpaired, forward"
bwa mem -M -t 10 \
-R $READGROUP \
$REF $FORWARD_UNPAIR | samtools view -b | samtools sort -T ${PREFIX} > ${PREFIX}_up_forward.bam

# Now we map unpaired, reverse reads
echo "Aligning $PREFIX, unpaired, reverse"
bwa mem -M -t 10 \
-R $READGROUP \
$REF $REVERSE_UNPAIR | samtools view -b | samtools sort -T ${PREFIX} > ${PREFIX}_up_reverse.bam

## Now we merge paired and unpaired files
echo "Merging $PREFIX bams"
samtools merge -rf ${PREFIX}_merge.bam \
${PREFIX}_pe.bam  ${PREFIX}_up_forward.bam ${PREFIX}_up_reverse.bam

# We sort the merged file
echo "Sorting merged $PREFIX bam"
samtools sort -T ${PREFIX} -o ${PREFIX}_merge_sort.bam ${PREFIX}_merge.bam

# rm unsorted
rm ${PREFIX}_merge.bam

# We now calculate depth
samtools depth ${PREFIX}.bam  |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' >  $PREFIX.depth.stats

# We now calculate mapping quality
samtools flagstat ${PREFIX}.bam > $PREFIX.mapqual.stats
```
