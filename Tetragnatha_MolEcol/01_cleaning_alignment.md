# 01 - We look at the quality of all fastq files using fastqc

F = Forward reads
R = Reverse reads

# did this for all 76 genomes
fastqc --nogroup -t 4 --outdir . $in/$F $in/$R

# 02 - We identify adapters for all the
AdapterRemoval --file1 $in/$F --file2 $in/$R --identify-adapters --threads 4

# 03 - and remove adapters
java -jar $TRIM PE -phred33 -threads 5 ${in}/${F} ${in}/${R} \
${sample}_F.cleaned.fastq.gz ${sample}_F.cleaned_unpaired.fastq.gz \
${sample}_R.cleaned.fastq.gz ${sample}_R.cleaned_unpaired.fastq.gz \
ILLUMINACLIP:${ADAPT_FOLDER}/${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

# 04 - now we align the reads
# Remember that now we have Forward and Reverse paired (both F and R were kept by trimmomatic), Forward unpaired (only F was kept by trimmomatic) and Reverse unpaired (only R was kept by trimmomatic).
# BWA needs an ID ... (sample name, e.g. T_anuenue_046)
ID=$PREFIX
# BWA needs a platform
PLATFORM=Illumina
# BWA needs a library number
LIBRARY=1
# extract PU data
PU_DATA=$(zcat $F | head -1 | cut -d ":" -f 3,4)

# construct read group
READGROUP="@RG\tID:${ID}\tPL:${PLATFORM}\tLB:${LIBRARY}\tSM:${PREFIX}\tPU:${PU_DATA}"

# First, we map paired ends (F and R paired)
echo "Aligning $PREFIX, paired-end"
bwa mem -M -t 10 \
-R $READGROUP \
$REF $FORWARD $REVERSE | samtools view -b | samtools sort -T ${PREFIX} > ${PREFIX}_pe.bam

# Then we map the unpaired, Forward
echo "Aligning $PREFIX, unpaired, forward"
bwa mem -M -t 10 \
-R $READGROUP \
$REF $FORWARD_UNPAIR | samtools view -b | samtools sort -T ${PREFIX} > ${PREFIX}_up_forward.bam

# Then we map the unpaired, Reverse
echo "Aligning $PREFIX, unpaired, reverse"
bwa mem -M -t 10 \
-R $READGROUP \
$REF $REVERSE_UNPAIR | samtools view -b | samtools sort -T ${PREFIX} > ${PREFIX}_up_reverse.bam

## merge paired and unpaired
echo "Merging $PREFIX bams"
samtools merge -rf ${PREFIX}_merge.bam \
${PREFIX}_pe.bam  ${PREFIX}_up_forward.bam ${PREFIX}_up_reverse.bam

# sort merged
echo "Sorting merged $PREFIX bam"
samtools sort -T ${PREFIX} -o ${PREFIX}_merge_sort.bam ${PREFIX}_merge.bam

# rm unsorted
rm ${PREFIX}_merge.bam


# 05 - Map quality stats now!
samtools flagstat ${PREFIX}.bam > ./mapqual/$PREFIX.MappingQuality.stats
