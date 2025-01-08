infile=$1 # e.g., /data/EGAD00001006421/EGAF00004691240/sortdonor-166_240114.bam
outfile=$2 # this shell script process unmapped bams, perform mapping and bqsr

ref= /disk/reference/genome_DNA/Homo_sapiens/hg38/hg38.fa 
prefix=`basename $infile .bam`
echo $prefix

samtools sort -n -o $prefix.sorted.unmapped.bam $infile
bedtools bamtofastq -i $prefix.sorted.unmapped.bam -fq $prefix.R1.fq -fq2 $prefix.R2.fq
bwa mem $ref $prefix.R1.fq $prefix.R2.fq -t 8 -k 30 -M > $prefix.mapped.sam
java -jar $EBROOTPICARD/picard.jar SortSam I=$prefix.mapped.sam O=$prefix.sorted.mapped.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT
java -jar $EBROOTPICARD/picard.jar BuildBamIndex I=$prefix.sorted.mapped.bam VALIDATION_STRINGENCY=LENIENT
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$prefix.sorted.mapped.bam O=$prefix.sorted.mapped.rmdup.bam R=$ref VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true M=$prefix.rmdup.stat
gatk BaseRecalibrator -R $ref -O $prefix.recal.table --known-sites /disk/database/gatk_bundle/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf --known-sites /disk/database/gatk_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf -I $prefix.sorted.mapped.rmdup.bam
gatk ApplyBQSR -R $ref -O $outfile --bqsr-recal-file $prefix.recal.table -I $prefix.sorted.mapped.rmdup.bam