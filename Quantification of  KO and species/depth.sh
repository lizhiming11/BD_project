sampleID = $1
bwa mem -t 20 -R "@RG\tID:$1\tLB:$1\tPL:ILLUMINA\tPM:HISEQ\tSM:$sampleID" UHGP.fa $sampleID.rm.1.fq.gz $sampleID.rm.2.fq.gz > $sampleID.sam
/ldfssz1/ST_HEALTH/Aging/lizhiming/bin/bin/samtools view -bS -@ 20 $sampleID.sam > $sampleID.bam 
/ldfssz1/ST_HEALTH/Aging/lizhiming/bin/bin/samtools sort $sampleID.bam $sampleID.sort 
sh mapping_rate.sh $sampleID
jgi_summarize_bam_contig_depths --outputDepth $sampleID.depth $sampleID.sort.bam
