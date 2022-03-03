sampleID = $1
echo -n -e  $sampleID"\t" > $sampleID_rate.list
samtools view -f 1 -F 12 $sampleID.sort.bam |cut -f 1|sort|uniq|wc -l >> $sampleID_rate.list
