fq1 = $1
fq2 = $2
outfile = $3
sample = $4

SOAPnuke filter -l 20 -q 0.5 -n 0.1 -d -m 5 -Q 2 -G -1 $fq1 -2 $fq2 -o $outfile -C $sample_clean.1.fq.gz -D $sample_clean.2.fq.gz
