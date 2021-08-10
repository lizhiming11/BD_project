sample = $1
sam = $2
hg38 = $3

bowtie2 -p 20 -x $3 -1 $1.fq.gz -2 $1.fq.gz -S ./1.Cleandata/$2.sam --un-conc ./1.Cleandata/$1.rm.fq
gzip -f 1.Cleandata/$1.rm.1.fq 1.Cleandata/$1.rm.2.fq
