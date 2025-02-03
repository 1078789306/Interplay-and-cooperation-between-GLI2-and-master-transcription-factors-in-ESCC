### ChIP-Seq and CUT&Tag data ananlysis

#Sequence alignment
file="H3K27ac"
/home/jyxyyxy/bowtie2-2.4.5-linux-x86_64/bowtie2 -p 100 -x /home/jyxyyxy/bowtie2index/hg19/hg19 ./$file/$file'.fastq.gz' -S ./$file/$file'.sam'
samtools view -@ 10 -b -S ./$file/$file'.sam' > ./$file/$file'.bam'
samtools sort -@ 10 ./$file/$file'.bam' -o ./$file/$file'.sort.bam'
samtools index ./$file/$file'.sort.bam'
file="INPUT"
/home/jyxyyxy/bowtie2-2.4.5-linux-x86_64/bowtie2 -p 100 -x /home/jyxyyxy/bowtie2index/hg19/hg19 ./$file/$file'.fastq.gz' -S ./$file/$file'.sam'
samtools view -@ 10 -b -S ./$file/$file'.sam' > ./$file/$file'.bam'
samtools sort -@ 10 ./$file/$file'.bam' -o ./$file/$file'.sort.bam'
samtools index ./$file/$file'.sort.bam'
#Removing PCR duplicates
file="H3K27ac"
picard MarkDuplicates \INPUT=./$file/$file'.sort.bam' \OUTPUT=./$file/picard2/$file'_picard_.sort.bam' \METRICS_FILE= ./$file/$file'_markDup1.metric' REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true 
file="INPUT"
picard MarkDuplicates \INPUT=./$file/$file'.sort.bam' \OUTPUT=./$file/picard2/$file'_picard_.sort.bam' \METRICS_FILE= ./$file/$file'_markDup1.metric' REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true 
#MACS2 call peak
macs2 callpeak  -t ./$file/picard2/$file'_H3K27ac_picard_.sort.bam' -c ./$file/picard2/$file'INPUT_picard_.sort.bam' -f BAM -g hs --keep-dup all --name ./$file/peak2/$file --bdg --broad --SPMR --nomodel --extsize 200 -q 0.000001
awk '{print $1,$2,$3,$4,$5}' ./$file/peak2/$file'_peaks.broadPeak' > ./$file/peak2/$file'_broadPeak.bed'
sed -i 's/ /\t/g' ./$file/peak2/$file'_broadPeak.bed'
awk 'BEGIN{OFS=" "}{$5= "@"}{$6= $4}{$7 = "."}{$8 = "@"}{$9 = "@"}{print $1,$4,$5,$2,$3,$8,$7,$9,$6}' ./$file/peak2/$file'_broadPeak.bed' > ./$file/peak2/$file'_broadPeak.gff'
sed -i 's/@//g' ./$file/peak2/$file'_broadPeak.gff'
sed -i 's/ /\t/g' ./$file/peak2/$file'_broadPeak.gff'
