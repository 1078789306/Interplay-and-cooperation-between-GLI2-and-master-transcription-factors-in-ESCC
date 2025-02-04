bowtie2 -p 100 -x ./bowtie2index/hg19/hg19 -1 $file'_2_clean.fq.gz' -S $file'_Input.sam'
grep -v 'chrM' $file'_Input.sam' >$file'_Input_chrM.sam'
samtools view -@ 10 -b -S $file'_Input_chrM.sam' > $file'_Input.bam'
picard MarkDuplicates \INPUT=$file'_Input.bam' \OUTPUT=$file'_dechrM_picard.sort.bam' \METRICS_FILE= $file'.markDup.metric' REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true 
samtools sort -@ 10 $file'_dechrM_picard.sort.bam' -o $file'_picard.sort.bam'
samtools index $file'_picard.sort.bam'
macs2 callpeak -t $file'_picard.sort.bam' -f BAM -g hs --keep-dup all --name $file --bdg --broad --SPMR --nomodel --extsize 200 -q 0.01