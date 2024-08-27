#bam to bw
bamCompare -p 10 -b1 ~/TE5/TP63/TE5_TP63_picard_.sort.bam -b2 ~/TE5/TE5_Input_picard_.sort.bam --skipNAs --scaleFactorsMethod readCount --operation subtract --outFileFormat bedgraph -o ~/TE5/TP63/$file'.bedgraph' --extendReads 200

awk '{if($4<0){$4=0};print}' ~/TE5/TP63/$file'.bedgraph' > ~/TE5/TP63/$file'_move0.bedgraph'
totalreads= awk '{sum=sum+$4}END{print sum}' ~/TE5/TP63/$file'_move0.bedgraph'
awk -v totalreads=6.97793e+07 '{$4=$4/totalreads*1000000;print}' ~/TE5/TP63/$file'_move0.bedgraph' > ~/TE5/TP63/$file'_rpm.bedgraph'
sort -k1,1 -k2,2n ~/TE5/TP63/$file'_rpm.bedgraph' > ~/TE5/TP63/$file'_sort.bedgraph'

/home/jyxyyxy/publicNAS1/bedGraphToBigWig ~/TE5/TP63/$file'_sort.bedgraph' /home/jyxyyxy/genome/hg19/hg19.fa.fai ~/TE5/TP63/$file'.bw'


