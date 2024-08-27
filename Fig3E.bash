########H3K27ac denisty
###Fig 3E Solo
bedtools intersect -a GLI2_CUT-Tag_peaks.bed -b H3K27ac_ChIP-Seq.broadPeak.bed  > GLI2_H3K27ac.broadPeak.bed
bedtools intersect -a TP63_CUT-Tag_peaks.bed -b H3K27ac_ChIP-Seq.broadPeak.bed  > TP63_H3K27ac.broadPeak.bed
bedtools intersect -a RUNX1_CUT-Tag_peaks.bed -b H3K27ac_ChIP-Seq.broadPeak.bed  > RUNX1_H3K27ac.broadPeak.bed
cat GLI2_H3K27ac.broadPeak.bed  TP63_H3K27ac.broadPeak.bed RUXN1_H3K27ac.broadPeak.bed > Solo_H3K27ac.broadPeak.bed
sort -k1,1 -k2,2n  Solo_H3K27ac.broadPeak.bed >  Solo_H3K27ac.broadPeak.sort.bed
bedtools merge -i  Solo_H3K27ac.broadPeak.sort.bed >  Solo_H3K27ac.broadPeak.sort.merge.bed
###Fig 3E Dual
bedtools intersect -a GLI2_H3K27ac.broadPeak.bed -b TP63_H3K27ac.broadPeak.bed  > GLI2_TP63_H3K27ac.broadPeak.bed
bedtools intersect -a GLI2_H3K27ac.broadPeak.bed -b RUNX1_H3K27ac.broadPeak.bed  > GLI2_RUXN1_H3K27ac.broadPeak.bed
bedtools intersect -a TP63_H3K27ac.broadPeak.bed -b RUNX1_H3K27ac.broadPeak.bed  > TP63_RUXN1_H3K27ac.broadPeak.bed
cat GLI2_TP63_H3K27ac.broadPeak.bed  GLI2_RUXN1_H3K27ac.broadPeak.bed TP63_RUXN1_H3K27ac.broadPeak.bed > Dual_H3K27ac.broadPeak.bed
sort -k1,1 -k2,2n Dual_H3K27ac.broadPeak.bed > Dual_H3K27ac.broadPeak.sort.bed
bedtools merge -i Dual_H3K27ac.broadPeak.sort.bed > Dual_H3K27ac.broadPeak.sort.merge.bed
###Fig 3E Trio
bedtools intersect -a GLI2_TP63_H3K27ac.broadPeak.bed -b RUNX1_H3K27ac.broadPeak.bed  > GLI2_TP63_RUNX1_H3K27ac.broadPeak.bed
sort -k1,1 -k2,2n GLI2_TP63_RUNX1_H3K27ac.broadPeak.bed > Trio_H3K27ac.broadPeak.bed
bedtools merge -i Trio_H3K27ac.broadPeak.sort.bed > Trio_H3K27ac.broadPeak.sort.merge.bed

##########ATAC denisty
###Fig 3E Solo
bedtools intersect -a GLI2_CUT-Tag_peaks.bed -b ATAC_ChIP-Seq.broadPeak.bed  > GLI2_ATAC.broadPeak.bed
bedtools intersect -a TP63_CUT-Tag_peaks.bed -b ATAC_ChIP-Seq.broadPeak.bed  > TP63_ATAC.broadPeak.bed
bedtools intersect -a RUNX1_CUT-Tag_peaks.bed -b ATAC_ChIP-Seq.broadPeak.bed  > RUNX1_ATAC.broadPeak.bed
cat GLI2_ATAC.broadPeak.bed  TP63_ATAC.broadPeak.bed RUXN1_ATAC.broadPeak.bed > Solo_ATAC.broadPeak.bed
sort -k1,1 -k2,2n  Solo_ATAC.broadPeak.bed >  Solo_ATAC.broadPeak.sort.bed
bedtools merge -i  Solo_ATAC.broadPeak.sort.bed >  Solo_ATAC.broadPeak.sort.merge.bed
###Fig 3E Dual
bedtools intersect -a GLI2_ATAC.broadPeak.bed -b TP63_ATAC.broadPeak.bed  > GLI2_TP63_ATAC.broadPeak.bed
bedtools intersect -a GLI2_ATAC.broadPeak.bed -b RUNX1_ATAC.broadPeak.bed  > GLI2_RUXN1_ATAC.broadPeak.bed
bedtools intersect -a TP63_ATAC.broadPeak.bed -b RUNX1_ATAC.broadPeak.bed  > TP63_RUXN1_ATAC.broadPeak.bed
cat GLI2_TP63_ATAC.broadPeak.bed  GLI2_RUXN1_ATAC.broadPeak.bed TP63_RUXN1_ATAC.broadPeak.bed > Dual_ATAC.broadPeak.bed
sort -k1,1 -k2,2n Dual_ATAC.broadPeak.bed > Dual_ATAC.broadPeak.sort.bed
bedtools merge -i Dual_ATAC.broadPeak.sort.bed > Dual_ATAC.broadPeak.sort.merge.bed
###Fig 3E Trio
bedtools intersect -a GLI2_TP63_ATAC.broadPeak.bed -b RUNX1_ATAC.broadPeak.bed  > GLI2_TP63_RUNX1_ATAC.broadPeak.bed
sort -k1,1 -k2,2n GLI2_TP63_RUNX1_ATAC.broadPeak.bed > Trio_ATAC.broadPeak.bed
bedtools merge -i Trio_ATAC.broadPeak.sort.bed > Trio_ATAC.broadPeak.sort.merge.bed






