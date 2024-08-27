#step1
bedtools intersect -a H3K4me3_broadPeak.bed -b H3K27Ac_peaks.bed > promoter.bed
H3K4me1_broadPeak.bed -b H3K27Ac_peaks.bed > enhancer.bed
#step2
bedtools intersect -a GLI2_CUT-Tag_peaks.bed -b promoter.bed  > GLI2_CUT-Tag-promoter.bed
bedtools intersect -a GLI2_CUT-Tag_peaks.bed -b enhancer.bed > GLI2_CUT-Tag-enhancer.bed
#step3
computeMatrix reference-point \
             -S GLI2_CUT-Tag.bw\
			    TP63_CUT-Tag.bw \
				RUNX1_CUT-Tag.bw \
				H3K27Ac.bw \
				ATAC.bw \
             -R GLI2_CUT-Tag-promoter.bed\
			    GLI2_CUT-Tag-enhancer.bed \
             --referencePoint center \
             -a 10000 -b 10000 -out heatmap.gz
 
plotHeatmap \
-m heatmap.gz \
-out heatmap.pdf 
 
 