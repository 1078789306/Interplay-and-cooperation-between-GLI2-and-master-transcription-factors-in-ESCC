#Identify super-enhancer
pwd="/home/jyxyyxy/publicNAS4/pjy"
file2="INPUT"
file="H3K27ac"
mkdir -p $pwd/ROSE/$file'rose'
python ROSE_main.py -g HG19 \
-i $pwd/peak2/$file'_peaks.broadPeak.gff' \
-r $pwd/picard2/$file'.picard.bam' \
-c $pwd/picard2/$file2'.picard.bam' \
-o $pwd/ROSE/$file'rose'/ -s 12500 -t 2000
#super-enhancer related genes annotation
python ROSE_geneMapper.py -g HG19 -w 100000 -w 100000 \
-i $pwd/ROSE/$file'rose'/$file'_peaks_AllEnhancers.table.txt'

python ROSE_geneMapper.py -g HG19 -w 100000 -w 100000 \
-i $pwd/ROSE/$file'rose'/$file'_peaks_SuperEnhancers.table.txt'