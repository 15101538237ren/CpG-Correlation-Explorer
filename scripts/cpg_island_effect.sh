BASE_DIR=/Users/emmanueldollinger/PycharmProjects/CpG-Correlation-Explorer
DATA_DIR=$BASE_DIR/data
Genomic_Regions=$DATA_DIR/Genomic_Context/Genomic_Regions
cd $DATA_DIR

K=$DATA_DIR/K.bed

CGI=$Genomic_Regions/CGI.bed
Promoter=$Genomic_Regions/Promoter.bed
Enhancer=$Genomic_Regions/Enhancer.bed
UTR5=$Genomic_Regions/5UTR.bed

INPUT=$K
OUTPUT1=$DATA_DIR/K_CGI_intersected.bed
bedtools intersect -a $INPUT -b $CGI -wa  -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50% | awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"1}' >$OUTPUT1
OUTPUT2=$DATA_DIR/K_CGI_non_intersected.bed
bedtools intersect -a $INPUT -b $CGI -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%| awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"0}' >$OUTPUT2

INPUT=$DATA_DIR/K_CGI.bed
cat $OUTPUT1 $OUTPUT2 | gsort -k 1,1 -k2,2n --parallel=8  -S 50% > $INPUT
rm $OUTPUT1 $OUTPUT2 

OUTPUT1=$DATA_DIR/K_Promoter_intersected.bed
bedtools intersect -a $INPUT -b $Promoter -wa  -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50% | awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"1}' >$OUTPUT1
OUTPUT2=$DATA_DIR/K_Promoter_non_intersected.bed
bedtools intersect -a $INPUT -b $Promoter -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%| awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"0}' >$OUTPUT2

INPUT=$DATA_DIR/K_CGI_Prom.bed
cat $OUTPUT1 $OUTPUT2 | gsort -k 1,1 -k2,2n --parallel=8  -S 50% > $INPUT
rm $OUTPUT1 $OUTPUT2

OUTPUT1=$DATA_DIR/K_Enhancer_intersected.bed
bedtools intersect -a $INPUT -b $Enhancer -wa  -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50% | awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"1}' >$OUTPUT1
OUTPUT2=$DATA_DIR/K_Enhancer_non_intersected.bed
bedtools intersect -a $INPUT -b $Enhancer -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%| awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"0}' >$OUTPUT2

INPUT=$DATA_DIR/K_CGI_Prom_Enh.bed
cat $OUTPUT1 $OUTPUT2 | gsort -k 1,1 -k2,2n --parallel=8  -S 50% > $INPUT
rm $OUTPUT1 $OUTPUT2

OUTPUT1=$DATA_DIR/K_5UTR_intersected.bed
bedtools intersect -a $INPUT -b $UTR5 -wa  -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50% | awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"1}' >$OUTPUT1
OUTPUT2=$DATA_DIR/K_5UTR_non_intersected.bed
bedtools intersect -a $INPUT -b $UTR5 -v -sorted| gsort -k 1,1 -k2,2n --parallel=8  -S 50%| awk 'BEGIN {FS="\t"; OFS=","} {print $0"\t"0}' >$OUTPUT2

INPUT=$DATA_DIR/K_CGI_Prom_Enh_5UTR.bed
cat $OUTPUT1 $OUTPUT2 | gsort -k 1,1 -k2,2n --parallel=8  -S 50% > $INPUT
rm $OUTPUT1 $OUTPUT2

