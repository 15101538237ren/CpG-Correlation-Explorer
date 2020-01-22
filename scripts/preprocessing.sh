BASE_DIR=/Users/emmanueldollinger/PycharmProjects/CpG-Correlation-Explorer
DATA_DIR=$BASE_DIR/data
PEAK_DATA=/Users/emmanueldollinger/PycharmProjects/prediction_by_k/DATA/PLOT_PREPARATION/DATA/PEAK_DATA
cd $DATA_DIR

K=$DATA_DIR/K.bed
declare -a HISTONES=('H3k4me1' 'H3k4me2' 'H3k4me3' 'H3k9me3' 'H3k9ac' 'H3k27ac' 'H3k27me3' 'H3k36me3' 'H4k20me' 'CTCF' 'P300')
OUTPUT=$DATA_DIR/K_intersect_more.bed
INPUT=$DATA_DIR/K_intersect.bed

# Repeat change the following histone marks for more columns in OUTPUT
histone="P300"
histone_file=$PEAK_DATA/$histone".bed"
echo $histone_file
bedtools intersect -a $INPUT -b $histone_file -wa -wb -loj | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$20}'> $OUTPUT


K_EXTEND=$DATA_DIR/K-100bp-local-window.bed
RADIUS=50
awk -v radius=$RADIUS 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2-radius"\t"$3+radius-1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17}'  $INPUT>$K_EXTEND

K_WITH_CPG_DENSITY=$DATA_DIR/K-100bp-CpGd.bed
bedtools nuc -fi $DATA_DIR/hg19.fa -bed $K_EXTEND | awk -v radius=$RADIUS 'BEGIN {FS="\t"; OFS=","} {if(NR!=1) {print $1"\t"$2+radius"\t"$3-radius+1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$19}}'>$K_WITH_CPG_DENSITY

chromHMM=$DATA_DIR/HmmH1hescHMM.bed
OUTPUT=$DATA_DIR/K-100bp-CpGd_HMM.bed
bedtools intersect -a $K_WITH_CPG_DENSITY -b $chromHMM -wa -wb -loj | head -n 10 | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$22}'> $OUTPUT






















