BASE_DIR=/Users/emmanueldollinger/PycharmProjects/CpG-Correlation-Explorer
DATA_DIR=$BASE_DIR/data
Genomic_Context=$DATA_DIR/Genomic_Context
cd $DATA_DIR

bigWigToWig DNase.bigWig DNase.wig
## convert wig to bed
wig2bed < DNase.wig > DNase.bed
gsort -k 1,1 -k2,2n --parallel=8 DNase.bed | awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$5}' > DNase.sorted.bed

K=$DATA_DIR/K.bed

DNASE=$DATA_DIR/DNase.bed
KDNASE=$DATA_DIR/K_DNase.bed
cat $DNASE| awk '/^chr[1-22|X|Y]/'> DNase.sorted.bed
DNASE=$DATA_DIR/DNase.sorted.bed
bedtools intersect -a $K -b $DNASE -wa -wb -loj -sorted| awk 'BEGIN {FS="\t"; OFS=","} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10}'> $KDNASE

paste K-100bp-CpGd_HMM.bed K_DNase.bed K_nucleo.bed | cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,26,33 > K-100bp-CpGd_HMM_DNase_Nucleo.bed


methy_and_k=$DATA_DIR/K-100bp-CpGd_HMM_DNase_Nucleo.bed
K_region_intersect=$DATA_DIR/K_region_intersect

mkdir -p $K_region_intersect
declare -a MARKERS=('Genomic_Regions' 'Histone_Modification' 'ChromHMM' 'TFBS') #
for region in "${MARKERS[@]}"; 
do
	WORK_DIR=$Genomic_Context/$region
	cd $WORK_DIR
	OUT_DIR=$K_region_intersect/$region
    mkdir -p $OUT_DIR
	#Do intersection
	for f in *.bed;
	do
		filename=${f%%.*}
		OUT_FILE=$OUT_DIR/$filename.bed
		echo $filename
		bedtools intersect -a $methy_and_k -b $f -wa -sorted>$OUT_FILE
	done
done


