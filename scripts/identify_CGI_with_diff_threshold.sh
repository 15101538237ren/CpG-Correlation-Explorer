BASE_DIR=/Users/emmanueldollinger/PycharmProjects/CpG-Correlation-Explorer
DATA_DIR=$BASE_DIR/data
CGI_OUTPUT=$DATA_DIR/CGI_identified_with_different_thereshold
mkdir -p $CGI_OUTPUT
cd $CGI_OUTPUT

cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 0.6 -noplot -outfile CGI_06_ratio.cpgplot -outfeat CGI_06_ratio.gff &
cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 0.3 -noplot -outfile CGI_03_ratio.cpgplot -outfeat CGI_03_ratio.gff &
cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 0.4 -noplot -outfile CGI_04_ratio.cpgplot -outfeat CGI_04_ratio.gff &
cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 0.5 -noplot -outfile CGI_05_ratio.cpgplot -outfeat CGI_05_ratio.gff &
cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 0.7 -noplot -outfile CGI_07_ratio.cpgplot -outfeat CGI_07_ratio.gff &
cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 0.8 -noplot -outfile CGI_08_ratio.cpgplot -outfeat CGI_08_ratio.gff &
cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 0.9 -noplot -outfile CGI_09_ratio.cpgplot -outfeat CGI_09_ratio.gff &
cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 1.2 -noplot -outfile CGI_12_ratio.cpgplot -outfeat CGI_12_ratio.gff &
cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 1.5 -noplot -outfile CGI_15_ratio.cpgplot -outfeat CGI_15_ratio.gff &
cpgplot -sequence ../hg19.fa -window 100 -minlen 200 -minpc 50 -minoe 1.8 -noplot -outfile CGI_18_ratio.cpgplot -outfeat CGI_18_ratio.gff &

K=$DATA_DIR/K.bed
declare -a CGI_RATIOS=("12" "15" "18") # "03" "04" "05" "06" "07" "08" "09" 
for cgr in "${CGI_RATIOS[@]}"; 
do
	INPUT=$CGI_OUTPUT/CGI_"$cgr"_ratio.gff    
	OUTPUT=$CGI_OUTPUT/CGI_"$cgr"_ratio.bed
	echo $INPUT
	awk '!/^ *#/ {print $1"\t"$4"\t"$5}' $INPUT | gsort -k 1,1 -k2,2n --parallel=8 > $OUTPUT
	INPUT=$OUTPUT
	OUTPUT=$CGI_OUTPUT/CGI_"$cgr"_K_intersected.bed
	bedtools intersect -a $K -b $INPUT -wa  -sorted | gsort -k 1,1 -k2,2n --parallel=8 | awk 'BEGIN {FS="\t"; OFS=","} {print;}' >$OUTPUT
done


# Second Trail
BASE_DIR=/Users/emmanueldollinger/PycharmProjects/CpG-Correlation-Explorer
DATA_DIR=$BASE_DIR/data
CGI_OUTPUT=$DATA_DIR/CGI_identified_with_different_thereshold
mkdir -p $CGI_OUTPUT
cd $CGI_OUTPUT
CGI_ORIGIN=$DATA_DIR/CGI.tsv
K=$DATA_DIR/K.bed
declare -a RATIOS=(0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3)
declare -a CGI_RATIOS=("06" "07" "08" "09" "10" "11" "12" "13")
for cgi in "${!CGI_RATIOS[@]}"; 
do
	cglabel="${CGI_RATIOS[$cgi]}"
	cgratio=${RATIOS[$cgi]}  
	OUTPUT=$CGI_OUTPUT/CGI_"$cglabel"_ratio.bed
	echo $cgratio
	awk -v ratio=$cgratio '!/^ *#/ && NR>1 { if($12>ratio) {print $2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$1}}' $CGI_ORIGIN | gsort -k 1,1 -k2,2n --parallel=8 > $OUTPUT
	INPUT=$OUTPUT
	OUTPUT=$CGI_OUTPUT/CGI_"$cglabel"_K_intersected.bed
	bedtools intersect -a $K -b $INPUT -wa -wb -sorted | gsort -k 1,1 -k2,2n --parallel=8>$OUTPUT
done







