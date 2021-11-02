#!/bin/sh
set -e -x

BAM="./../data/chromosom2_exons.bed"
INDEX="/home/udo/worms/"
NOW=$( date '+%F' )
OUTPUT="./../analysis/exon_quantification/"
ANNOTATION=""

#test worked:
#bathometer index G1-N2-.bam.idx /home/udo/worms/G1-N2-.bam 

###############make indexes for all bam files
#for i in /home/udo/worms/*.bam; do
#	echo ${i##*/}.idx
#	bathometer index ${i##*/}.idx $i
#ydone 
#worked!

################count exon coverage for all exons in bed file
mkdir -p $OUTPUT
bathometer fpkm ${INDEX}*.idx ::: $BAM > ${OUTPUT}exon_coverage_chromosom2_exons_FPKM.tsv
# only counts as output 


###############combine counts with annotation file 





